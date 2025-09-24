package raj.cbm.alg;

import raj.cbm.SimulationConfig;

import java.util.*;
import java.util.stream.Collectors;

/**
 * NSGAII (Non-dominated Sorting Genetic Algorithm II) — real implementation for CBM mapping
 * ----------------------------------------------------------------------------------------
 * Representation: integer array of length N (N = #cloudlets), where gene i ∈ [0, VM-1]
 * is the VM index chosen for cloudlet i.
 *
 * Objectives (minimize): {energy (kWh), P95 latency (ms), SLO-miss (%)} evaluated using a
 * fast mapping-dependent surrogate (same family used by MOPSO/MOGA/MOACO in this project).
 *
 * NSGA-II pipeline:
 *  • Initialization: random population of mappings
 *  • Fast non-dominated sorting → fronts F1, F2, ...
 *  • Crowding distance within fronts
 *  • Parent selection: binary tournament on (rank, crowding)
 *  • Variation: uniform crossover + per-gene random-reset mutation
 *  • Environmental selection: combine parents+offspring → keep best by fronts + crowding
 *
 * Returned {@link Algorithm.Result} contains a deployable mapping (from the archive's
 * best-by SLO→P95→energy), the Pareto front, and simple HV/IGD indicators.
 *
 * Notes:
 *  • This class is self-contained; no extra files required. If you later want
 *    ground-truth CloudSim evaluation for candidates, replace evaluate(...) body
 *    with a call to SimulationRunner and consider caching.
 *
 * @author rchoudhary
 */
public class NSGAII implements Algorithm {
    private final SimulationConfig sc;
    private final Random rnd;

    // Hyperparameters (tweak as needed)
    private int popSize;
    private int maxGen;
    private double pc = 0.9;    // crossover probability
    private double pm = 0.06;   // mutation probability per gene

    public NSGAII(SimulationConfig sc, long seed) {
        this.sc = sc;
        this.rnd = new Random(seed ^ 7737619851417673738L);
        int n = Math.max(1, sc.totalCloudlets());
        this.popSize = Math.max(40, Math.min(120, 16 + n / 4));
        this.maxGen  = Math.max(30, Math.min(120, 20 + n / 8));
    }

    @Override
    public Result optimize(String outPrefix) {
        final int N  = Math.max(1, sc.totalCloudlets());
        final int VM = Math.max(1, sc.VMS);
        final int[] cls = buildClassVector();
        final int[] dlMs = new int[]{ sc.DL_LOW_MS, sc.DL_MED_MS, sc.DL_HIGH_MS, sc.DL_CRITICAL_MS };

        // Initialize population
        List<Solution> pop = new ArrayList<>(popSize);
        for (int i = 0; i < popSize; i++) pop.add(randomSol(N, VM, cls, dlMs));

        // Archive holds a rolling non-dominated set (for convenience/return)
        List<Solution> archive = new ArrayList<>();
        updateArchive(archive, pop);

        int convergeIter = 0; double lastHV = hypervolumeMC(archive, 2000);

        for (int gen = 1; gen <= maxGen; gen++) {
            // Parent selection
            List<Solution> parents = new ArrayList<>(pop);
            // Variation: create offspring via tournament + crossover + mutation
            List<Solution> offspring = new ArrayList<>(popSize);
            while (offspring.size() < popSize) {
                Solution p1 = tournament(parents);
                Solution p2 = tournament(parents);
                Solution c1 = cloneOf(p1);
                Solution c2 = cloneOf(p2);
                if (rnd.nextDouble() < pc) uniformCrossover(c1, c2);
                mutate(c1, VM); mutate(c2, VM);
                evaluate(c1, cls, dlMs); evaluate(c2, cls, dlMs);
                offspring.add(c1); offspring.add(c2);
            }
            if (offspring.size() > popSize) offspring = offspring.subList(0, popSize);

            // Combine and perform environmental selection via non-dom sorting + crowding
            List<Solution> union = new ArrayList<>(pop.size() + offspring.size());
            union.addAll(pop); union.addAll(offspring);
            pop = selectByNondomAndCrowding(union, popSize);

            // Update archive and convergence indicator
            updateArchive(archive, pop);
            double hv = hypervolumeMC(archive, 2000);
            if (hv > lastHV + 1e-6) { lastHV = hv; convergeIter = gen; }
        }

        // Choose a deploy mapping from the archive with practical tie-breaks
        Solution deploy = archive.stream()
                .sorted(Comparator
                        .comparingDouble((Solution s) -> s.fit.slo)
                        .thenComparingDouble(s -> s.fit.p95)
                        .thenComparingDouble(s -> s.fit.energy))
                .findFirst()
                .orElse(pop.get(0));

        Result r = new Result();
        r.mapping = deploy.x.clone();
        r.paretoFront = archive.stream()
                .map(s -> new double[]{ s.fit.energy, s.fit.p95, s.fit.slo })
                .collect(Collectors.toList());
        r.hv = hypervolumeMC(archive, 8000);
        r.igd = igdToIdeal(archive);
        r.convergeIters = convergeIter;
        return r;
    }

    // ===== Solution type =====
    private static class Solution {
        int[] x;        // mapping (cloudlet -> VM)
        Fitness fit;    // objective values
        int rank;       // non-dominated rank (1 = best/front 1)
        double crowd;   // crowding distance
        Solution(int[] x){ this.x = x; }
    }

    private Solution randomSol(int N, int VM, int[] cls, int[] dl){
        int[] x = new int[N];
        for (int i = 0; i < N; i++) x[i] = rnd.nextInt(VM);
        Solution s = new Solution(x);
        evaluate(s, cls, dl);
        return s;
    }

    private static Solution cloneOf(Solution s){ Solution c = new Solution(s.x.clone()); c.fit = s.fit; c.rank = s.rank; c.crowd = s.crowd; return c; }

    // ===== Tournament selection on (rank, crowding) =====
    private Solution tournament(List<Solution> pool){
        Solution a = pool.get(rnd.nextInt(pool.size()));
        Solution b = pool.get(rnd.nextInt(pool.size()));
        if (dominatesByRankCrowding(a, b)) return a;
        if (dominatesByRankCrowding(b, a)) return b;
        return rnd.nextBoolean() ? a : b;
    }

    private static boolean dominatesByRankCrowding(Solution a, Solution b){
        if (a.rank < b.rank) return true;
        if (a.rank > b.rank) return false;
        return a.crowd > b.crowd; // larger crowding preferred
    }

    // ===== Variation operators =====
    private void uniformCrossover(Solution a, Solution b){
        for (int i = 0; i < a.x.length; i++) if (rnd.nextBoolean()) { int t = a.x[i]; a.x[i] = b.x[i]; b.x[i] = t; }
    }

    private void mutate(Solution s, int VM){
        for (int i = 0; i < s.x.length; i++) if (rnd.nextDouble() < pm) s.x[i] = rnd.nextInt(VM);
    }

    // ===== Environmental selection (NSGA-II) =====
    private List<Solution> selectByNondomAndCrowding(List<Solution> union, int k){
        List<List<Solution>> fronts = fastNonDominatedSort(union);
        List<Solution> next = new ArrayList<>(k);
        for (List<Solution> F : fronts){
            if (next.size() + F.size() <= k) {
                computeCrowding(F);
                next.addAll(F);
            } else {
                computeCrowding(F);
                F.sort((s1, s2) -> Double.compare(s2.crowd, s1.crowd)); // desc by crowding
                for (Solution s : F){ if (next.size() < k) next.add(s); }
                break;
            }
        }
        return next;
    }

    private List<List<Solution>> fastNonDominatedSort(List<Solution> P){
        // Simple O(MN^2) implementation (M=3 objectives)
        List<List<Solution>> fronts = new ArrayList<>();
        Map<Solution, List<Solution>> S = new HashMap<>();
        Map<Solution, Integer> n = new HashMap<>();
        List<Solution> F1 = new ArrayList<>();
        for (Solution p : P){ S.put(p, new ArrayList<>()); n.put(p, 0); }
        for (int i = 0; i < P.size(); i++){
            Solution p = P.get(i);
            for (int j = 0; j < P.size(); j++) if (i != j){
                Solution q = P.get(j);
                if (dominates(p.fit, q.fit)) S.get(p).add(q);
                else if (dominates(q.fit, p.fit)) n.put(p, n.get(p)+1);
            }
            if (n.get(p) == 0){ p.rank = 1; F1.add(p); }
        }
        fronts.add(F1);
        int k = 0;
        while (k < fronts.size()){
            List<Solution> F = fronts.get(k);
            List<Solution> H = new ArrayList<>();
            for (Solution p : F){
                for (Solution q : S.get(p)){
                    n.put(q, n.get(q)-1);
                    if (n.get(q) == 0){ q.rank = k + 2; H.add(q); }
                }
            }
            if (!H.isEmpty()) fronts.add(H);
            k++;
        }
        return fronts;
    }

    private void computeCrowding(List<Solution> F){
        if (F.isEmpty()) return;
        for (Solution s : F) s.crowd = 0.0;
        crowdingFor(F, Comparator.comparingDouble(s -> s.fit.energy), s -> s.fit.energy);
        crowdingFor(F, Comparator.comparingDouble(s -> s.fit.p95),    s -> s.fit.p95);
        crowdingFor(F, Comparator.comparingDouble(s -> s.fit.slo),    s -> s.fit.slo);
    }
    private interface Obj { double v(Solution s); }
    private void crowdingFor(List<Solution> F, Comparator<Solution> cmp, Obj f){
        F.sort(cmp);
        F.get(0).crowd = F.get(F.size()-1).crowd = Double.POSITIVE_INFINITY;
        double min = f.v(F.get(0)), max = f.v(F.get(F.size()-1));
        double range = Math.max(1e-12, max - min);
        for (int i = 1; i < F.size()-1; i++) {
            F.get(i).crowd += (f.v(F.get(i+1)) - f.v(F.get(i-1))) / range;
        }
    }

    // ===== Archive maintenance =====
    private void updateArchive(List<Solution> archive, List<Solution> cand){
        for (Solution s : cand) addToArchive(archive, s);
        // Optional cap
        if (archive.size() > 200) archive.sort((a,b) -> a.rank!=b.rank ? Integer.compare(a.rank,b.rank) : Double.compare(b.crowd, a.crowd));
        while (archive.size() > 200) archive.remove(archive.size()-1);
    }

    private void addToArchive(List<Solution> arch, Solution s){
        arch.removeIf(o -> dominates(s.fit, o.fit) && !equalsFit(s.fit, o.fit));
        for (Solution o : arch) if (dominates(o.fit, s.fit) && !equalsFit(o.fit, s.fit)) return;
        arch.add(s);
    }

    // ===== Fitness & helpers (surrogate; same family as other algos) =====
    private void evaluate(Solution s, int[] cls, int[] dlMs){ s.fit = evaluate(s.x, cls, dlMs); }

    private int[] buildClassVector(){
        final int N = Math.max(1, sc.totalCloudlets());
        int[] cls = new int[N]; int idx = 0;
        idx = fill(cls, idx, sc.A_COUNT, 3);
        idx = fill(cls, idx, sc.B_COUNT, 2);
        idx = fill(cls, idx, sc.C_COUNT, 1);
        fill(cls, idx, sc.D_COUNT, 0);
        return cls;
    }

    private Fitness evaluate(int[] mapping, int[] cls, int[] dlMs){
        final int N = mapping.length;
        final int VM = Math.max(1, sc.VMS);
        final double jobServiceSec = sc.CLOUDLET_LENGTH / (double)(Math.max(1, sc.VM_PES) * Math.max(1, sc.HOST_MIPS));

        int[] counts = new int[VM];
        double[] lastFinish = new double[VM];
        double[] compMs = new double[N];
        for (int i = 0; i < N; i++){
            int vm = Math.floorMod(mapping[i], VM);
            int k = counts[vm]++;
            double finishMs = (k + 1) * jobServiceSec * 1000.0; // sequential per VM
            compMs[i] = lastFinish[vm] = finishMs;
        }
        double[] sorted = compMs.clone(); Arrays.sort(sorted);
        double p95 = sorted[(int)Math.floor(0.95 * (sorted.length - 1))];

        int viol = 0; for (int i = 0; i < N; i++) if (compMs[i] > dlMs[cls[i]]) viol++;
        double sloPct = 100.0 * viol / Math.max(1, N);

        int HOSTS = Math.max(1, sc.HOSTS), HOST_PES = Math.max(1, sc.HOST_PES);
        double makespanSec = Arrays.stream(lastFinish).max().orElse(0.0) / 1000.0;
        double totalCoreSeconds = N * (sc.CLOUDLET_LENGTH / (double)Math.max(1, sc.HOST_MIPS));
        double avgUtil = Math.min(1.0, totalCoreSeconds / (HOSTS * HOST_PES * Math.max(1e-9, makespanSec)));
        double pIdle = sc.STATIC_POWER_W, pMax = sc.MAX_POWER_W;
        double pMeanPerHost = pIdle + (pMax - pIdle) * avgUtil;
        double energyKwh = (HOSTS * pMeanPerHost * makespanSec) / 3600.0 / 1000.0;

        return new Fitness(energyKwh, p95, sloPct);
    }

    private static class Fitness { final double energy, p95, slo; Fitness(double e,double p,double s){ energy=e; p95=p; slo=s; } }

    private static boolean dominates(Fitness a, Fitness b){
        boolean be = a.energy <= b.energy && a.p95 <= b.p95 && a.slo <= b.slo;
        boolean sb = a.energy < b.energy || a.p95 < b.p95 || a.slo < b.slo;
        return be && sb;
    }
    private static boolean equalsFit(Fitness a, Fitness b){
        return Math.abs(a.energy-b.energy) < 1e-12 && Math.abs(a.p95-b.p95) < 1e-12 && Math.abs(a.slo-b.slo) < 1e-12;
    }

    private static int fill(int[] a, int start, int count, int val){ for(int i=0;i<count && start+i<a.length;i++) a[start+i]=val; return Math.min(a.length, start+count); }

    // ===== Indicators (simple approximations) =====
    private static double hypervolumeMC(List<Solution> front, int samples){
        if (front.isEmpty()) return 0.0;
        double minE=1e9,minP=1e9,minS=1e9,maxE=0,maxP=0,maxS=0;
        for (Solution s : front){ minE=Math.min(minE,s.fit.energy); maxE=Math.max(maxE,s.fit.energy);
                                  minP=Math.min(minP,s.fit.p95);    maxP=Math.max(maxP,s.fit.p95);
                                  minS=Math.min(minS,s.fit.slo);    maxS=Math.max(maxS,s.fit.slo); }
        double refE=maxE*1.2+1e-6, refP=maxP*1.2+1e-6, refS=maxS*1.2+1e-6;
        double vol=(refE-minE)*(refP-minP)*(refS-minS); if (vol<=0) return 0.0; int hit=0; Random r = new Random(2027);
        for (int i=0;i<samples;i++){
            double e=minE+r.nextDouble()*(refE-minE);
            double p=minP+r.nextDouble()*(refP-minP);
            double s=minS+r.nextDouble()*(refS-minS);
            if (isDominatedByFront(e,p,s, front)) hit++;
        }
        return (hit/(double)samples) * vol;
    }

    private static boolean isDominatedByFront(double e,double p,double s, List<Solution> front){
        for (Solution x: front){
            if (x.fit.energy<=e && x.fit.p95<=p && x.fit.slo<=s && (x.fit.energy<e || x.fit.p95<p || x.fit.slo<s)) return true;
        }
        return false;
    }

    private static double igdToIdeal(List<Solution> front){
        if (front.isEmpty()) return 0.0;
        double minE=1e9,minP=1e9,minS=1e9; for (Solution s: front){ minE=Math.min(minE,s.fit.energy); minP=Math.min(minP,s.fit.p95); minS=Math.min(minS,s.fit.slo);} double acc=0;
        for (Solution s: front){ double de=s.fit.energy-minE, dp=s.fit.p95-minP, ds=s.fit.slo-minS; acc += Math.sqrt(Math.max(0, de*de + dp*dp + ds*ds)); }
        return acc / front.size();
    }
}
