package raj.cbm.alg;

import raj.cbm.SimulationConfig;

import java.util.*;
import java.util.stream.Collectors;

/**
 * SPEAII (Strength Pareto Evolutionary Algorithm II) — real implementation for CBM mapping
 * ----------------------------------------------------------------------------------------
 * Representation: integer array of length N (N = #cloudlets), gene i ∈ [0, VM-1] gives the VM index
 * selected for cloudlet i (job→VM mapping).
 *
 * Objectives (minimize): {energy (kWh), P95 latency (ms), SLO-miss (%)} evaluated using a
 * fast mapping-dependent surrogate (same family used by MOPSO/MOGA/MOACO/NSGAII in this project).
 *
 * SPEA2 pipeline (per Zitzler et al., 2001/2002):
 *  • Maintain a population P and an external archive A (bounded).
 *  • Compute fitness for the union U=P∪A using strength S, raw fitness R, density D (k-NN), F=R+D.
 *  • Environmental selection builds next archive A' of fixed size:
 *      - Insert all nondominated (R=0) solutions; if overflow, truncate by crowding (k-NN).
 *      - If underfull, fill with dominated with smallest F.
 *  • Mating selection: binary tournament biased to lower F; apply crossover & mutation to create new P.
 *  • Iterate for maxGen generations.
 *
 * Returned {@link Algorithm.Result} contains a deployable mapping (best archive member by
 * SLO→P95→Energy), the Pareto front, and simple HV/IGD indicators with a convergence iteration.
 *
 * Notes:
 *  • Self-contained; no extra files needed. To evaluate via CloudSim, replace evaluate(...)
 *    with SimulationRunner-based scoring and consider caching.
 *
 * @author rchoudhary
 */
public class SPEAII implements Algorithm {
    private final SimulationConfig sc;
    private final Random rnd;

    // Hyperparameters
    private int popSize;          // population size
    private int archiveSize = 120; // bounded external archive
    private int maxGen;           // generations
    private double pc = 0.9;      // crossover probability
    private double pm = 0.06;     // per-gene mutation probability

    public SPEAII(SimulationConfig sc, long seed) {
        this.sc = sc;
        this.rnd = new Random(seed ^ 1954196260373900364L);
        int n = Math.max(1, sc.totalCloudlets());
        this.popSize = Math.max(40, Math.min(120, 16 + n / 4));
        this.maxGen  = Math.max(30, Math.min(120, 20 + n / 8));
        if (archiveSize < popSize) archiveSize = popSize; // typical SPEA2 choice
    }

    @Override
    public Result optimize(String outPrefix) {
        final int N  = Math.max(1, sc.totalCloudlets());
        final int VM = Math.max(1, sc.VMS);
        final int[] cls = buildClassVector();
        final int[] dlMs = new int[]{ sc.DL_LOW_MS, sc.DL_MED_MS, sc.DL_HIGH_MS, sc.DL_CRITICAL_MS };

        // Initialize population and empty archive
        List<Solution> P = new ArrayList<>(popSize);
        for (int i = 0; i < popSize; i++) P.add(randomSol(N, VM, cls, dlMs));
        List<Solution> A = new ArrayList<>();

        int convergeIter = 0; double lastHV = 0.0;

        for (int gen = 1; gen <= maxGen; gen++) {
            // 1) Fitness assignment on union U = P ∪ A
            List<Solution> U = new ArrayList<>(P.size() + A.size());
            U.addAll(P); U.addAll(A);
            assignSpea2Fitness(U);

            // 2) Environmental selection → build next archive A'
            A = environmentalSelection(U, archiveSize);

            // 3) Mating selection + variation → new P
            P = reproduce(A, VM, cls, dlMs);

            // Convergence via HV over archive
            double hv = hypervolumeMC(A, 2000);
            if (hv > lastHV + 1e-6) { lastHV = hv; convergeIter = gen; }
        }

        // Final union to make sure archive is up to date
        List<Solution> U = new ArrayList<>(P.size() + A.size());
        U.addAll(P); U.addAll(A);
        assignSpea2Fitness(U);
        A = environmentalSelection(U, archiveSize);

        // Deploy mapping: prefer by SLO → P95 → Energy
        Solution deploy = A.stream()
                .sorted(Comparator
                        .comparingDouble((Solution s)->s.fit.slo)
                        .thenComparingDouble(s->s.fit.p95)
                        .thenComparingDouble(s->s.fit.energy))
                .findFirst()
                .orElse(P.get(0));

        Result r = new Result();
        r.mapping = deploy.x.clone();
        r.paretoFront = A.stream().map(s -> new double[]{s.fit.energy, s.fit.p95, s.fit.slo}).collect(Collectors.toList());
        r.hv = hypervolumeMC(A, 8000);
        r.igd = igdToIdeal(A);
        r.convergeIters = convergeIter;
        return r;
    }

    // ===== Solution type =====
    private static class Solution {
        int[] x;            // mapping (cloudlet -> VM)
        Fitness fit;        // objective values
        // SPEA2 fitness components
        double strength;    // S(i)
        double raw;         // R(i)
        double density;     // D(i) = 1/(sigma_k + 2)
        double f;           // final fitness F(i) = R(i) + D(i) (smaller is better)
        Solution(int[] x){ this.x = x; }
    }

    private Solution randomSol(int N, int VM, int[] cls, int[] dl){
        int[] x = new int[N];
        for (int i=0;i<N;i++) x[i] = rnd.nextInt(VM);
        Solution s = new Solution(x);
        evaluate(s, cls, dl);
        return s;
    }

    // ===== SPEA2 core steps =====

    /** Assign strength, raw, density and final fitness for all individuals in U. */
    private void assignSpea2Fitness(List<Solution> U){
        // 1) Strength S(i): number of solutions i dominates
        for (Solution s : U) s.strength = 0.0;
        for (int i=0;i<U.size();i++){
            Solution a = U.get(i);
            for (int j=0;j<U.size();j++) if (i!=j){
                Solution b = U.get(j);
                if (dominates(a.fit, b.fit)) a.strength += 1.0;
            }
        }
        // 2) Raw fitness R(i): sum of strengths of solutions dominating i
        for (Solution s : U) s.raw = 0.0;
        for (int i=0;i<U.size();i++){
            Solution a = U.get(i);
            for (int j=0;j<U.size();j++) if (i!=j){
                Solution b = U.get(j);
                if (dominates(b.fit, a.fit)) a.raw += b.strength;
            }
        }
        // 3) Density D(i): based on k-th nearest neighbor in objective space (normalized)
        int kNN = Math.max(1, (int)Math.round(Math.sqrt(U.size())));
        // Precompute normalized objective vectors
        double[] min = new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        double[] max = new double[]{-1, -1, -1};
        for (Solution s: U){
            min[0]=Math.min(min[0], s.fit.energy); max[0]=Math.max(max[0], s.fit.energy);
            min[1]=Math.min(min[1], s.fit.p95);    max[1]=Math.max(max[1], s.fit.p95);
            min[2]=Math.min(min[2], s.fit.slo);    max[2]=Math.max(max[2], s.fit.slo);
        }
        double[] span = new double[]{Math.max(1e-12, max[0]-min[0]), Math.max(1e-12, max[1]-min[1]), Math.max(1e-12, max[2]-min[2])};
        List<double[]> norm = new ArrayList<>(U.size());
        for (Solution s: U){
            norm.add(new double[]{ (s.fit.energy-min[0])/span[0], (s.fit.p95-min[1])/span[1], (s.fit.slo-min[2])/span[2] });
        }
        for (int i=0;i<U.size();i++){
            double[] ai = norm.get(i);
            double[] dists = new double[U.size()-1];
            int t=0;
            for (int j=0;j<U.size();j++) if (i!=j){
                double[] bj = norm.get(j);
                double de = ai[0]-bj[0], dp = ai[1]-bj[1], ds = ai[2]-bj[2];
                dists[t++] = Math.sqrt(Math.max(0, de*de + dp*dp + ds*ds));
            }
            Arrays.sort(dists);
            double sigmaK = dists[Math.min(kNN-1, dists.length-1)];
            U.get(i).density = 1.0 / (sigmaK + 2.0);
        }
        // 4) Final fitness
        for (Solution s : U) s.f = s.raw + s.density;
    }

    /** Environmental selection: build new archive A' of fixed size. */
    private List<Solution> environmentalSelection(List<Solution> U, int cap){
        // Split nondominated (raw==0) and dominated
        List<Solution> nonDom = U.stream().filter(s -> s.raw == 0.0).collect(Collectors.toList());
        // De-duplicate identical objective points (optional)
        // Fill archive
        List<Solution> A = new ArrayList<>(cap);
        if (nonDom.size() >= cap){
            // Truncate by k-NN crowding until size=cap
            A.addAll(nonDom);
            while (A.size() > cap){ truncateOnce(A); }
        } else {
            A.addAll(nonDom);
            // Add dominated by increasing final fitness f
            List<Solution> dominated = U.stream().filter(s -> s.raw > 0.0).collect(Collectors.toList());
            dominated.sort(Comparator.comparingDouble(s -> s.f));
            for (Solution s: dominated){ if (A.size() < cap) A.add(s); else break; }
        }
        if (A.isEmpty()){
            // Fallback: pick the single best by f
            U.sort(Comparator.comparingDouble(s -> s.f));
            A.add(U.get(0));
        }
        return A;
    }

    /** Remove one individual from A using SPEA2 truncation (closest pair / nearest neighbor). */
    private void truncateOnce(List<Solution> A){
        int m = A.size();
        if (m <= 1) return;
        // Build normalized objective matrix
        double minE=Double.POSITIVE_INFINITY, minP=Double.POSITIVE_INFINITY, minS=Double.POSITIVE_INFINITY;
        double maxE=-1, maxP=-1, maxS=-1;
        for (Solution s: A){
            minE=Math.min(minE,s.fit.energy); maxE=Math.max(maxE,s.fit.energy);
            minP=Math.min(minP,s.fit.p95);    maxP=Math.max(maxP,s.fit.p95);
            minS=Math.min(minS,s.fit.slo);    maxS=Math.max(maxS,s.fit.slo);
        }
        double spanE=Math.max(1e-12, maxE-minE), spanP=Math.max(1e-12, maxP-minP), spanS=Math.max(1e-12, maxS-minS);
        double[][] norm = new double[m][3];
        for (int i=0;i<m;i++){
            Solution s = A.get(i);
            norm[i][0]=(s.fit.energy-minE)/spanE; norm[i][1]=(s.fit.p95-minP)/spanP; norm[i][2]=(s.fit.slo-minS)/spanS;
        }
        // Compute all pairwise distances and find the minimal nearest-neighbor distance
        double bestDist = Double.POSITIVE_INFINITY; int removeIdx = -1;
        for (int i=0;i<m;i++){
            double nn = Double.POSITIVE_INFINITY;
            for (int j=0;j<m;j++) if (i!=j){
                double de = norm[i][0]-norm[j][0], dp = norm[i][1]-norm[j][1], ds = norm[i][2]-norm[j][2];
                double d = Math.sqrt(Math.max(0, de*de + dp*dp + ds*ds));
                if (d < nn) nn = d;
            }
            if (nn < bestDist){ bestDist = nn; removeIdx = i; }
        }
        if (removeIdx >= 0) A.remove(removeIdx);
    }

    /** Build next generation via tournament selection on A, then crossover/mutation. */
    private List<Solution> reproduce(List<Solution> A, int VM, int[] cls, int[] dl){
        List<Solution> next = new ArrayList<>(popSize);
        if (A.isEmpty()) return next;
        // Tournament favors smaller F; if tie, random
        while (next.size() < popSize){
            Solution p1 = tournamentByF(A);
            Solution p2 = tournamentByF(A);
            Solution c1 = cloneOf(p1);
            Solution c2 = cloneOf(p2);
            if (rnd.nextDouble() < pc) uniformCrossover(c1, c2);
            mutate(c1, VM); mutate(c2, VM);
            evaluate(c1, cls, dl); evaluate(c2, cls, dl);
            next.add(c1); if (next.size() < popSize) next.add(c2);
        }
        return next;
    }

    private Solution tournamentByF(List<Solution> pool){
        Solution a = pool.get(rnd.nextInt(pool.size()));
        Solution b = pool.get(rnd.nextInt(pool.size()));
        if (a.f < b.f) return a; if (b.f < a.f) return b; return rnd.nextBoolean()?a:b;
    }

    private static Solution cloneOf(Solution s){ Solution c = new Solution(s.x.clone()); c.fit = s.fit; c.strength=s.strength; c.raw=s.raw; c.density=s.density; c.f=s.f; return c; }

    // ===== Variation operators =====
    private void uniformCrossover(Solution a, Solution b){
        for (int i=0;i<a.x.length;i++) if (rnd.nextBoolean()){ int t=a.x[i]; a.x[i]=b.x[i]; b.x[i]=t; }
    }
    private void mutate(Solution s, int VM){
        for (int i=0;i<s.x.length;i++) if (rnd.nextDouble() < pm) s.x[i] = rnd.nextInt(VM);
    }

    // ===== Fitness surrogate (shared family) =====
    private void evaluate(Solution s, int[] cls, int[] dlMs){ s.fit = evaluate(s.x, cls, dlMs); }

    private int[] buildClassVector(){
        final int N = Math.max(1, sc.totalCloudlets());
        int[] cls = new int[N]; int idx=0;
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
        for (int i=0;i<N;i++){
            int vm = Math.floorMod(mapping[i], VM);
            int k = counts[vm]++;
            double finishMs = (k+1) * jobServiceSec * 1000.0; // sequential per VM
            compMs[i] = lastFinish[vm] = finishMs;
        }
        double[] sorted = compMs.clone(); Arrays.sort(sorted);
        double p95 = sorted[(int)Math.floor(0.95 * (sorted.length - 1))];

        int viol=0; for (int i=0;i<N;i++) if (compMs[i] > dlMs[cls[i]]) viol++;
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

    // ===== Dominance helpers =====
    private static boolean dominates(Fitness a, Fitness b){
        boolean be = a.energy <= b.energy && a.p95 <= b.p95 && a.slo <= b.slo;
        boolean sb = a.energy < b.energy || a.p95 < b.p95 || a.slo < b.slo;
        return be && sb;
    }
    private static boolean equalsFit(Fitness a, Fitness b){
        return Math.abs(a.energy-b.energy) < 1e-12 && Math.abs(a.p95-b.p95) < 1e-12 && Math.abs(a.slo-b.slo) < 1e-12;
    }

    private static int fill(int[] a,int start,int count,int val){ for(int i=0;i<count && start+i<a.length;i++) a[start+i]=val; return Math.min(a.length, start+count); }

    // ===== Indicators (simple approximations) =====
    private static double hypervolumeMC(List<Solution> front, int samples){
        if (front.isEmpty()) return 0.0;
        double minE=1e9,minP=1e9,minS=1e9,maxE=0,maxP=0,maxS=0;
        for (Solution s : front){ minE=Math.min(minE,s.fit.energy); maxE=Math.max(maxE,s.fit.energy);
                                  minP=Math.min(minP,s.fit.p95);    maxP=Math.max(maxP,s.fit.p95);
                                  minS=Math.min(minS,s.fit.slo);    maxS=Math.max(maxS,s.fit.slo); }
        double refE=maxE*1.2+1e-6, refP=maxP*1.2+1e-6, refS=maxS*1.2+1e-6;
        double vol=(refE-minE)*(refP-minP)*(refS-minS); if (vol<=0) return 0.0; int hit=0; Random r = new Random(2112);
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
