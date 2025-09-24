package raj.cbm.alg;

import raj.cbm.SimulationConfig;

import java.util.*;
import java.util.stream.Collectors;

/**
 * MOACO — Multi-Objective Ant Colony Optimization (real implementation)
 * <p>
 * Builds job→VM mappings for Cloud-CBM using a pheromone/heuristic construction and a
 * non-dominated external archive. Objectives (all minimized):
 * <ul>
 *   <li>Energy (kWh)</li>
 *   <li>P95 latency (ms)</li>
 *   <li>SLO-miss % (unweighted by class)</li>
 * </ul>
 * Representation: int[] mapping of length N (#cloudlets); gene i ∈ [0, VM-1].
 * <p>
 * Pipeline per iteration:
 * <ol>
 *   <li>Each ant constructs a mapping probabilistically with probability
 *       P(vm=j|task=i) ∝ τ[i][j]^α · η[i][j]^β, where η favors balancing and soft-deadline respect.</li>
 *   <li>Evaluate objectives via a fast surrogate (sequential per-VM service).</li>
 *   <li>Insert candidates into a bounded non-dominated archive.</li>
 *   <li>Pheromone evaporation and deposit from elite archive members.</li>
 * </ol>
 * <p>
 * Notes:
 * <ul>
 *   <li>This class is self-contained. To use CloudSim for ground truth, replace {@link #evaluate(int[], int[], int[])}.</li>
 *   <li>Deploy mapping is chosen by SLO → P95 → Energy to prioritize deadline satisfaction.</li>
 * </ul>
 *
 * @author rchoudhary
 */
public class MOACO implements Algorithm {
    private final SimulationConfig sc;
    private final Random rnd;

    // ===== Hyperparameters =====
    private int numAnts;          // ants per iteration
    private int maxIter;          // iterations
    private int archiveSize = 150;
    private double alpha = 1.2;   // pheromone exponent
    private double beta  = 3.0;   // heuristic exponent (favor heuristic strongly)
    private double rho   = 0.15;  // evaporation rate
    private double q     = 0.8;   // deposit magnitude scaler
    private int eliteK   = 24;    // #archive elites to deposit
    private double gammaSlack = 2.0; // deadline soft-penalty in heuristic
    private double mutateProb = 0.03; // small random reassignment per gene

    public MOACO(SimulationConfig sc, long seed) {
        this.sc = sc;
        this.rnd = new Random(seed ^ 3323661488414983760L);
        int n = Math.max(1, sc.totalCloudlets());
        this.numAnts = Math.max(24, Math.min(80, 12 + n / 6));
        this.maxIter = Math.max(30, Math.min(120, 20 + n / 8));
    }

    @Override
    public Result optimize(String outPrefix) {
        final int N  = Math.max(1, sc.totalCloudlets());
        final int VM = Math.max(1, sc.VMS);

        // Build class vector in the order A, B, C, D as produced by CloudResourceFactory;
        // we encode as 3=A, 2=B, 1=C, 0=D and map deadlines accordingly.
        final int[] cls = new int[N];
        int idx = 0;
        idx = fill(cls, idx, sc.A_COUNT, 3);
        idx = fill(cls, idx, sc.B_COUNT, 2);
        idx = fill(cls, idx, sc.C_COUNT, 1);
        fill(cls, idx, sc.D_COUNT, 0);

        final int[] dlMs = new int[]{ sc.DL_LOW_MS, sc.DL_MED_MS, sc.DL_HIGH_MS, sc.DL_CRITICAL_MS };

        // Pheromone matrix per task×VM
        double[][] tau = new double[N][VM];
        double tau0 = 1.0; // initial uniform pheromone
        for (int i = 0; i < N; i++) Arrays.fill(tau[i], tau0);

        // External archive
        Archive archive = new Archive(archiveSize);

        int convergeIter = 0; double lastHV = 0.0;

        for (int iter = 1; iter <= maxIter; iter++) {
            List<Solution> batch = new ArrayList<>(numAnts);

            for (int a = 0; a < numAnts; a++) {
                int[] mapping = constructMapping(tau, VM, cls, dlMs);
                // Light mutation to avoid stagnation
                mutateMapping(mapping, VM);
                Fitness f = evaluate(mapping, cls, dlMs);
                batch.add(new Solution(mapping, f));
            }

            // Push batch to archive
            for (Solution s : batch) archive.add(s);

            // Pheromone update
            evaporate(tau, rho);
            depositFromArchive(tau, archive, cls, dlMs);

            // Convergence heuristic via HV
            double hv = hypervolumeMC(archive.front, 2000);
            if (hv > lastHV + 1e-6) { lastHV = hv; convergeIter = iter; }
        }

        // Choose deploy mapping SLO→P95→Energy
        Solution deploy = archive.front.stream()
                .sorted(Comparator
                        .comparingDouble((Solution s)->s.fit.slo)
                        .thenComparingDouble(s->s.fit.p95)
                        .thenComparingDouble(s->s.fit.energy))
                .findFirst()
                .orElseGet(() -> new Solution(randomMapping(N, VM), evaluate(randomMapping(N, VM), cls, dlMs)));

        Result r = new Result();
        r.mapping = deploy.mapping.clone();
        r.paretoFront = archive.front.stream()
                .map(s -> new double[]{ s.fit.energy, s.fit.p95, s.fit.slo })
                .collect(Collectors.toList());
        r.hv = hypervolumeMC(archive.front, 8000);
        r.igd = igdToIdeal(archive.front);
        r.convergeIters = convergeIter;
        return r;
    }

    // ===== Construction =====

    /** Probabilistic mapping construction using τ^α · η^β (with online load-aware heuristic). */
    private int[] constructMapping(double[][] tau, int VM, int[] cls, int[] dlMs){
        final int N = tau.length;
        int[] map = new int[N];
        int[] load = new int[VM];
        double serviceMs = sc.CLOUDLET_LENGTH * 1000.0 / (double)(Math.max(1, sc.VM_PES) * Math.max(1, sc.HOST_MIPS));
        for (int i = 0; i < N; i++){
            // Build η for this task based on current loads and soft-deadline slack
            double[] eta = new double[VM];
            int c = cls[i];
            double dl = dlMs[c];
            for (int v = 0; v < VM; v++){
                double finish = (load[v] + 1) * serviceMs; // sequential proxy
                double slack = Math.max(0.0, finish - dl);
                // Heuristic: prefer low load and low deadline violation
                eta[v] = 1.0 / (1.0 + load[v]) * Math.exp(-gammaSlack * (slack / Math.max(1.0, dl)));
            }
            // Roulette wheel with τ^α · η^β
            int vm = sampleVM(tau[i], eta);
            map[i] = vm;
            load[vm]++;
        }
        return map;
    }

    private int sampleVM(double[] tau_i, double[] eta_i){
        double sum = 0.0; double[] w = new double[tau_i.length];
        for (int v = 0; v < tau_i.length; v++){
            double val = Math.pow(Math.max(1e-12, tau_i[v]), alpha) * Math.pow(Math.max(1e-12, eta_i[v]), beta);
            w[v] = val; sum += val;
        }
        if (sum <= 0) return rnd.nextInt(tau_i.length);
        double r = rnd.nextDouble() * sum, acc = 0.0;
        for (int v = 0; v < w.length; v++){ acc += w[v]; if (acc >= r) return v; }
        return w.length - 1;
    }

    private void evaporate(double[][] tau, double rho){
        for (double[] row : tau) for (int j = 0; j < row.length; j++) row[j] *= (1.0 - rho);
    }

    /** Deposit pheromones from up to {@code eliteK} archive members (quality-weighted). */
    private void depositFromArchive(double[][] tau, Archive archive, int[] cls, int[] dlMs){
        if (archive.front.isEmpty()) return;
        // Normalize objectives among archive for quality scoring
        double minE=1e9,minP=1e9,minS=1e9,maxE=0,maxP=0,maxS=0;
        for (Solution s : archive.front){
            minE=Math.min(minE,s.fit.energy); maxE=Math.max(maxE,s.fit.energy);
            minP=Math.min(minP,s.fit.p95);    maxP=Math.max(maxP,s.fit.p95);
            minS=Math.min(minS,s.fit.slo);    maxS=Math.max(maxS,s.fit.slo);
        }
        double spanE=Math.max(1e-9, maxE-minE), spanP=Math.max(1e-9, maxP-minP), spanS=Math.max(1e-9, maxS-minS);
        // Pick top-K by SLO→P95→Energy (deadline-first)
        List<Solution> elites = new ArrayList<>(archive.front);
        elites.sort(Comparator
                .comparingDouble((Solution s)->s.fit.slo)
                .thenComparingDouble(s->s.fit.p95)
                .thenComparingDouble(s->s.fit.energy));
        if (elites.size() > eliteK) elites = elites.subList(0, eliteK);
        for (Solution s : elites){
            double eN = (s.fit.energy - minE) / spanE;
            double pN = (s.fit.p95   - minP) / spanP;
            double sN = (s.fit.slo   - minS) / spanS;
            double quality = 1.0 / (1e-9 + 0.5*eN + 0.8*pN + 1.3*sN); // weight deadlines more
            // Deposit onto edges (task i, vm = mapping[i])
            for (int i = 0; i < s.mapping.length; i++){
                int v = s.mapping[i];
                tau[i][v] += q * quality;
            }
        }
        // Optional: cap pheromones to avoid explosion
        double cap = 10.0;
        for (double[] row : tau) for (int j = 0; j < row.length; j++) if (row[j] > cap) row[j] = cap;
    }

    private void mutateMapping(int[] map, int VM){
        for (int i = 0; i < map.length; i++) if (rnd.nextDouble() < mutateProb) map[i] = rnd.nextInt(VM);
    }

    private int[] randomMapping(int N, int VM){ int[] m = new int[N]; for(int i=0;i<N;i++) m[i]=rnd.nextInt(VM); return m; }

    // ===== Fitness surrogate (same family as other algos, but unweighted SLO) =====

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
            double finishMs = (k + 1) * jobServiceSec * 1000.0; // sequential per VM
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

    // ===== Archive (non-dominated set with crowding) =====

    private static class Solution {
        final int[] mapping; final Fitness fit; double crowd;
        Solution(int[] m, Fitness f){ this.mapping=m; this.fit=f; }
    }

    private static class Archive {
        final int cap; final List<Solution> front = new ArrayList<>();
        Archive(int cap){ this.cap = cap; }
        void add(Solution s){
            front.removeIf(o -> dominates(s.fit, o.fit) && !equalsFit(s.fit, o.fit));
            for (Solution o: front) if (dominates(o.fit, s.fit) && !equalsFit(o.fit, s.fit)) return;
            front.add(s);
            if (front.size() > cap) pruneByCrowding();
        }
        private void pruneByCrowding(){ computeCrowding(); Solution worst = front.stream().min(Comparator.comparingDouble(x->x.crowd)).orElse(null); if (worst!=null) front.remove(worst); }
        private void computeCrowding(){ if (front.size()<3){ for(Solution s: front) s.crowd=Double.POSITIVE_INFINITY; return; }
            crowd(Objective.ENERGY); crowd(Objective.P95); crowd(Objective.SLO); }
        enum Objective{ENERGY,P95,SLO}
        private void crowd(Objective o){ Comparator<Solution> cmp = switch(o){
                case ENERGY -> Comparator.comparingDouble(s->s.fit.energy);
                case P95    -> Comparator.comparingDouble(s->s.fit.p95);
                default     -> Comparator.comparingDouble(s->s.fit.slo);
            }; front.sort(cmp);
            front.get(0).crowd=front.get(front.size()-1).crowd=Double.POSITIVE_INFINITY;
            double min, max; switch (o){
                case ENERGY -> {min = front.get(0).fit.energy; max=front.get(front.size()-1).fit.energy;}
                case P95    -> {min = front.get(0).fit.p95;    max=front.get(front.size()-1).fit.p95;}
                default     -> {min = front.get(0).fit.slo;    max=front.get(front.size()-1).fit.slo;}
            }
            double range = Math.max(1e-12, max-min);
            for (int i=1;i<front.size()-1;i++){
                Solution prev=front.get(i-1), next=front.get(i+1), cur=front.get(i);
                double delta = switch (o){
                    case ENERGY -> (next.fit.energy - prev.fit.energy)/range;
                    case P95    -> (next.fit.p95    - prev.fit.p95)/range;
                    default     -> (next.fit.slo    - prev.fit.slo)/range;
                }; if(!Double.isFinite(cur.crowd)) cur.crowd=0; cur.crowd += delta;
            }
        }
    }

    // ===== Dominance & helpers =====

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
        double vol=(refE-minE)*(refP-minP)*(refS-minS); if (vol<=0) return 0.0; int hit=0; Random r = new Random(2025);
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
