package raj.cbm.alg;

import raj.cbm.SimulationConfig;

import java.util.*;

/**
 * MOPSO (deadline-aware, priority-weighted) for CBM mapping
 * ---------------------------------------------------------
 * Encoding:
 *   • Particle position = double[N] in [0, VM); decode by floor() → VM index per cloudlet.
 *
 * Objectives (minimize):
 *   • energy (kWh), p95 latency (ms), weighted SLO-miss (%).
 *     Weighted SLO counts misses from higher-priority classes more:
 *       D=1.0, C=1.5, B=2.5, A=4.0  (tunable via CLASS_WEIGHTS)
 *
 * Tailoring vs. vanilla MOPSO:
 *   • Fitness: SLO is priority-weighted.
 *   • gBest: sampled deadline-aware (strong bias to lower weighted SLO, then lower P95),
 *            still diversity-aware (crowding as multiplier).
 *   • pBest: when neither dominates, lexicographic tie-break SLO → P95 → Energy.
 *
 * Notes:
 *   • Uses a fast mapping-dependent surrogate for speed. You can swap evaluate(...) to call
 *     SimulationRunner for ground-truth (much slower) if desired.
 *
 * @author rchoudhary
 */
public class MOPSO implements Algorithm {
    private final SimulationConfig sc;
    private final Random rnd;

    // ===== Hyperparameters =====
    private int swarmSize;
    private int maxIters;
    private double wStart = 0.72;   // inertia start
    private double wEnd   = 0.32;   // inertia end
    private double c1 = 1.7;        // cognitive
    private double c2 = 1.7;        // social
    private int archiveCapacity = 120;
    private double vMax;            // velocity clamp per dim
    private double mutationProb = 0.06;

    // Class urgency weights for weighted SLO: indices 0..3 = D,C,B,A (see cls[] mapping below)
    private static final double[] CLASS_WEIGHTS = {1.0, 1.5, 2.5, 4.0};

    public MOPSO(SimulationConfig sc, long seed) {
        this.sc = sc;
        this.rnd = new Random(seed ^ 3800162996182131818L);
        int n = Math.max(1, sc.totalCloudlets());
        int vm = Math.max(1, sc.VMS);
        this.swarmSize = Math.max(24, Math.min(64, 8 + n / 8));
        this.maxIters  = Math.max(30, Math.min(100, 20 + n / 10));
        this.vMax = Math.max(1.0, vm - 1.0);
    }

    @Override
    public Result optimize(String outPrefix) {
        final int N  = Math.max(1, sc.totalCloudlets());
        final int VM = Math.max(1, sc.VMS);

        // Build class vector aligned with CloudResourceFactory ordering:
        // A then B then C then D → we store as 3=A, 2=B, 1=C, 0=D
        final int[] cls = new int[N];
        int idx = 0;
        idx = fill(cls, idx, sc.A_COUNT, 3);
        idx = fill(cls, idx, sc.B_COUNT, 2);
        idx = fill(cls, idx, sc.C_COUNT, 1);
        fill(cls, idx, sc.D_COUNT, 0);

        // Deadlines array remapped to indices 0..3 matching cls[] above (D,C,B,A)
        final int[] dlMs = new int[]{ sc.DL_LOW_MS, sc.DL_MED_MS, sc.DL_HIGH_MS, sc.DL_CRITICAL_MS };

        // Swarm init
        List<Particle> swarm = new ArrayList<>(swarmSize);
        for (int i = 0; i < swarmSize; i++) {
            Particle p = new Particle(N, VM, rnd);
            p.fit = evaluate(p.mapping(VM), cls, dlMs);
            p.pBestPos = p.pos.clone();
            p.pBestFit = p.fit;
            swarm.add(p);
        }

        Archive archive = new Archive(archiveCapacity);
        for (Particle p : swarm) archive.add(p.copySolution(VM));

        double lastHV = hypervolumeMC(archive.front, 2000);
        int convergeIter = 0;

        // Main loop
        for (int iter = 1; iter <= maxIters; iter++) {
            double w = wStart + (wEnd - wStart) * (iter / (double) maxIters);

            // Deadline-aware gBest (biased to low weighted-SLO, then low P95, x diversity)
            Solution gBest = archive.sampleDeadlineAware(rnd);
            double[] gPos = (gBest != null) ? gBest.pos : swarm.get(rnd.nextInt(swarm.size())).pos;

            for (Particle p : swarm) {
                // Velocity/position update
                for (int d = 0; d < N; d++) {
                    double r1 = rnd.nextDouble(), r2 = rnd.nextDouble();
                    double cognitive = c1 * r1 * (p.pBestPos[d] - p.pos[d]);
                    double social    = c2 * r2 * (gPos[d]     - p.pos[d]);
                    p.vel[d] = clamp(w * p.vel[d] + cognitive + social, -vMax, vMax);
                    p.pos[d] = clamp(p.pos[d] + p.vel[d], 0.0, VM - 1e-6);

                    // Light mutation/turbulence
                    if (rnd.nextDouble() < mutationProb) {
                        p.pos[d] = clamp(p.pos[d] + rnd.nextGaussian() * 0.5, 0.0, VM - 1e-6);
                    }
                }

                // Evaluate and update pBest with deadline-aware tie-break
                Fitness f = evaluate(p.mapping(VM), cls, dlMs);
                p.fit = f;
                if (dominates(f, p.pBestFit) || equalsFit(f, p.pBestFit)
                    || (!dominates(p.pBestFit, f) && lexBetterDeadlineFirst(f, p.pBestFit))) {
                    p.pBestFit = f;
                    p.pBestPos = p.pos.clone();
                }

                // Push to archive
                archive.add(new Solution(p.pos.clone(), p.mapping(VM), f));
            }

            // Convergence monitor via HV
            double hv = hypervolumeMC(archive.front, 2000);
            if (hv > lastHV + 1e-6) { lastHV = hv; convergeIter = iter; }
        }

        // Pick deploy candidate: lowest weighted-SLO, then lowest P95, then lowest energy
        Solution deploy = archive.front.stream()
                .min((a, b) -> {
                    int c = Double.compare(a.fit.slo, b.fit.slo);
                    if (c != 0) return c;
                    c = Double.compare(a.fit.p95, b.fit.p95);
                    if (c != 0) return c;
                    return Double.compare(a.fit.energy, b.fit.energy);
                })
                .orElseGet(() -> swarmBest(swarm, VM));

        // Pack results
        Result r = new Result();
        r.mapping = deploy.mapping;
        r.paretoFront = new ArrayList<>(archive.front.size());
        for (Solution s : archive.front) r.paretoFront.add(new double[]{s.fit.energy, s.fit.p95, s.fit.slo});
        r.hv = hypervolumeMC(archive.front, 8000);
        r.igd = igdToIdeal(archive.front);
        r.convergeIters = convergeIter;
        return r;
    }

    // ===== Fitness & evaluation (priority-weighted SLO) =====

    private Fitness evaluate(int[] mapping, int[] cls, int[] dlMs) {
        final int N = mapping.length;
        final int VM = Math.max(1, sc.VMS);
        final double jobServiceSec =
                sc.CLOUDLET_LENGTH / (double) (Math.max(1, sc.VM_PES) * Math.max(1, sc.HOST_MIPS));

        // Per-VM sequential queues (proxy)
        int[] counts = new int[VM];
        double[] lastFinish = new double[VM];
        double[] compMs = new double[N];
        for (int i = 0; i < N; i++) {
            int vm = Math.floorMod(mapping[i], VM);
            int k = counts[vm]++;
            double finishMs = (k + 1) * jobServiceSec * 1000.0;
            compMs[i] = lastFinish[vm] = finishMs;
        }

        // P95
        double[] sorted = compMs.clone();
        Arrays.sort(sorted);
        double p95 = sorted[(int) Math.floor(0.95 * (sorted.length - 1))];

        // Weighted SLO (%): A,B miss cost > C,D
        double missWeighted = 0.0, weightSum = 0.0;
        for (int i = 0; i < N; i++) {
            int c = cls[i]; // 3=A .. 0=D
            double w = CLASS_WEIGHTS[c];
            weightSum += w;
            if (compMs[i] > dlMs[c]) missWeighted += w;
        }
        double sloWeightedPct = 100.0 * missWeighted / Math.max(1e-9, weightSum);

        // Energy via mean utilization over makespan
        int HOSTS = Math.max(1, sc.HOSTS);
        int HOST_PES = Math.max(1, sc.HOST_PES);
        double makespanSec = Arrays.stream(lastFinish).max().orElse(0.0) / 1000.0;
        double totalCoreSeconds = N * (sc.CLOUDLET_LENGTH / (double) Math.max(1, sc.HOST_MIPS));
        double avgUtil = Math.min(1.0, totalCoreSeconds / (HOSTS * HOST_PES * Math.max(1e-9, makespanSec)));
        double pIdle = sc.STATIC_POWER_W, pMax = sc.MAX_POWER_W;
        double pMeanPerHost = pIdle + (pMax - pIdle) * avgUtil;
        double energyKwh = (HOSTS * pMeanPerHost * makespanSec) / 3600.0 / 1000.0;

        return new Fitness(energyKwh, p95, sloWeightedPct);
    }

    // ===== Archive & dominance =====

    private static class Fitness {
        final double energy, p95, slo;
        Fitness(double e, double p, double s) { this.energy = e; this.p95 = p; this.slo = s; }
    }

    private static class Solution {
        final double[] pos;   // continuous position (for PSO)
        final int[] mapping;  // decoded mapping
        final Fitness fit;
        double crowd;         // crowding distance
        Solution(double[] pos, int[] map, Fitness f) { this.pos = pos; this.mapping = map; this.fit = f; }
    }

    private static class Archive {
        final int capacity;
        final List<Solution> front = new ArrayList<>();
        Archive(int capacity) { this.capacity = capacity; }

        void add(Solution s) {
            // remove any solutions dominated by s
            front.removeIf(o -> dominates(s.fit, o.fit) && !equalsFit(s.fit, o.fit));
            // only insert if not dominated by someone inside
            for (Solution o : front) {
                if (dominates(o.fit, s.fit) && !equalsFit(o.fit, s.fit)) return;
            }
            front.add(s);
            if (front.size() > capacity) pruneByCrowding();
        }

        /** Diversity-only sampler (unused now, kept for reference). */
        Solution sampleByCrowding(Random rnd) {
            if (front.isEmpty()) return null;
            computeCrowding();
            double sum = 0.0;
            for (Solution s : front) sum += s.crowd > 0 ? s.crowd : 0.0;
            if (sum <= 0) return front.get(rnd.nextInt(front.size()));
            double r = rnd.nextDouble() * sum, acc = 0;
            for (Solution s : front) { acc += Math.max(0.0, s.crowd); if (acc >= r) return s; }
            return front.get(front.size() - 1);
        }

        /** Deadline-aware gBest sampler: bias to low weighted-SLO, then low P95; crowding boosts diversity. */
        Solution sampleDeadlineAware(Random rnd) {
            if (front.isEmpty()) return null;
            computeCrowding();

            double minS = Double.POSITIVE_INFINITY, maxS = 0, minP = Double.POSITIVE_INFINITY, maxP = 0;
            for (Solution s : front) {
                minS = Math.min(minS, s.fit.slo);  maxS = Math.max(maxS, s.fit.slo);
                minP = Math.min(minP, s.fit.p95);  maxP = Math.max(maxP, s.fit.p95);
            }
            double rangeS = Math.max(1e-9, maxS - minS);
            double rangeP = Math.max(1e-9, maxP - minP);

            final double alpha = 4.0, beta = 1.0; // emphasis on deadlines vs latency
            double[] score = new double[front.size()];
            double sum = 0.0;
            for (int i = 0; i < front.size(); i++) {
                Solution s = front.get(i);
                double sloN = (s.fit.slo - minS) / rangeS; // 0 best
                double p95N = (s.fit.p95 - minP) / rangeP; // 0 best
                double crowdBoost = 1.0 + Math.max(0.0, s.crowd);
                double sc = crowdBoost * Math.exp(-alpha * sloN) * Math.exp(-beta * p95N);
                score[i] = sc; sum += sc;
            }
            double r = rnd.nextDouble() * (sum > 0 ? sum : 1.0);
            double acc = 0.0;
            for (int i = 0; i < front.size(); i++) { acc += score[i]; if (acc >= r) return front.get(i); }
            return front.get(front.size() - 1);
        }

        private void pruneByCrowding() {
            computeCrowding();
            Solution worst = front.stream().min(Comparator.comparingDouble(s -> s.crowd)).orElse(null);
            if (worst != null) front.remove(worst);
        }

        private void computeCrowding() {
            if (front.size() < 3) { for (Solution s : front) s.crowd = Double.POSITIVE_INFINITY; return; }
            computeCrowdingFor(Objective.ENERGY);
            computeCrowdingFor(Objective.P95);
            computeCrowdingFor(Objective.SLO);
        }

        enum Objective { ENERGY, P95, SLO }

        private void computeCrowdingFor(Objective obj) {
            Comparator<Solution> cmp = switch (obj) {
                case ENERGY -> Comparator.comparingDouble(s -> s.fit.energy);
                case P95    -> Comparator.comparingDouble(s -> s.fit.p95);
                case SLO    -> Comparator.comparingDouble(s -> s.fit.slo);
            };
            front.sort(cmp);
            front.get(0).crowd = front.get(front.size() - 1).crowd = Double.POSITIVE_INFINITY;

            double min, max;
            switch (obj) {
                case ENERGY -> { min = front.get(0).fit.energy; max = front.get(front.size() - 1).fit.energy; }
                case P95    -> { min = front.get(0).fit.p95;    max = front.get(front.size() - 1).fit.p95; }
                default     -> { min = front.get(0).fit.slo;    max = front.get(front.size() - 1).fit.slo; }
            }
            double range = Math.max(1e-12, max - min);
            for (int i = 1; i < front.size() - 1; i++) {
                Solution prev = front.get(i - 1), next = front.get(i + 1), cur = front.get(i);
                double delta = switch (obj) {
                    case ENERGY -> (next.fit.energy - prev.fit.energy) / range;
                    case P95    -> (next.fit.p95    - prev.fit.p95)    / range;
                    default     -> (next.fit.slo    - prev.fit.slo)    / range;
                };
                if (!Double.isFinite(cur.crowd)) cur.crowd = 0;
                cur.crowd += delta;
            }
        }
    }

    // ===== Particle =====

    private static class Particle {
        final double[] pos;   // [0, VM)
        final double[] vel;   // bounded by vMax
        Fitness fit;
        double[] pBestPos;
        Fitness pBestFit;
        private final Random rnd;

        Particle(int dims, int VM, Random rnd) {
            this.rnd = rnd;
            this.pos = new double[dims];
            this.vel = new double[dims];
            for (int d = 0; d < dims; d++) {
                pos[d] = rnd.nextDouble() * Math.max(1, VM) - 1e-6;
                vel[d] = rnd.nextGaussian() * 0.5;
            }
        }

        int[] mapping(int VM) {
            int[] m = new int[pos.length];
            for (int i = 0; i < pos.length; i++) {
                int v = (int) Math.floor(pos[i]);
                if (v < 0) v = 0; else if (v >= VM) v = VM - 1;
                m[i] = v;
            }
            return m;
        }

        Solution copySolution(int VM) {
            return new Solution(pos.clone(), mapping(VM), fit);
        }
    }

    // ===== Utilities =====

    private static int fill(int[] a, int start, int count, int val) {
        for (int i = 0; i < count && start + i < a.length; i++) a[start + i] = val;
        return Math.min(a.length, start + count);
    }

    private static boolean dominates(Fitness a, Fitness b) {
        boolean betterOrEqual = a.energy <= b.energy && a.p95 <= b.p95 && a.slo <= b.slo;
        boolean strictlyBetter = a.energy < b.energy || a.p95 < b.p95 || a.slo < b.slo;
        return betterOrEqual && strictlyBetter;
    }

    private static boolean equalsFit(Fitness a, Fitness b) {
        return Math.abs(a.energy - b.energy) < 1e-12
            && Math.abs(a.p95 - b.p95) < 1e-12
            && Math.abs(a.slo - b.slo) < 1e-12;
    }

    private static boolean lexBetterDeadlineFirst(Fitness a, Fitness b) {
        int c = Double.compare(a.slo, b.slo);  if (c != 0) return c < 0; // lower weighted SLO
        c = Double.compare(a.p95, b.p95);      if (c != 0) return c < 0; // then lower P95
        return a.energy < b.energy;                                     // then lower energy
    }

    private static double clamp(double x, double lo, double hi) {
        return (x < lo) ? lo : (x > hi) ? hi : x;
    }

    private static Solution swarmBest(List<Particle> swarm, int VM) {
        Particle best = null;
        for (Particle p : swarm) {
            if (best == null) { best = p; continue; }
            Fitness a = p.fit, b = best.fit;
            if (dominates(a, b)) best = p;
            else if (!dominates(b, a)) {
                // tie-break SLO→P95→Energy
                int c = Double.compare(a.slo, b.slo);
                if (c < 0) best = p;
                else if (c == 0) {
                    c = Double.compare(a.p95, b.p95);
                    if (c < 0) best = p;
                    else if (c == 0 && a.energy < b.energy) best = p;
                }
            }
        }
        if (best == null) return null;
        return new Solution(best.pos.clone(), best.mapping(Math.max(1, VM)), best.fit);
    }

    // ===== Indicators (simple approximations) =====

    /** Monte-Carlo HV (reference point slightly worse than the worst observed). */
    private static double hypervolumeMC(List<Solution> front, int samples) {
        if (front.isEmpty()) return 0.0;
        double minE = Double.POSITIVE_INFINITY, minP = Double.POSITIVE_INFINITY, minS = Double.POSITIVE_INFINITY;
        double maxE = 0, maxP = 0, maxS = 0;
        for (Solution s : front) {
            minE = Math.min(minE, s.fit.energy); maxE = Math.max(maxE, s.fit.energy);
            minP = Math.min(minP, s.fit.p95);    maxP = Math.max(maxP, s.fit.p95);
            minS = Math.min(minS, s.fit.slo);    maxS = Math.max(maxS, s.fit.slo);
        }
        double refE = maxE * 1.2 + 1e-6, refP = maxP * 1.2 + 1e-6, refS = maxS * 1.2 + 1e-6;
        double vol = (refE - minE) * (refP - minP) * (refS - minS);
        if (vol <= 0) return 0.0;
        int hit = 0;
        Random r = new Random(1337);
        for (int i = 0; i < samples; i++) {
            double e = minE + r.nextDouble() * (refE - minE);
            double p = minP + r.nextDouble() * (refP - minP);
            double s = minS + r.nextDouble() * (refS - minS);
            if (isDominatedByFront(e, p, s, front)) hit++;
        }
        return (hit / (double) samples) * vol;
    }

    private static boolean isDominatedByFront(double e, double p, double s, List<Solution> front) {
        for (Solution x : front) {
            if (x.fit.energy <= e && x.fit.p95 <= p && x.fit.slo <= s
                    && (x.fit.energy < e || x.fit.p95 < p || x.fit.slo < s)) return true;
        }
        return false;
    }

    /** IGD to ideal point (min(E), min(P95), min(SLO)) — simple, bounded > 0. */
    private static double igdToIdeal(List<Solution> front) {
        if (front.isEmpty()) return 0.0;
        double minE = Double.POSITIVE_INFINITY, minP = Double.POSITIVE_INFINITY, minS = Double.POSITIVE_INFINITY;
        for (Solution s : front) {
            minE = Math.min(minE, s.fit.energy);
            minP = Math.min(minP, s.fit.p95);
            minS = Math.min(minS, s.fit.slo);
        }
        double acc = 0.0;
        for (Solution s : front) {
            double de = s.fit.energy - minE, dp = s.fit.p95 - minP, ds = s.fit.slo - minS;
            acc += Math.sqrt(Math.max(0, de * de + dp * dp + ds * ds));
        }
        return acc / front.size();
    }
}
