package raj.cbm;

import raj.cbm.alg.*;
import raj.cbm.metrics.MetricsLogger;

import java.nio.file.Paths;
import java.util.List;
import java.util.Random;

/**
 * SingleRun
 * ---------
 * Code summary:
 * Orchestrates ONE full experiment for a given (scenario, algorithm, run, seed).
 * Steps:
 *   1) Build a {@link SimulationConfig} from the YAML-loaded {@link ScenarioConfig}.
 *   2) Create the chosen algorithm (MOPSO/MOGA/NSGA-II/SPEA-II/MOACO/Default).
 *   3) Run {@code optimize()} to get the best mapping & Pareto archive.
 *   4) Execute that mapping in CloudSim Plus ({@link SimulationRunner}) to measure metrics.
 *   5) Save the Pareto archive and return a compact {@link Result} for logging/CSV.
 *
 * Output layout (per run):
 *   <outRoot>/<scenario>/<algo>/run_<NN>/
 *     ├─ archive.csv              (Pareto front: energy/p95/slo and decision vars)
 *     ├─ cloudlets.csv            (finished cloudlets table)
 *     └─ cloudlets_metrics.csv    (derived metrics per cloudlet + energy totals)
 *
 * Notes:
 *   • Metrics returned: energy (kWh), P95 latency (ms), SLO-miss (%), HV, IGD, iterations to converge.
 *   • The "Default" path runs the baseline/non-optimized mapping.
 *
 * @author rchoudhary
 */
public class SingleRun {

    /** Minimal container with everything the caller typically wants to log. */
    public static class Result {
        public String scenarioName;
        public String algorithm;
        public int run;
        public long seed;
        public double energyKwh;
        public double p95Ms;
        public double sloPct;
        public double hv;
        public double igd;
        public int convergeIters;
        public int cloudletCount;
        public int vmCount;
        public int archiveSize;
        public String outDir;
    }

    /**
     * Execute one scenario–algorithm–run–seed combination end-to-end.
     *
     * @param cfg          Global scenario configuration loaded from YAML.
     * @param scenarioName Scenario key to select workload/limits/power.
     * @param run          1-based run index (for folder naming/logging).
     * @param seed         RNG seed for reproducibility.
     * @param algoName     One of: "mopso","moga","nsgaii","speaii","moaco","default".
     * @param outRoot      Root output folder for CSVs.
     * @return             A {@link Result} summarizing this run.
     * @throws Exception   If any file I/O or simulation error occurs.
     */
    public static Result execute(ScenarioConfig cfg, String scenarioName, int run, long seed, String algoName, String outRoot) throws Exception {
        // 1) Build a SimulationConfig from the scenario.
        SimulationConfig sc = SimulationConfig.fromScenario(cfg, scenarioName);

        // 2) Pick algorithm implementation.
        Algorithm algo = switch (algoName) {
            case "mopso" -> new MOPSO(sc, seed);
            case "moga" -> new MOGA(sc, seed);
            case "nsgaii" -> new NSGAII(sc, seed);
            case "speaii" -> new SPEAII(sc, seed);
            case "moaco" -> new MOACO(sc, seed);
            default -> new DefaultAlgo(sc, seed);
        };

        // Output folder: <outRoot>/<scenario>/<algo>/run_XX
        String prefix = String.format("%s/%s/%s/run_%02d", outRoot, scenarioName, algoName, run);
        java.nio.file.Files.createDirectories(java.nio.file.Paths.get(prefix));

        // 3) Optimize to get mapping + Pareto archive + indicators.
        Algorithm.Result opt = algo.optimize(prefix);

        // 4) Execute mapping in CloudSim Plus and collect runtime metrics.
        SimulationRunner.RunMetrics rm = SimulationRunner.runSimulation(sc, true, opt.mapping, prefix);

        // 5) Persist Pareto archive (for plots/HV/IGD audits).
        MetricsLogger.saveArchive(prefix + "/archive.csv", opt.paretoFront);

        // Build and return the summary object.
        Result r = new Result();
        r.scenarioName = scenarioName;
        r.algorithm = algoName;
        r.run = run;
        r.seed = seed;
        r.energyKwh = rm.energyKwh;
        r.p95Ms = rm.p95Ms;
        r.sloPct = rm.sloPct;
        r.hv = opt.hv;
        r.igd = opt.igd;
        r.convergeIters = opt.convergeIters;
        r.cloudletCount = sc.totalCloudlets();
        r.vmCount = sc.VMS;
        r.archiveSize = opt.paretoFront.size();
        r.outDir = prefix;
        return r;
    }
}
