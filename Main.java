/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package raj.cbm;

/**
 * Main
 * ----
 * Code summary:
 * This is the launcher for all experiments. It reads a scenario file and a list
 * of random seeds, then runs **every algorithm** on **every scenario** for a
 * fixed number of runs. After each run it appends one summary row to CSV files
 * (per-scenario and an all-scenarios file).
 *
 * Command-line flags (with defaults):
 *   --scenarios <path>   YAML with workloads & run counts (default: config/scenarios.yaml)
 *   --seeds     <path>   Text file of seeds (default: config/seeds.txt)
 *   --out       <dir>    Output folder for CSVs (default: results/)
 *
 * Output files:
 *   results/summary_all.csv
 *   results/summary_<scenarioName>.csv
 * Each row has metrics like energy (kWh), P95 latency (ms), and SLO-miss (%).
 *
 * Example:
 *   java -jar target/cloudcbm.jar \
 *     --scenarios config/scenarios.yaml \
 *     --seeds config/seeds.txt \
 *     --out results/
 */

/**
 *
 * @author rchoudhary
 */

import raj.cbm.util.Seeds;
import raj.cbm.metrics.MetricsLogger;
import raj.cbm.alg.*;
import org.yaml.snakeyaml.Yaml;

import java.io.*;
import java.nio.file.*;
import java.time.LocalDateTime;
import java.util.*;

/**
 * Main entry point: loops scenarios × algorithms × runs, executes simulations,
 * and appends a summary row for each run.
 *
 * Usage:
 *   java -jar target/cloudcbm.jar --scenarios config/scenarios.yaml --seeds config/seeds.txt --out results/
 */
public class Main {

    public static void main(String[] args) throws Exception {
        // Parse CLI flags or fall back to sensible defaults.
        Map<String, String> cli = parseArgs(args);
        final String cfgPath  = cli.getOrDefault("--scenarios", "config/scenarios.yaml");
        final String seedsPath = cli.getOrDefault("--seeds", "config/seeds.txt");
        final String outDir   = cli.getOrDefault("--out", "results");

        // Ensure output directory exists.
        Files.createDirectories(Paths.get(outDir));

        // Load scenario configuration and seeds matrix: [scenarioIndex][runIndex].
        ScenarioConfig config = ScenarioConfig.load(cfgPath);
        long[][] seeds = Seeds.load(seedsPath, config.workloads.size(), config.runsPerCombo);

        // Algorithms to evaluate (names used by SingleRun).
        String[] algos = new String[]{"default", "mopso", "moga", "nsgaii", "speaii", "moaco"};

        // Prepare CSVs: one global file and one per scenario.
        String summaryAll = Paths.get(outDir, "summary_all.csv").toString();
        MetricsLogger.initSummary(summaryAll);

        int scenarioIndex = 0;
        for (ScenarioConfig.Workload wl : config.workloads) {
            final String scenarioName = wl.name;
            final String scenarioSummary = Paths.get(outDir, "summary_" + scenarioName + ".csv").toString();
            MetricsLogger.initSummary(scenarioSummary);

            // Repeat each scenario for the configured number of runs.
            for (int run = 1; run <= config.runsPerCombo; run++) {
                long seed = seeds[scenarioIndex][run - 1];

                // Try every algorithm for this scenario+run.
                for (String algoName : algos) {
                    SingleRun.Result rr = SingleRun.execute(config, scenarioName, run, seed, algoName, outDir);

                    // Append results to both CSVs.
                    MetricsLogger.appendSummary(scenarioSummary, rr);
                    MetricsLogger.appendSummary(summaryAll, rr);

                    System.out.printf(
                        "Completed %s | %s | run %d | seed=%d | E=%.3f kWh, P95=%.1f ms, SLO=%.2f%%%n",
                        scenarioName, algoName, run, seed, rr.energyKwh, rr.p95Ms, rr.sloPct
                    );
                }
            }
            scenarioIndex++;
        }

        System.out.println("All runs completed at " + LocalDateTime.now());
    }

    /**
     * Minimal flag parser: looks for pairs like ["--out", "results"] and
     * stores them in a map. Unknown flags are ignored. If a flag is provided
     * without a value, it is skipped.
     */
    private static Map<String, String> parseArgs(String[] args) {
        Map<String, String> m = new HashMap<>();
        for (int i = 0; i < args.length - 1; i++) {
            if (args[i].startsWith("--")) m.put(args[i], args[i + 1]);
        }
        return m;
    }
}
