package raj.cbm;

import org.yaml.snakeyaml.Yaml;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * ScenarioConfig
 * --------------
 * Code summary:
 * Loads experiment settings from a YAML file and exposes them as simple fields.
 * Think of it as the “what to run” file:
 *   • which workloads to simulate (name, #jobs/cloudlets, #VMs),
 *   • how urgent each class is (priority weights),
 *   • per-class deadlines (in milliseconds),
 *   • host power model numbers (Watts),
 *   • and how many repeated runs to do per (scenario × algorithm) combo.
 *
 * Expected YAML structure (example):
 *
 * workloads:
 *   - name: medium
 *     cloudlets: 120
 *     vms: 24
 *   - name: heavy
 *     cloudlets: 200
 *     vms: 40
 * priorities:
 *   A: 4.0
 *   B: 2.0
 *   C: 1.25
 *   D: 1.0
 * deadlines_ms:
 *   A: 200
 *   B: 300
 *   C: 500
 *   D: 800
 * power_model:
 *   max_power_w: 200.0
 *   static_power_w: 50.0
 * runs_per_combo: 10
 *
 * Notes:
 * - Units: deadlines in ms; power in Watts; counts are integers.
 * - Power sanity: max_power_w should be > static_power_w (idle).
 *
 * @author rchoudhary
 */
public class ScenarioConfig {

    /** One workload “scenario”: name + how many cloudlets and VMs to use. */
    public static class Workload {
        public String name;
        public int cloudlets;
        public int vms;
    }

    /** Host power model parameters (Watts). */
    public static class PowerModel {
        public double max_power_w;
        public double static_power_w;
    }

    /** List of workloads to run (parsed from YAML 'workloads'). */
    public List<Workload> workloads;

    /** Priority weights per class label (e.g., A,B,C,D) from YAML 'priorities'. */
    public Map<String, Double> priorities;

    /** Deadlines (ms) per class label from YAML 'deadlines_ms'. */
    public Map<String, Integer> deadlines_ms;

    /** Power model settings for hosts from YAML 'power_model'. */
    public PowerModel power_model;

    /** Number of repeats for each (scenario × algorithm) from YAML 'runs_per_combo'. */
    public int runsPerCombo;

    /**
     * Load a ScenarioConfig from a YAML file at 'path'.
     * Throws an exception if the file can’t be read or fields are missing/mismatched.
     */
    @SuppressWarnings("unchecked")
    public static ScenarioConfig load(String path) throws Exception {
        Yaml yaml = new Yaml();
        try (InputStream is = Files.newInputStream(Paths.get(path))) {
            Map<String, Object> map = yaml.load(is);

            ScenarioConfig cfg = new ScenarioConfig();
            cfg.workloads = new ArrayList<>();

            // workloads: list of {name, cloudlets, vms}
            List<Map<String, Object>> wls = (List<Map<String, Object>>) map.get("workloads");
            for (Map<String, Object> w : wls) {
                Workload wl = new Workload();
                wl.name      = String.valueOf(w.get("name"));
                wl.cloudlets = ((Number) w.get("cloudlets")).intValue();
                wl.vms       = ((Number) w.get("vms")).intValue();
                cfg.workloads.add(wl);
            }

            // priorities: map like {"A": 4.0, "B": 2.0, ...}
            cfg.priorities = (Map<String, Double>) map.get("priorities");

            // deadlines_ms: map like {"A": 200, "B": 300, ...} in milliseconds
            cfg.deadlines_ms = (Map<String, Integer>) map.get("deadlines_ms");

            // power_model: {max_power_w: 200.0, static_power_w: 50.0}
            Map<String, Object> pm = (Map<String, Object>) map.get("power_model");
            PowerModel p = new PowerModel();
            p.max_power_w    = ((Number) pm.get("max_power_w")).doubleValue();
            p.static_power_w = ((Number) pm.get("static_power_w")).doubleValue();
            cfg.power_model = p;

            // runs_per_combo: integer
            cfg.runsPerCombo = ((Number) map.get("runs_per_combo")).intValue();

            return cfg;
        }
    }
}
