package raj.cbm;

import java.util.Map;

/**
 * SimulationConfig
 * ----------------
 * Code summary:
 * A single place to hold all the knobs that define one simulation setup:
 *   • Host (physical server) shape: cores, CPU speed (MIPS), RAM, bandwidth, storage, power.
 *   • VM shape and counts.
 *   • Cloudlet (job) sizes.
 *   • How many jobs belong to each priority class (A/B/C/D = critical/high/medium/low).
 *   • Per-class deadlines and a convenience SLO threshold.
 *
 * Units:
 *   • RAM/Storage in MB, Bandwidth in Mbps, Power in Watts, CPU in MIPS, time in ms.
 *
 * How it’s used:
 *   • Call {@link #fromScenario(ScenarioConfig, String)} to build a ready-to-use config
 *     for a given scenario name from the YAML-loaded {@link ScenarioConfig}.
 *   • Then pass this config into factories (e.g., CloudResourceFactory) to build hosts, VMs, and jobs.
 *
 * Notes & safety:
 *   • Power: STATIC_POWER_W (idle) must be strictly less than MAX_POWER_W (full load),
 *     or construction should fail.
 *   • Priority counts: This code expects cfg.priorities to provide FRACTIONS (e.g., 0.2, 0.3, 0.3, 0.2)
 *     that sum to ~1.0. If your YAML stores WEIGHTS (e.g., 4.0, 2.0, 1.25, 1.0), convert to fractions first.
 *   • Rounding: counts are rounded; D_COUNT gets the remainder to ensure totals match the requested cloudlets.
 *
 * @author rchoudhary
 */
public class SimulationConfig {
    // Hosts (physical servers)
    public int HOSTS = 1;
    public int HOST_PES = 8;            // cores per host
    public int HOST_MIPS = 1000;        // per-core MIPS
    public int HOST_RAM = 16384;        // MB
    public long HOST_BW = 100_000;      // Mbps
    public long HOST_STORAGE = 1_000_000; // MB

    // VMs (virtual machines)  -> ADD THESE KNOBS
    public int VMS = 50;
    public int VM_PES = 2;              // vCPU per VM
    public int VM_MIPS = 1000;          // per-vCPU MIPS (used by new VmSimple(VM_MIPS, VM_PES))
    public int VM_RAM = 2048;           // MB
    public long VM_BW = 10_000;         // Mbps (10 Gbps)

    // Cloudlets (jobs)
    public int CLOUDLET_PES = 1;
    public int CLOUDLET_LENGTH = 10_000;

    // Priority counts
    public int A_COUNT, B_COUNT, C_COUNT, D_COUNT;

    // Deadlines
    public int DL_CRITICAL_MS, DL_HIGH_MS, DL_MED_MS, DL_LOW_MS;
    public double SLO_THRESHOLD_MS; // currently unused

    public static SimulationConfig fromScenario(ScenarioConfig cfg, String scenarioName) {
        SimulationConfig sc = new SimulationConfig();

        // Workload selection
        ScenarioConfig.Workload wl = cfg.workloads.stream()
                .filter(w -> w.name.equalsIgnoreCase(scenarioName))
                .findFirst().orElseThrow();

        sc.VMS = wl.vms;

        // Split cloudlets
        int n = wl.cloudlets;
        sc.A_COUNT = (int)Math.round(n * cfg.priorities.get("critical"));
        sc.B_COUNT = (int)Math.round(n * cfg.priorities.get("high"));
        sc.C_COUNT = (int)Math.round(n * cfg.priorities.get("medium"));
        sc.D_COUNT = n - sc.A_COUNT - sc.B_COUNT - sc.C_COUNT;

        // Deadlines
        sc.DL_CRITICAL_MS = cfg.deadlines_ms.get("critical");
        sc.DL_HIGH_MS     = cfg.deadlines_ms.get("high");
        sc.DL_MED_MS      = cfg.deadlines_ms.get("medium");
        sc.DL_LOW_MS      = cfg.deadlines_ms.get("low");

        // Power
        sc.MAX_POWER_W    = cfg.power_model.max_power_w;
        sc.STATIC_POWER_W = cfg.power_model.static_power_w;
        if (sc.STATIC_POWER_W >= sc.MAX_POWER_W) {
            throw new IllegalArgumentException(String.format(
                "Bad power params: max=%.1fW, static=%.1fW", sc.MAX_POWER_W, sc.STATIC_POWER_W));
        }

        // Cloudlet size per scenario
        if ("low".equalsIgnoreCase(scenarioName))      sc.CLOUDLET_LENGTH = 8_000;
        else if ("medium".equalsIgnoreCase(scenarioName)) sc.CLOUDLET_LENGTH = 10_000;
        else                                           sc.CLOUDLET_LENGTH = 12_000;

        // OPTIONAL: scale host count for large VM counts (keeps ~50 VMs/host)
        if (sc.VMS > 100 && sc.HOSTS == 1) {
            sc.HOSTS = Math.max( (int)Math.ceil(sc.VMS / 50.0), 2 );
        }

        return sc;
    }

    public int totalCloudlets() {
        return A_COUNT + B_COUNT + C_COUNT + D_COUNT;
    }

    // Power model fields
    public double MAX_POWER_W;
    public double STATIC_POWER_W;
}
