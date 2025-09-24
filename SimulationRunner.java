package raj.cbm;

import org.cloudsimplus.brokers.DatacenterBrokerSimple;
import org.cloudsimplus.cloudlets.Cloudlet;
import org.cloudsimplus.core.CloudSimPlus;
import org.cloudsimplus.datacenters.Datacenter;
import org.cloudsimplus.vms.Vm;
import raj.cbm.metrics.MetricsLogger;

import java.util.List;

/**
 * SimulationRunner
 * ----------------
 * Code summary:
 * Runs one complete CloudSim Plus experiment using the sizes/counts from
 * {@link SimulationConfig}. Optionally applies a job→VM mapping you provide.
 * After the run, it writes CSVs and returns a few easy-to-read metrics.
 *
 * What happens:
 *   1) Build datacenter/VMs/cloudlets via {@link CloudResourceFactory}.
 *   2) If a mapping is provided, bind each cloudlet to a specific VM.
 *   3) Start the simulation and wait until everything finishes.
 *   4) Save per-cloudlet results to CSV and compute summary metrics.
 *
 * Units:
 *   • energyKwh in kWh, p95Ms in milliseconds, sloPct in percent (0–100).
 *
 * Notes on mapping:
 *   • If {@code useMapping == true} and {@code mapping != null}, cloudlet i is
 *     bound to VM index {@code floorMod(mapping[i], vmCount)}. Negative or
 *     oversized indices are safely wrapped by {@code floorMod}.
 *
 * @author rchoudhary
 */
public class SimulationRunner {

    /** Minimal result bundle with the key metrics we plot/compare. */
    public static class RunMetrics {
        /** Total datacenter energy in kWh. */
        public double energyKwh;
        /** 95th percentile (tail) latency across cloudlets in ms. */
        public double p95Ms;
        /** Percent of cloudlets that missed their deadline (0–100). */
        public double sloPct;
    }

    /**
     * Execute one simulation.
     *
     * @param sc         The simulation shape (hosts, VMs, cloudlets, power).
     * @param useMapping If true, bind cloudlets to VMs using {@code mapping}.
     * @param mapping    An int array of VM indices (one per cloudlet). Values are wrapped
     *                   with {@code floorMod} so negatives/large numbers won’t crash.
     * @param prefix     Output folder prefix for CSVs.
     * @return           Aggregated metrics (energy, P95, SLO%).
     */
    public static RunMetrics runSimulation(SimulationConfig sc, boolean useMapping, int[] mapping, String prefix) {
        // 1) Bootstrap simulator and build resources from config.
        CloudSimPlus simulation = new CloudSimPlus();
        CloudResourceFactory factory = new CloudResourceFactory(sc);
        Datacenter datacenter = factory.createDatacenter(simulation);
        DatacenterBrokerSimple broker = new DatacenterBrokerSimple(simulation);
        List<Vm> vmList = factory.createVms();
        List<Cloudlet> cloudletList = factory.createCloudlets();

        // Register VMs with the broker before submitting cloudlets.
        broker.submitVmList(vmList);

        // 2) Optional: apply a fixed mapping cloudlet i -> VM mapping[i] (wrapped to range).
        if (useMapping && mapping != null) {
            for (int i = 0; i < cloudletList.size(); i++) {
                int vmIndex = Math.floorMod(mapping[i], vmList.size());
                broker.bindCloudletToVm(cloudletList.get(i), vmList.get(vmIndex));
            }
        }

        // Submit cloudlets (either free to any VM, or pre-bound above).
        broker.submitCloudletList(cloudletList);

        // 3) Run the simulation to completion.
        simulation.start();

        // 4) Persist detailed results and compute summary metrics.
        var finished = broker.getCloudletFinishedList();
        MetricsLogger.saveCloudletsTable(finished, prefix + "/cloudlets.csv");
        // Build a deadline map using cloudlet priority (4,3,2,1 -> critical, high, medium, low)
        java.util.Map<Long,Integer> deadlineMsById = new java.util.HashMap<>();
        for (var c : finished) {
            int prio = c.getPriority();
            int dl;
            if (prio >= 4) dl = sc.DL_CRITICAL_MS;
            else if (prio == 3) dl = sc.DL_HIGH_MS;
            else if (prio == 2) dl = sc.DL_MED_MS;
            else dl = sc.DL_LOW_MS;
            deadlineMsById.put(c.getId(), dl);
        }
        MetricsLogger.saveCloudletsMetricsCsv(finished, datacenter, prefix + "/cloudlets_metrics.csv", deadlineMsById);

        // Build and return the summary object used by higher-level runners.
        RunMetrics rm = new RunMetrics();
        MetricsLogger.Summary s = MetricsLogger.computeSummary(finished, datacenter, deadlineMsById);
        rm.energyKwh = s.energyKwh;
        rm.p95Ms = s.p95Ms;
        rm.sloPct = s.sloPct;
        return rm;
    }
}
