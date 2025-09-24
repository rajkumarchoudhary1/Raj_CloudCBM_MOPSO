package raj.cbm;

import org.cloudsimplus.cloudlets.Cloudlet;
import org.cloudsimplus.cloudlets.CloudletSimple;
import org.cloudsimplus.core.CloudSimPlus;
import org.cloudsimplus.datacenters.Datacenter;
import org.cloudsimplus.datacenters.DatacenterSimple;
import org.cloudsimplus.hosts.Host;
import org.cloudsimplus.hosts.HostSimple;
import org.cloudsimplus.power.models.PowerModelHostSimple;
import org.cloudsimplus.resources.Pe;
import org.cloudsimplus.resources.PeSimple;
import org.cloudsimplus.utilizationmodels.UtilizationModelDynamic;
import org.cloudsimplus.vms.Vm;
import org.cloudsimplus.vms.VmSimple;

import java.util.ArrayList;
import java.util.List;

/**
 * CloudResourceFactory
 * --------------------
 * Plain-English summary:
 * Builds a small “mini-cloud” for the simulation using values from {@link SimulationConfig}.
 * It creates:
 *   • A datacenter: a list of physical servers (hosts) with CPU cores and a power model
 *   • A pool of VMs: virtual machines that will run the jobs
 *   • A list of jobs (cloudlets): grouped by urgency A/B/C/D with different average CPU use
 *
 * Key ideas:
 *   • Host = physical server (RAM/BW/Storage in MB/Mbps/MB; CPU speed in MIPS per core/PE)
 *   • VM   = a slice of a host that runs cloudlets
 *   • Cloudlet = one job; we set its priority (4..1) and typical CPU utilization (e.g., 80%)
 *
 * Units & safety:
 *   • RAM/Storage are MB; Bandwidth is Mbps; Power is Watts; CPU is MIPS.
 *   • Ensure sc.MAX_POWER_W > sc.STATIC_POWER_W, or CloudSim will throw an error.
 *
 * @author rchoudhary
 */
public class CloudResourceFactory {

    /** All sizes/counts (hosts, VMs, MIPS, power, etc.) come from this config. */
    private final SimulationConfig sc;

    public CloudResourceFactory(SimulationConfig sc) {
        this.sc = sc;
    }

    /**
     * Create a datacenter by instantiating 'sc.HOSTS' physical servers.
     */
    public Datacenter createDatacenter(CloudSimPlus sim) {
        final var hostList = new ArrayList<Host>(sc.HOSTS);
        for (int i = 0; i < sc.HOSTS; i++) {
            final var host = createHost();
            hostList.add(host);
        }
        return new DatacenterSimple(sim, hostList);
    }

    /**
     * Build one physical server (Host):
     *  - 'sc.HOST_PES' CPU cores (PEs), each with 'sc.HOST_MIPS'
     *  - RAM/BW/Storage from config (MB/Mbps/MB)
     *  - Power model with max and static (idle) power in Watts
     */
    private Host createHost() {
        final List<Pe> peList = new ArrayList<>(sc.HOST_PES);
        for (int i = 0; i < sc.HOST_PES; i++) {
            peList.add(new PeSimple(sc.HOST_MIPS));
        }
        Host host = new HostSimple(sc.HOST_RAM, sc.HOST_BW, sc.HOST_STORAGE, peList);
        host.setPowerModel(new PowerModelHostSimple(sc.MAX_POWER_W, sc.STATIC_POWER_W));
        return host;
    }

    /**
     * Create 'sc.VMS' virtual machines.
     * Note: VM CPU uses 'sc.HOST_MIPS' per core here; if you have a separate VM_MIPS,
     * swap it in to decouple VM speed from host core speed.
     */
    public List<Vm> createVms() {
        final var vmList = new ArrayList<Vm>(sc.VMS);
        for (int i = 0; i < sc.VMS; i++) {
            final var vm = new VmSimple(sc.HOST_MIPS, sc.VM_PES);
            vm.setRam(1024).setBw(2000).setSize(20_000); // MB, Mbps, MB
            vmList.add(vm);
        }
        return vmList;
    }

    /**
     * Build cloudlets in four urgency classes (A..D) with different average CPU use.
     * Then sort so higher-priority jobs are submitted first.
     */
    public List<Cloudlet> createCloudlets() {
        var list = new ArrayList<Cloudlet>();
        list.addAll(makeClassCloudlets(sc.A_COUNT, 4, 0.80)); // A: highest priority (~80% CPU)
        list.addAll(makeClassCloudlets(sc.B_COUNT, 3, 0.60));
        list.addAll(makeClassCloudlets(sc.C_COUNT, 2, 0.40));
        list.addAll(makeClassCloudlets(sc.D_COUNT, 1, 0.20)); // D: lowest priority (~20% CPU)
        list.sort((c1, c2) -> Integer.compare(c2.getPriority(), c1.getPriority()));
        return list;
    }

    /**
     * Helper to create 'count' cloudlets of the same class.
     * - CLOUDLET_LENGTH: total work (in MI)
     * - CLOUDLET_PES:    number of cores the job needs concurrently
     * - UtilizationModelDynamic(utilFrac): average CPU fraction over time (not a hard cap)
     */
  //  private List<Cloudlet> makeClassCloudlets(int count, int priority, double utilFrac) {
  //      var list = new ArrayList<Cloudlet>(count);
  //      var util = new UtilizationModelDynamic(utilFrac);
  //      for (int i = 0; i < count; i++) {
  //          var c = new CloudletSimple(sc.CLOUDLET_LENGTH, sc.CLOUDLET_PES, util);
  //          c.setSizes(1024);        // simple placeholder for I/O sizes (MB)
  //          c.setPriority(priority); // 4 = highest ... 1 = lowest
   //         list.add(c);
   //     }
   //     return list;
  //  }
    
    
    private List<Cloudlet> makeClassCloudlets(int count, int priority, double cpuFrac) {
    var list = new ArrayList<Cloudlet>(count);

    // CPU (fraction of VM MIPS) – keep reasonably high; time-sharing handles this
    var utilCpu = new UtilizationModelDynamic(cpuFrac);

    // RAM/BW fractions MUST be small so that (concurrency_per_vm * frac) < 1.0
    // With 10k cloudlets / 500 VMs ≈ 20 cloudlets per VM on average.
    // Target sum ~ 0.6 → use ≈ 0.6 / 20 = 0.03 per cloudlet (3%).
    double ramFrac = switch (priority) {
        case 4 -> 0.030;   // critical
        case 3 -> 0.025;   // high
        case 2 -> 0.020;   // medium
        default -> 0.015;  // low
    };
    double bwFrac = ramFrac; // mirror RAM for simplicity

    var utilRam = new UtilizationModelDynamic(ramFrac);
    var utilBw  = new UtilizationModelDynamic(bwFrac);

    for (int i = 0; i < count; i++) {
        var c = new CloudletSimple(sc.CLOUDLET_LENGTH, sc.CLOUDLET_PES, utilCpu);
        c.setSizes(1024);
        c.setPriority(priority);
        c.setUtilizationModelRam(utilRam);
        c.setUtilizationModelBw(utilBw);
        // (Optional) smear arrivals to avoid a giant spike:
        // c.setSubmissionDelay(i * 0.0005); // 0.5 ms apart within class
        list.add(c);
    }
    return list;
}

    
    
}
