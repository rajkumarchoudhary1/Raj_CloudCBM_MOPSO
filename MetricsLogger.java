package raj.cbm.metrics;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.cloudsimplus.cloudlets.Cloudlet;
import org.cloudsimplus.datacenters.Datacenter;
import org.cloudsimplus.hosts.Host;
import org.cloudsimplus.utilizationmodels.UtilizationModel;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Utility helpers to export CloudSim+ run data to CSV and to compute coarse
 * aggregate metrics used by the experiments (energy, P95 latency, SLO%).
 *
 * <h3>What this class does</h3>
 * <ul>
 *   <li>Writes a per-cloudlet table (arrival/start/finish/response time, status).</li>
 *   <li>Writes a per-cloudlet metrics CSV (response time, an approximate energy figure,
 *       and a placeholder SLO flag).</li>
 *   <li>Computes a quick, order-of-magnitude summary:
 *       <em>energy kWh (approx)</em>, <em>P95 latency (ms)</em>, and <em>SLO% (placeholder)</em>.</li>
 *   <li>Manages CSVs for run summaries and Pareto archives.</li>
 * </ul>
 *
 * <h4>Important approximations</h4>
 * <ul>
 *   <li><b>Energy</b> is computed via a crude average of idle and max power per host,
 *       multiplied by the simulation makespan. For publication-grade accuracy, integrate
 *       the power model over time instead.</li>
 *   <li><b>SLO miss</b> at per-cloudlet level is a placeholder here because actual
 *       deadlines are class-dependent; attach class/deadline metadata to each cloudlet
 *       to compute true SLOs.</li>
 * </ul>
 *
 * @author rchoudhary
 */
public class MetricsLogger {

    /** Deadline map type: Cloudlet ID -> deadline (ms). */
    public static class DeadlineMaps {
        public final java.util.Map<Long,Integer> byId;
        public DeadlineMaps(java.util.Map<Long,Integer> byId){ this.byId = byId; }
    }


    /**
     * Small carrier for aggregate metrics of a run.
     * <ul>
     *   <li>{@code energyKwh} — total energy consumed (kWh, approximate).</li>
     *   <li>{@code p95Ms} — 95th percentile of response time across finished cloudlets (ms).</li>
     *   <li>{@code sloPct} — percentage of SLO violations (placeholder here; see notes).</li>
     * </ul>
     */
    public static class Summary {
        /** Total energy in kWh (approximate). */
        public double energyKwh;
        /** 95th percentile of per-cloudlet response time in milliseconds. */
        public double p95Ms;
        /** SLO miss rate in percent (placeholder unless deadlines are attached). */
        public double sloPct;
    }

    /**
     * Writes a per-cloudlet table to CSV with timing and status information.
     *
     * <p><b>Columns</b>:
     * <ul>
     *   <li>CloudletId</li>
     *   <li>VmId</li>
     *   <li>ArrivalTime (sim seconds)</li>
     *   <li>StartTime (sim seconds)</li>
     *   <li>FinishTime (sim seconds)</li>
     *   <li>ResponseTimeMs (milliseconds)</li>
     *   <li>Status</li>
     * </ul>
     *
     * @param list the finished cloudlets to log
     * @param path the CSV file path to create/overwrite
     * @throws RuntimeException if writing fails
     */
    public static void saveCloudletsTable(List<Cloudlet> list, String path) {
        try (var out = Files.newBufferedWriter(Paths.get(path));
             var printer = new CSVPrinter(out, CSVFormat.DEFAULT
                     .withHeader("CloudletId","VmId","ArrivalTime","StartTime","FinishTime","ResponseTimeMs","Status"))) {
            for (var c : list) {
                long vmId = c.getVm() == null ? -1 : c.getVm().getId();
                double responseMs = (c.getFinishTime() - c.getArrivalTime()) * 1000.0;
                printer.printRecord(c.getId(), vmId, c.getArrivalTime(), c.getExecStartTime(), c.getFinishTime(), responseMs, c.getStatus());
            }
        } catch (IOException e) {
            throw new RuntimeException("CSV write error", e);
        }
    }

    /**
     * Writes a per-cloudlet metrics CSV including response time, a run-level approximate
     * energy (kWh) copied onto each row, and a placeholder SLO flag.
     *
     * <p><b>Columns</b>:
     * <ul>
     *   <li>CloudletId</li>
     *   <li>VmId</li>
     *   <li>ResponseTimeMs</li>
     *   <li>EnergyKwhApprox — same value for all rows; comes from {@link #computeSummary(List, Datacenter)}</li>
     *   <li>SloMiss — currently 0 for all rows (class-deadline info not available here)</li>
     * </ul>
     *
     * @param list finished cloudlets
     * @param dc   the datacenter used (for accessing host power models)
     * @param path output CSV path
     * @throws RuntimeException if writing fails
     */
    public static void saveCloudletsMetricsCsv(List<Cloudlet> list, Datacenter dc, String path) {
        saveCloudletsMetricsCsv(list, dc, path, null);
    }

    /**
     * Deadline-aware variant: if {@code deadlineMsById} is non-null, computes true SLO misses per row
     * and in the aggregate summary.
     */
    public static void saveCloudletsMetricsCsv(List<Cloudlet> list, Datacenter dc, String path, java.util.Map<Long,Integer> deadlineMsById) {
        try (var out = Files.newBufferedWriter(Paths.get(path));
             var printer = new CSVPrinter(out, CSVFormat.DEFAULT
                     .withHeader("CloudletId","VmId","ResponseTimeMs","EnergyKwhApprox","SloMiss"))) {

            Summary s = computeSummary(list, dc);
            // For per-cloudlet, mark SLO miss if response time exceeds a simple global threshold proxy (class deadlines vary by priority in full model)
            for (var c : list) {
                long vmId = c.getVm() == null ? -1 : c.getVm().getId();
                double responseMs = (c.getFinishTime() - c.getArrivalTime()) * 1000.0;
                int sloMiss = 0;
                if (deadlineMsById != null) {
                    Integer dl = deadlineMsById.get(c.getId());
                    if (dl != null && responseMs > dl) sloMiss = 1;
                }
                printer.printRecord(c.getId(), vmId, responseMs, s.energyKwh, sloMiss);
            }
        } catch (IOException e) {
            throw new RuntimeException("CSV write error", e);
        }
    }

    /**
     * Computes a coarse aggregate summary across finished cloudlets.
     *
     * <ul>
     *   <li><b>Energy (kWh)</b>: for each host, take power at 0% and 100% utilization, average them,
     *       multiply by makespan (max finish time across cloudlets), and sum. This is rough—prefer
     *       integrating host power over time for better accuracy.</li>
     *   <li><b>P95 latency (ms)</b>: 95th percentile of cloudlet response times
     *       (finish − arrival) × 1000.</li>
     *   <li><b>SLO%</b>: placeholder = 0.0; attach class/deadline metadata to compute real values.</li>
     * </ul>
     *
     * @param list finished cloudlets
     * @param dc   datacenter (for host list and power models)
     * @return a populated {@link Summary}
     */
    public static Summary computeSummary(List<Cloudlet> list, Datacenter dc) { return computeSummary(list, dc, null); }

    /**
     * Deadline-aware variant: if {@code deadlineMsById} is non-null, computes true SLO%% as
     * (#cloudlets with responseMs > deadline) / total.
     */
    public static Summary computeSummary(List<Cloudlet> list, Datacenter dc, java.util.Map<Long,Integer> deadlineMsById) {
        Summary s = new Summary();
        // Energy: use host power model * makespan (very rough; replace by time integration if needed)
        double makespanSec = list.stream().mapToDouble(Cloudlet::getFinishTime).max().orElse(0.0);
        double energyWh = 0.0;
        for (Host h : dc.getHostList()) {
            double idle = h.getPowerModel().getPower(0);
            double max = h.getPowerModel().getPower(1);
            double avg = (idle + max) / 2.0; // crude average
            energyWh += avg * (makespanSec / 3600.0);
        }
        s.energyKwh = energyWh / 1000.0;

        // P95
        List<Double> rt = list.stream().map(c -> (c.getFinishTime() - c.getArrivalTime()) * 1000.0).sorted().collect(Collectors.toList());
        int idx = (int)Math.ceil(0.95 * rt.size()) - 1;
        s.p95Ms = rt.isEmpty() ? 0.0 : rt.get(Math.max(0, Math.min(rt.size()-1, idx)));

        int misses = 0;
        if (deadlineMsById != null && !deadlineMsById.isEmpty()) {
            for (var c : list) {
                double responseMs = (c.getFinishTime() - c.getArrivalTime()) * 1000.0;
                Integer dl = deadlineMsById.get(c.getId());
                if (dl != null && responseMs > dl) misses++;
            }
            s.sloPct = list.isEmpty() ? 0.0 : (100.0 * misses / list.size());
        } else {
            s.sloPct = 0.0; // legacy behavior
        }
        return s;
    }

    /**
     * Creates the summary CSV (with header) if it does not already exist.
     *
     * <p><b>Columns</b>: scenario, algorithm, run, seed, energy_kwh, p95_ms, slo_pct,
     * hv, igd, converge_iters, cloudlets, vm_count, archive_size</p>
     *
     * @param path path to the summary CSV
     * @throws IOException if the header cannot be written
     */
    public static void initSummary(String path) throws IOException {
        if (Files.exists(Paths.get(path))) return;
        try (var out = Files.newBufferedWriter(Paths.get(path));
             var printer = new CSVPrinter(out, CSVFormat.DEFAULT
                     .withHeader("scenario","algorithm","run","seed","energy_kwh","p95_ms","slo_pct","hv","igd","converge_iters","cloudlets","vm_count","archive_size"))) {
        }
    }

    /**
     * Appends a single run’s results into the summary CSV (one line per run).
     *
     * @param path path to the summary CSV (must already have a header; see {@link #initSummary(String)})
     * @param r    the run result to append
     * @throws IOException if appending fails
     */
    public static void appendSummary(String path, raj.cbm.SingleRun.Result r) throws IOException {
        try (var out = Files.newBufferedWriter(Paths.get(path), java.nio.file.StandardOpenOption.APPEND);
             var printer = new CSVPrinter(out, CSVFormat.DEFAULT)) {
            printer.printRecord(r.scenarioName, r.algorithm, r.run, r.seed, r.energyKwh, r.p95Ms, r.sloPct, r.hv, r.igd, r.convergeIters, r.cloudletCount, r.vmCount, r.archiveSize);
        }
    }

    /**
     * Writes the current Pareto archive (objective triples) to CSV.
     *
     * <p><b>Columns</b>: energy_kwh, p95_ms, slo_pct</p>
     *
     * @param path    output CSV path
     * @param archive list of objective triples (energy, p95, slo)
     * @throws RuntimeException if writing fails
     */
    public static void saveArchive(String path, List<double[]> archive) {
        try (var out = Files.newBufferedWriter(Paths.get(path));
             var printer = new CSVPrinter(out, CSVFormat.DEFAULT.withHeader("energy_kwh","p95_ms","slo_pct"))) {
            for (double[] f : archive) {
                printer.printRecord(f[0], f[1], f[2]);
            }
        } catch (IOException e) {
            throw new RuntimeException("CSV write error", e);
        }
    }
}
