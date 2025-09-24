package raj.cbm.alg;

import raj.cbm.SimulationConfig;
import java.util.*;

/**
 * DefaultAlgo
 * -----------
 * Code summary:
 * Baseline scheduler with **no multi-objective optimization**.
 * It assigns jobs to VMs in a simple **round-robin** pattern:
 * cloudlet 0 → VM 0, cloudlet 1 → VM 1, ... wrapping around at vmCount.
 *
 * Why it exists:
 * - Provides a simple, deterministic baseline to compare against MOPSO/MOEA methods.
 * - Useful to show the benefit of intelligent scheduling on energy, P95, and SLO%.
 *
 * Outputs (placeholders for consistency with other algos):
 * - paretoFront: a single dummy point
 * - hv / igd   : dummy indicators (since no real Pareto search is performed)
 *
 * @author rchoudhary
 */
public class DefaultAlgo implements Algorithm {
    private final SimulationConfig sc;
    private final Random rnd; // Not used for mapping; kept for API symmetry and potential extensions.

    public DefaultAlgo(SimulationConfig sc, long seed) {
        this.sc = sc;
        this.rnd = new Random(seed ^ 0x9E3779B97F4A7C15L);
    }

    @Override
    public Result optimize(String outPrefix) {
        Result r = new Result();

        // Build a round-robin mapping: i-th cloudlet -> (i mod vmCount)
        int n = sc.totalCloudlets();
        r.mapping = new int[n];
        int vmCount = Math.max(1, sc.VMS);
        for (int i = 0; i < n; i++) {
            r.mapping[i] = i % vmCount;
        }

        // Dummy outputs for compatibility with reporting/plots.
        r.paretoFront = List.of(new double[]{10.0, 250.0, 8.0}); // {energy, p95, slo%} example
        r.hv = 0.5;
        r.igd = 0.1;
        r.convergeIters = 0;

        return r;
    }
}
