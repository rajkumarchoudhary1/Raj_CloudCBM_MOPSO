Cloud-CBM (Raj_CloudCBM_MOPSO)
================================

Repository
----------
GitHub: https://github.com/rajkumarchoudhary1/Raj_CloudCBM_MOPSO

Clone:
  git clone https://github.com/rajkumarchoudhary1/Raj_CloudCBM_MOPSO
  cd Raj_CloudCBM_MOPSO

What is this code?
------------------
This repository contains a reproducible research package for Cloud-based Condition-Based Maintenance (CBM) scheduling.
It models CBM workloads in CloudSim Plus (v8.0.0) and compares a domain-adapted, deadline-/priority-aware Multi-Objective
Particle Swarm Optimizer (MOPSO) against established baselines (NSGA-II, SPEA-II, MOGA, MOACO, and a Round-Robin control).

Why we wrote this code
----------------------
Aircraft CBM analytics must be processed quickly and efficiently in the cloud.
We formulate scheduling as a tri-objective optimization problem:
  • Minimize datacenter ENERGY (kWh)
  • Minimize tail latency P95 (ms)
  • Minimize SLO-MISS rate (%)
The code provides a fully reproducible pipeline—fixed seeds, pinned dependencies, and standard outputs—so results can be
validated, compared, and extended with new algorithms or workloads.

What’s included
---------------
- Java 17 project (Maven) pinned to CloudSim Plus 8.0.0
- A common experiment harness that runs (workload × algorithm × runs) and writes deterministic CSVs
- Algorithm adapters:
    * MOPSO (proposed, deadline-/priority-aware)
    * NSGA-II, SPEA-II, MOGA (genetic/strength Pareto baselines)
    * MOACO (multi-objective ant colony)
    * DefaultRoundRobin (non-optimizing control)
- Configuration files:
    * config/scenarios.yaml  — workloads, VM counts, power model, deadlines, repetitions
    * config/seeds.txt       — 30 deterministic seeds per workload
- Analysis scripts (optional, outside Java): generate Pareto plots, confidence-interval bars, and ANOVA tables
- Results directory schema with per-run artifacts (cloudlets.csv, archive.csv, optional power traces)

Typical repository layout
-------------------------
.
├─ pom.xml
├─ src/main/java/raj/cbm/
│  ├─ Main.java                 (orchestrator)
│  ├─ ScenarioConfig.java       (YAML loader)
│  ├─ SimulationConfig.java     (derive simulator knobs)
│  ├─ CloudResourceFactory.java (Datacenter/Host/VM/Cloudlet builders)
│  ├─ SimulationRunner.java     (execute mapping in CloudSim Plus)
│  ├─ metrics/MetricsLogger.java(CSV schemas & I/O)
│  └─ alg/
│     ├─ Algorithm.java         (common interface: optimize(prefix) → Result)
│     ├─ MOPSO.java             (proposed)
│     ├─ NSGAII.java            (baseline)
│     ├─ SPEAII.java            (baseline)
│     ├─ MOGA.java              (baseline)
│     ├─ MOACO.java             (baseline)
│     └─ DefaultRoundRobin.java (control)
├─ config/
│  ├─ scenarios.yaml
│  └─ seeds.txt
├─ results/                     (auto-created at runtime)
└─ analysis/                    (optional R scripts for plots & stats)

Prerequisites
-------------
- Java 17 (LTS)
- Maven 3.8+
- Internet access on first build to fetch dependencies
- Optional: R (v4+) with tidyverse/ggplot2 for analysis scripts

How to build and run
--------------------
1) Build the shaded JAR:
   mvn -q clean package

2) Execute the experiment matrix:

Rs-MacBook-Pro:Raj_CloudCBM_MOPSO rchoudhary$ mvn -U -DskipTests package

Rs-MacBook-Pro:Raj_CloudCBM_MOPSO rchoudhary$ mvn -q -Dexec.mainClass=raj.cbm.Main exec:java -- \
>   --scenarios config/scenarios.yaml \
>   --seeds     config/seeds.txt \
>   --out       results

Expected outputs:
- results/summary_all.csv                (540 rows for full matrix: 3 workloads × 6 algorithms × 30 runs)
- results/summary_low.csv                (by scenario)
- results/summary_medium.csv
- results/summary_high.csv
- results/<scenario>/<algorithm>/run_<r>/
    cloudlets.csv, cloudlets_metrics.csv, archive.csv, (optional) power_trace.csv

Configuration
-------------
config/scenarios.yaml controls workloads and environment. Example keys:
  workloads: [{{name: low, cloudlets: 100, vms: 50}}, {{name: medium, cloudlets: 1000, vms: 150}}, {{name: high, cloudlets: 10000, vms: 500}}]
  priorities: {{critical: 0.20, high: 0.30, medium: 0.30, low: 0.20}}
  deadlines_ms: {{critical: 300, high: 600, medium: 900, low: 1200}}
  power_model: {{max_power_w: 200.0, static_power_w: 120.0}}
  runs_per_combo: 30

config/seeds.txt has 3 lines (one per workload), each with 30 integers (one seed per run). The harness maps (scenario, run) → seed.

What each output file contains
------------------------------
summary_*.csv columns:
  scenario,algorithm,run,seed,energy_kwh,p95_ms,slo_pct,hv,igd,converge_iters,cloudlets,vm_count,archive_size

Per-run artifacts (results/<scenario>/<algorithm>/run_<r>/):
  cloudlets.csv          id, priority, vm, length_mi, start, finish, response_ms, deadline_ms, slo_miss
  cloudlets_metrics.csv  per-cloudlet metrics + run-level energy_kwh for quick joins
  archive.csv            Pareto points (energy_kwh, p95_ms, slo_pct)
  power_trace.csv        OPTIONAL: time_s,host_id,util,power_w

Algorithms & default parameters
-------------------------------
MOPSO (proposed): w=0.7, c1=c2=2.0, N=30, T=50, archive=100
NSGA-II / MOGA:   N=60, G=50, SBX pc=0.9, poly mutation pm=0.02
SPEA-II:          N=60, archive=60, G=50
MOACO:            ants=30, ρ=0.1, α=1, β=2, iterations=50
DefaultRoundRobin: no parameters (simple cyclic assignment)

How to add a new algorithm
--------------------------
1) Implement raj.cbm.alg.Algorithm with method:
   Result optimize(String outputPrefix);
   returning: mapping:int[], paretoFront:double[][3], hv:double, igd:double, converge_iters:int

2) Create src/main/java/raj/cbm/alg/MyAlgo.java and register it in Main.java where algorithms are enumerated.

3) Rebuild & run the matrix:
   mvn -q clean package
   java -jar target/cloudcbm-1.0.0-shaded.jar --scenarios config/scenarios.yaml --seeds config/seeds.txt --out results/

Reproducibility checklist
-------------------------
[ ] Java 17; Maven build OK
[ ] pom.xml pins CloudSim Plus 8.0.0
[ ] config/scenarios.yaml and config/seeds.txt present
[ ] static_power_w < max_power_w (guardrail enforced)
[ ] results/summary_all.csv has expected rows (full matrix: 540)
[ ] analysis scripts read results/ and regenerate figures/tables

Troubleshooting
---------------
- Build fails: confirm Java 17 and Maven are installed; run `mvn -q dependency:tree` to verify CloudSim Plus 8.0.0.
- Guardrail abort: if static_power_w >= max_power_w, fix values in config/scenarios.yaml.
- Very slow runs: reduce runs_per_combo in scenarios.yaml or test only the 'low' workload during development.
- Empty figures: ensure results/summary_all.csv exists before running the R scripts in analysis/.

License / Citation / Contact
----------------------------
- License: add MIT/Apache-2.0 (or as you prefer) to the repo.
- Cite: your thesis/paper if used in publications.
- Maintainer: RajKumar Choudhary (GitHub: https://github.com/rajkumarchoudhary1)
