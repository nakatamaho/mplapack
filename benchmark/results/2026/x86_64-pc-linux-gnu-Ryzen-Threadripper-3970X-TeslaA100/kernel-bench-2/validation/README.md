# Validation and execution receipts

These are historical receipts, not a claim that every preliminary attempt
passed. Early build failures, interrupted attempts, cold-launch investigation,
and a corrected postprocessing-script syntax error are deliberately retained.
No numerical production correction was needed.

Final acceptance:

- `item4-review.log.exit`:0; final autotools regeneration/build and CMake builds,
  six policy tests, CPU CLI checks for four backends in both builds, and the
  last304 CPU/GPU operand comparisons. Individual `item4-review-*.log.command`,
  `.log`, `.log.exit` files contain commands and results.
- `item4-gpu-*.log`: required sequential focused GPU tests, then inverse512
  memcheck.12 tests passed; memcheck reported0 errors. The before/after receipts
  record shared-device conditions. `item4-gpu-acceptance.log.exit`:0.
- `item4-cmake-bf-host.log` and `item4-cmake-mpfr-host.log`:9/9 and6/6 focused
  host tests. The corresponding exit receipts are0.
- `item4-symbols-reviewed.log`:16 reused CPU/CUDA archives, exit0; exact archive
  arguments are in its `.command` receipt.
- `item4-architectures-reviewed.log`: sm80/sm90 cubins in all eight GPU
  benchmark executables. Only sm80 hardware was executed.
- `item4-plot-regression-packaged.log`: actual historical MPFR512 CSV redraws
  byte-identical with the previous/current plotting defaults, light and dark.
- `item4-summary-packaged.log`: final derived coverage, ratios and accounting.
- `item4-plots-caption-reviewed.log.exit`:0, all76 graphs generated.
- `item4-sweep-primed.log.exit` and `item4-extensions.log.exit`:0. Performance
  finished before final review tests began. The initial CPU interruption has
  no recoverable parent exit status; that limitation is not converted to PASS.

The six sibling `checks*/checks.csv` files each record152 comparisons; all912
are component-identical. Their bundled `runs.log` files preserve all1824
CPU/GPU command/output/exit receipts. Temporary operand dumps were compared
then removed by the test harness, not retained as performance data.

`../raw/` retains8534 point command/output/exit receipts, including all excluded
and superseded measurements. Bundling did not change `../kernels.csv`; its
pre-bundle SHA256 is retained here. Remaining partial interruption logs are
kept separately. The source snapshot is identified by the prerequisite commits
and the benchmark commit containing this directory.

Reproduce derived tables and the legacy-render check from the repository root:

```sh
D=benchmark/results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/kernel-bench-2
python3 "$D/validation/derive_tables.py"
PYTHONPATH=/home/docker/plotlib python3 "$D/validation/legacy_plot_check.py"
```

These two helpers do not execute GPU or numerical benchmark programs.
