# Event, dimension and kernel scaling

Tested runtime source: `1e9b01b54cf8d8fc89a0915ee68c931f9d20d39d`.
The accompanying changes concern measurements and documentation; the likelihood
and optimizer implementations are unchanged.

## Scope

- Public numeric `logLik()` and `logLik(infer = TRUE)`: 1,000, 10,000 and
  100,000 events; 1, 4 and 8 dimensions; one or two kernels per event type.
  The component matrix has dimensions times kernels entries, so the largest
  case has 8 rows and 16 columns (128 components). Decay rates differ by column.
- Complete public `hfit()`: simulated 4-dimensional paths with 1,000 and 10,000
  events, 4 or 8 kernel columns and three shared parameters. Both engines use
  identical starts, explicit initial components, positivity constraints, BFGS,
  numerical differentiation and the default maxLik Hessian. The fit fixture
  shares one decay parameter across all components; doubling kernel columns
  measures recursion cost without introducing another decay timescale.
- Bid-ask: 1,000, 5,000 and 10,000 events, four dimensions and four kernel
  columns. Compare the same reset rule receiving either current components or
  full component history. Measure fitting-objective preparation and its first
  evaluation separately, rather than a full optimization at every size.

The numeric likelihood fixtures use seeded exponential intervals and independent
event types. They test evaluation scaling, not recovery of a Hawkes generating
model. Full-fit and Bid-ask paths are simulated from their supplied models.
The initial zero row is excluded from the reported event count.

## Measurement

Windows 11, Intel Core Ultra 7 356H, about 32 GiB physical memory, R 4.6.0,
maxLik 1.5-2.2 and Rcpp 1.1.2.
Each timing or memory sample runs in a fresh R process pinned to logical CPU 2
(affinity mask 4). Engine order alternates over three timing rounds. Tables use
elapsed-time medians. Loading, simulation, a small warm-up and an explicit GC
are excluded. Short numeric evaluations are timed in batches; per-call elapsed
and CPU times, batch sizes and warnings are retained in the CSV files.

Memory is measured in a separate process, with one evaluation and no timing
claim. Full-fit allocation traces can produce gigabytes of temporary logs and
take substantially longer than normal fitting; logs are removed after parsing. `Rprofmem()` totals its numeric R-allocation events (not small-block reuse or
`new page` notices); this is recorded cumulative allocation,
not peak memory, and it excludes C++ standard-library buffers. The recorded
R-heap peak is the sum of the post-call GC maxima for node and vector heaps;
these maxima need not occur simultaneously. Its increment is measured above
the GC baseline immediately before evaluation. Windows peak working set records
the whole process, including libraries, loaded fixtures and warm-up. It is read
before allocation-log parsing. Neither measure is an isolated native-buffer peak.

The Bid-ask auto column calls the native fitting adapter directly; standalone
`logLik()` still selects R for function-valued models. The R reference restores
parameters and calls the R recursion with pre-normalized observations. Objective
preparation excludes the common observation normalization and parameter-map
construction. The objective is fresh for every measured evaluation, so its
one-point cache cannot turn a recorded evaluation into a cache hit.

## Results

All 464 sample files were checked: the 392 numerical-result records were
exactly identical to their R references, including memory-profiled runs; the
72 preparation-only records have no numerical result. Maximum absolute
difference was zero. All 32 full-fit samples converged with code 0 and had identical
estimates, log likelihoods, gradients, Hessians and iteration counts. No warnings
were recorded. See [checks.csv](scaling-results/checks.csv),
[fit_fields.csv](scaling-results/fit_fields.csv) and
[estimates.csv](scaling-results/estimates.csv).

### Numeric evaluation at 100,000 events

Elapsed seconds per public call; each cell is the median of three batched samples.
R uses the same installed source with `emhawkes.engine = "r"`.

| Dimensions | Kernel columns | Likelihood R | Auto | Speedup | Inference R | Auto | Speedup |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | 0.398 | 0.016 | 24.88x | 1.665 | 0.025 | 66.60x |
| 4 | 8 | 1.330 | 0.264 | 5.04x | 3.975 | 0.340 | 11.69x |
| 8 | 16 | 1.698 | 0.992 | 1.71x | 3.515 | 1.235 | 2.85x |

The dense 8-by-16 case benefits less than the univariate case. Native execution
still improves its time and allocations; these results do not establish a
uniform speedup across dimensions. Inference includes the returned intensity,
component and residual paths, whose output storage grows with data size.

### Complete estimation

Four-dimensional simulations with three shared parameters. These times include
public setup, numerical optimization and the default final Hessian.

| Events | Kernel columns | R seconds | Auto seconds | Speedup | R recorded MiB | Auto recorded MiB | R heap peak MiB | Auto heap peak MiB |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 4 | 1.03 | 0.28 | 3.68x | 85.97 | 0.81 | 119.0 | 101.3 |
| 10,000 | 4 | 9.54 | 2.17 | 4.40x | 823.48 | 7.58 | 145.3 | 106.4 |
| 1,000 | 8 | 1.08 | 0.39 | 2.77x | 443.84 | 0.92 | 135.7 | 101.5 |
| 10,000 | 8 | 20.14 | 4.84 | 4.16x | 6389.72 | 8.76 | 145.4 | 106.7 |

The large allocation reduction reflects fewer intermediate R objects; it is not
an equivalent reduction in peak RAM. For 10,000 events and eight kernel columns,
R records about 6.24 GiB cumulatively while its heap peak is only about 145 MiB.
Auto records about 8.76 MiB cumulatively and peaks at about 107 MiB.

### Bid-ask history cost

The three interfaces return exactly the same likelihood in both engines at all
three sizes. The ordinary full-history callback reads the current row and does
not retain the argument. The retained-history variant additionally saves the
argument whenever `n %% 100 == 0`, without changing its returned impact. This
retention scenario was measured in a follow-up pass.

At 10,000 events, one fresh fitting-objective evaluation gives:

| Interface | R seconds | Auto seconds | R recorded MiB | Auto recorded MiB | R heap peak MiB | Auto heap peak MiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Current components | 0.30 | 0.23 | 5.48 | 0.65 | 107.8 | 105.2 |
| Full history, no retention | 0.38 | 0.29 | 6.37 | 1.54 | 109.2 | 106.6 |
| Full history, retain every 100 indices | 0.49 | 0.37 | 128.46 | 123.63 | 227.2 | 227.3 |

For the native adapter, recorded allocation grows as follows:

| Events | Current components MiB | Full history, no retention MiB | Retain every 100 indices MiB |
| ---: | ---: | ---: | ---: |
| 1,000 | 0.38 | 0.17 | 1.39 |
| 5,000 | 0.50 | 0.78 | 31.30 |
| 10,000 | 0.65 | 1.54 | 123.63 |

In this R version, reading a full history without retaining it did not trigger
repeated full-matrix copies. Its main additional cost is storing a history with
`(events + 1) * dimensions * kernel_columns` entries: about 1.22 MiB in the
10,000-event, 4-by-4 case. Declaring a history alone does not imply quadratic
copying.

When the callback retains histories, later writes must preserve those snapshots.
The number of saved full matrices grows with events, so retaining every 100
indices yields approximately quadratic allocation and retained-storage growth.
Keeping every event would be much more expensive; that case was not run at
10,000 events. The copy is necessary for callback semantics. No unsafe in-place
mutation was introduced.

All 10, 50 and 100 retained matrices were also compared independently with public
R likelihood evaluation. R and C++ snapshots were exactly identical, future rows
remained zero and their likelihoods were identical. See
[bid-snapshots.csv](scaling-results/bid-snapshots.csv) and
[bid-variants.csv](scaling-results/bid-variants.csv).

The GitHub Markdown Bid-ask example now requests `lambda_component_n`, matching
the extended vignette. This avoids an unnecessary full history while preserving
the reset rule. Help and examples explain when history copying occurs. No
likelihood, optimizer or automatic-dispatch change was needed.

All combinations, per-sample elapsed/CPU times, recorded allocations, largest
blocks and whole-process peaks are in [summary.csv](scaling-results/summary.csv),
[timings.csv](scaling-results/timings.csv) and [memory.csv](scaling-results/memory.csv).
Three samples give a limited timing estimate. Small calls remain sensitive to
clock resolution, and CPU frequency/background work affect longer calls.
These measurements describe this machine and fixtures.

## Reproduction

Install the tested source in a private library. From the repository root, use a
new output directory and run the commands below. On Windows optionally set
`EMHAWKES_BENCH_AFFINITY=4` first; use an allowed mask on the measurement host.

```
Rscript --vanilla tools/benchmark-scaling.R prepare LIB OUT
Rscript --vanilla tools/benchmark-scaling.R run LIB OUT 3
Rscript --vanilla tools/benchmark-scaling.R compare LIB OUT
Rscript --vanilla tools/benchmark-scaling.R snapshots LIB OUT
```

The run resumes existing samples. Do not reuse an output directory after changing
the source, fixtures, worker code or CPU policy. Raw fixture and result RDS files
remain local; the report includes the CSV summaries and session information.
