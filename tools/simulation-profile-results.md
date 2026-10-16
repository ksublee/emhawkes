# Simulation profiling

## Findings

The first priority is full-history copying in the R event loop. A numeric
Hawkes simulation copies a full component-history-sized block about once per
event. Eager callback context construction adds copies even when the callback
declares only `lambda_component_n`. The state/history boundary and callback
context preparation should be addressed before choosing the native boundary.

This commit adds measurement tools and diagnostic method copies. The installed
`hsim()`, arrival algorithms, public interface and defaults are unchanged.
The diagnostic copies are R code, not a native simulator or a production speedup.

## Protocol

Baseline: `ba7965384a6c9d4dbd6c5aa9059e057223790b39`. Windows 11, Intel Core Ultra
7 356H, R 4.6.0, Rcpp 1.1.2. Each timing uses a fresh R process pinned to affinity
mask 4. Workers run sequentially. Three seeded 100-event warm-ups precede each
measurement. Setup, serialization checks and RSS queries are outside timings.
Timing, allocation tracing, sampling and RNG-call probes are separate runs.

Twelve fixed fixtures cover dimensions 1, 2, 4 and 8, one or two kernel columns
per dimension, equal/unequal within-row decays, marks, parameter-only baselines,
current-state impacts, observed-count baselines, ignored dots and signed
components. Exact parameters are in [the fixtures](simulation-profile-cases.R).
All baselines are positive. Every request generates the requested
events; `size = events + 1` includes the initial row. FRPP, zero-baseline
extinction, arbitrary retained histories and native simulation are outside this study.

Ordinary timing rows use three runs. The 8-dimension/16-column reference is
bounded at 3,000 events after a single 10,000-event stress run took **101.33 s**
(100.57 s of process CPU). That stress row has one sample, not a three-run median.
An interrupted larger-reference run is excluded. There is no 100,000-event
reference simulation. Diagnostic-only paths are measured at 100,000 events.

Rprof sampling requests 1 ms with source references and GC samples enabled.
Sample percentages, rather than nominal profiler seconds, are interpreted.
Short profiles have few samples; longer 100k diagnostic profiles support the
post-copy-removal attribution. Rprofmem counts cumulative R allocation bytes;
these are not peak RAM. Process peak RSS includes loading and warm-up.

## Reference timings

| Fixture | D | Kernel columns | Events | Median seconds |
| --- | ---: | ---: | ---: | ---: |
| univariate | 1 | 1 | 10,000 | 0.48 |
| multivariate | 4 | 4 | 10,000 | 1.75 |
| multi_equal | 4 | 8 | 10,000 | 4.44 |
| multi_unequal | 4 | 8 | 10,000 | 4.56 |
| large_equal | 8 | 16 | 3,000 | 1.16 |
| large_unequal | 8 | 16 | 3,000 | 1.18 |
| marked | 4 | 8 | 10,000 | 13.61 |
| parameter_mu | 4 | 8 | 10,000 | 4.51 |
| state_current | 4 | 8 | 10,000 | 13.45 |
| state_dots | 4 | 8 | 10,000 | 13.40 |
| state_observed | 4 | 8 | 10,000 | 13.86 |
| signed | 2 | 4 | 10,000 | 1.25 |

The parameter-only baseline resolves once in `setting()`; it does not add an
event-level baseline callback. Its generated path matches the numeric baseline.
The two current-state impact fixtures use the same model and seed; adding ignored
dots also leaves their paths unchanged.

## Allocation and source attribution

At D8/K16 and 3,000 events, the returned object is **6.71 MiB**, while cumulative
R allocation is **8.66 GiB**. A component-history-sized allocation is
`8 * (events + 1) * D * K + 48` bytes on this runtime: **3,073,072 bytes**,
observed **3,003 times**. At 1,000 events the corresponding block occurs 1,003
times and cumulative allocation is about 0.98 GiB. Block size and block count
both grow with events, giving quadratic allocation growth.

These block counts include initial/final allocations. They are shape-based
evidence rather than object identities. Source sampling and controlled input
changes corroborate the copying mechanism:

- [Arrival call](https://github.com/ksublee/emhawkes/blob/ba7965384a6c9d4dbd6c5aa9059e057223790b39/R/hsim.R#L273) passes the whole post-event history although the helper reads only its previous row.
- [Post-event history write](https://github.com/ksublee/emhawkes/blob/ba7965384a6c9d4dbd6c5aa9059e057223790b39/R/hsim.R#L354) triggers the large copy after that sharing.
- `context_at()` builds all histories before `call_component()` selects its declared arguments. Unused histories therefore become shared too.

In the D4/K8 numeric 10k reference profile, line 354 accounts for **58.2%**
of self samples. GC accounts for **26.4%**, while arrival generation occupies
only **7.3%** of the disjoint phase samples. The current-state and ignored-dots
fixtures have identical baseline allocation: declaring fewer arguments alone
does not fix eager context construction in the current implementation.

## Controlled R diagnostics

The probes evaluate a private copy of the method body. They do not register a
replacement S4 method or alter the installed namespace:

- `arrival_row`: pass only the previous history row to the Hawkes arrival helper. Its private row index is 2; public callback indices are unchanged.
- `selective_context`: evaluate only the values selected by the existing callback plan. Dots still receive the full available context.
- `both`: combine these two changes.

The probes bypass one S4 dispatch per path; they are mechanism checks rather
than released-package timings. Output and final RNG equality are verified.

| Fixture | Events | Reference | Arrival row | Selective context | Both |
| --- | ---: | ---: | ---: | ---: | ---: |
| multi_equal | 10,000 | 4.44 s | 0.38 s | 4.28 s | 0.38 s |
| large_equal | 3,000 | 1.16 s | 0.16 s | 1.10 s | 0.15 s |
| state_current | 10,000 | 13.45 s | 8.57 s | 4.64 s | 0.48 s |
| state_dots | 10,000 | 13.40 s | 8.36 s | 14.17 s | 8.92 s |
| state_observed | 10,000 | 13.86 s | 8.84 s | 6.31 s | 1.14 s |

Cumulative allocation and peak process memory from separate workers:

| Fixture | Events | Allocated GiB: reference / both | Peak RSS MiB: reference / both |
| --- | ---: | ---: | ---: |
| multi_equal | 10,000 | 23.903 / 0.064 | 163.4 / 141.8 |
| large_equal | 3,000 | 8.660 / 0.081 | 165.0 / 158.4 |
| state_current | 10,000 | 58.193 / 0.078 | 157.5 / 140.9 |
| state_dots | 10,000 | 58.193 / 34.354 | 157.5 / 166.6 |
| state_observed | 10,000 | 58.942 / 3.805 | 158.0 / 161.8 |

Peak RSS does not consistently decrease. The large change is repeated allocation
and copying, not the size of the returned history. Dots keep full histories
available, and observed-baseline/marked callbacks still receive the requested
full count/mark vectors. Their remaining quadratic costs cannot be removed by
silently replacing those public arguments with scalars or rows.

The numeric 10k stress path takes 101.33 s once, while the matching diagnostic
`both` path takes 0.50 s by three-run median. Its output and final RNG state
match. The reference stress sample is insufficient for a stable speedup estimate.

At 100k events, the diagnostic D8/K16 numeric path takes **9.38 s** and the
D4/K8 current-state path **9.15 s** by three-run medians. Returned objects are
222.81 and 64.10 MiB. No 100k reference timing or full-reference RNG comparison
is claimed. Public prefixes of 3k and 10k events, respectively, match exactly.

## Remaining Rcpp opportunity

After both diagnostic changes, the longer D8/K16 numeric profile places about
**30.4%** of samples in arrival generation, **60.2%** in the event loop and other
work, **4.6%** in validation and **4.8%** in GC. The current-state profile places
about **30.9%** in arrival generation, **44.6%** in the event loop and other work,
**15.9%** in callback dispatch, **4.7%** in validation and **4.0%** in GC.

This supports moving the numeric event loop as a unit after removing copies,
rather than accelerating only the arrival helper. Native speedup is not measured.
The unequal-decay probes reject about 23.3% of proposals at D4/K8 and 19.6% at
D8/K16 in these 1k fixtures. Distribution-call counts record requested variates,
not internal uniform RNG draws.

## Validation

All **202 completed samples** match the exact published design. Within each
model/size, simulated data fields (excluding the stored model object) and final
`.Random.seed` match across modes and diagnostic variants. This includes three
timing rounds, separate memory runs, source profiles and traced RNG-call probes.
At 100k only diagnostic variants exist, so identities there are between those
diagnostic runs.

All 12 fixtures also pass an independent 1k replay check: public and diagnostic
data/RNG are identical; counts match independent cumulative sums; both the R
and C++ likelihood/inference engines reproduce pre/post intensities and kernel
components. The maximum scaled replay discrepancy and per-event likelihood
difference in this check are **zero**. Inference does not change the RNG state.
The two larger diagnostic paths pass the public-prefix comparisons described
above. All simulator source files match the baseline commit.

No production runtime or public example changes are included. Timing figures
are host-specific; this is not a simulation distribution or coverage study.

## Next implementation

1. Separate the current post-event component matrix from output histories; pass it to the Hawkes arrival calculation.
2. Reuse the existing callback-plan/prepared-call mechanism so unrequested histories are not materialized. Preserve full arguments requested explicitly or through dots.
3. Verify complete output layout, seeds, mark/callback order, retained snapshots, signed components, empty paths and extinction before shipping. Avoid unnecessary final full-row copies when no trimming is needed, subject to the same checks.
4. Re-profile the production R path. Then start a native event loop with numeric parameters, positive baselines and equal within-row decays; expand thinning and callbacks in later changes.

## Reproduce

Install the same source into a normal library and a profiling library using
`R CMD INSTALL --with-keep.source`. From the repository root:

```sh
Rscript --vanilla tools/run-simulation-profile.R plan LIB PROFILE_LIB RAW R/hsim.R
Rscript --vanilla tools/run-simulation-profile.R run LIB PROFILE_LIB RAW R/hsim.R
Rscript --vanilla tools/check-simulation-profile.R RAW PUBLISHED
Rscript --vanilla tools/check-simulation-replay.R LIB PUBLISHED R/hsim.R
```

Set `EMHAWKES_BENCH_AFFINITY=4` on Windows for the recorded affinity. The optional
`stress` argument adds the single expensive 10k D8/K16 reference run. Keep raw
profiles, RDS files and immutable worker copies outside the repository. The
runner validates existing worker-source and installed-code hashes before resuming.

See [design](simulation-profile-results/design.csv),
[timing samples](simulation-profile-results/timings.csv),
[timing medians/ranges](simulation-profile-results/timing-summary.csv),
[memory](simulation-profile-results/memory.csv),
[history-sized allocations](simulation-profile-results/history-sized-allocations.csv),
[allocation attribution](simulation-profile-results/allocation-attribution.csv),
[profile phases](simulation-profile-results/profile-phases.csv),
[RNG distribution calls](simulation-profile-results/rng-calls.csv),
[identity checks](simulation-profile-results/identity-checks.csv),
[independent replay](simulation-profile-results/replay-checks.csv),
[large prefixes](simulation-profile-results/large-prefix-checks.csv), and
[session information](simulation-profile-results/session-info.txt),
[manifest](simulation-profile-results/manifest.txt), and
[source verification](simulation-profile-results/source-verification.txt).
