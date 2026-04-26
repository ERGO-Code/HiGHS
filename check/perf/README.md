# HiGHS hot-kernel microbenchmarks

Opt-in microbenchmarks for the kernels that dominate simplex/MIP runtime:

- `HFactor::build` — basis refactorisation
- `HFactor::ftranCall` / `btranCall` at hyper-sparse / sparse / dense input densities
- `HighsSparseMatrix::product` / `productTranspose` (constraint-matrix SpMV)
- `HighsSparseMatrix::priceByColumn` / `priceByRow` (the PRICE step)
- `HVectorBase<double>::saxpy` at three input densities

Output is a CSV (`min_ns`, `median_ns`, `mean_ns`, `stddev_ns`, `reps`) plus a
human-readable stdout summary, designed for before/after diffing of CPU
optimisation work.

## Build

```sh
cmake -S . -B build -DHIGHS_PERF=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j --target highs_perf
```

`HIGHS_PERF=OFF` by default; existing builds and CI are unaffected.

## Run

```sh
# Full baseline across the bundled instance set.
./build/check/perf/highs_perf

# Restrict instances or kernels.
./build/check/perf/highs_perf --instances=adlittle.mps,e226.mps --kernels=factor

# Pin to a core for stable timings, optionally under perf stat.
taskset -c 2 ./build/check/perf/highs_perf --instances=cplex1.mps
perf stat -e cache-misses,cache-references,instructions,cycles \
    taskset -c 2 ./build/check/perf/highs_perf \
        --instances=cplex1.mps --kernels=factor
```

## Notes on the numbers

- **Min** is the cleanest signal of best-case CPU performance; **median** is
  noise-robust; **stddev** flags timing instability (thermal, scheduling).
- Each kernel runs for at least `min_wall_ns` (default 100 ms) to amortise
  clock noise, capped by `max_reps`.
- Fixture setup (read MPS, run simplex, build initial factor) is done once per
  instance and excluded from kernel timings.
- For `ftranCall` / `btranCall`, the input HVector is restored after each
  timed call; the restore is outside the `steady_clock` bracket so it doesn't
  inflate the kernel timing.
- `HFactor::build` is benched on a *fresh* `HFactor` per repetition because
  `build()` mutates internal state.

## Adding a kernel

1. Add the timed function in `perf_<area>.cpp` and a header.
2. Wire it from `perf_main.cpp::main` and the kernel-name dispatch.
3. Add it to `CMakeLists.txt`.
