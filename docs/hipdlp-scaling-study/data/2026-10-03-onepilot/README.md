# One-pilot study evidence

`training.json` freezes the name-only family partition, label provenance,
thresholds and tree. `model.txt` is its study-driver serialization.
`static-ablation.json` confirms that omitting pilot progress does not change
this tree's validation choices. Runtime seconds are predicted from verified
CPU iteration labels and live pilot throughput; timeout times are not labels.

`features/` contains 200-iteration P observations for the 89 verified label
triples. `cpu-validation/` contains the 13 validation names **plus stair as a
training-sample diagnostic**, at 30 million iterations / 120 seconds.
`cpu-family-validation-summary.json` excludes stair; `stair-diagnostic.json`
contains it separately. `cpu-validation-summary.json` is the combined batch.
`cpu-heldout/` is the separate 26-name batch at 500,000 iterations / 60 seconds.
CPU elapsed times in raw records are not policy comparisons; use iteration
summaries and their explicitly censored names.

Each batch retains commands, driver/model hashes, full result records,
compressed complete logs, immutable source/payload hashes and post-run source
integrity checks. Benchmark payloads and executables remain in `/tmp`.
`build-manifest.json`, `hardware.txt` and validation logs document the frozen
build and pre-existing full-suite failure parity.


`cuda-common/` compares P/L/A/M serially at 180 seconds on nine known large
solves. `cuda-common-repeat/` repeats P/M; both batches verify all nine under
P/M and the tree proposes P on all of them. `cuda-heldout-pilot/` compares P/F
at 200 iterations on four independent large heldouts and measures observation
cost only. `cuda-extended/` runs P/M with 240 seconds on two previous timeouts;
cont11 is now verified, neso-2005 remains censored, and no P/L/A winning label
is assigned without longer L/A measurements.

`cuda-feature-coverage.json` reconstructs proposals and flags unused setup
estimates: the frozen driver fills that diagnostic for P leaves too, where it
does not affect the decision and does not represent extra setup.
