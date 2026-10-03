18 completed runs across 9 instances.

| Configuration | Verified | Measured total (s) | PAR2 (s) | Common solved total (s) |
|---|---:|---:|---:|---:|
| Ruiz+PC | 9/9 | 220.811 | 220.811 | 220.811 |
| PC-only | 9/9 | 194.602 | 194.602 | 194.602 |

Common verified subset: 9 instances: fome13, graph40-40, neos1, neso-2245, pds-100, pypsa-eur-240t5n, savsched1, steiner_t, watson_1.

Measured total includes failed runs at their observed duration. PAR2 charges
each failure twice its solver time limit. Input I/O and backend warmup are
excluded; scaling, pilot work and the solver final KKT check are included.
