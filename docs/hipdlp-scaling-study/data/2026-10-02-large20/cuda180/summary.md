80 completed runs across 20 instances.

| Configuration | Verified | Measured total (s) | PAR2 (s) | Common solved total (s) |
|---|---:|---:|---:|---:|
| auto (guard + CR/Ruiz/PC) | 8/20 | 2367.954 | 4523.697 | 203.697 |
| Ruiz+PC | 9/20 | 2210.761 | 4186.325 | 197.134 |
| PC-only | 9/20 | 2178.684 | 4154.405 | 166.886 |
| two-pilot selector | 9/20 | 2232.881 | 4208.689 | 214.896 |

Common verified subset: 8 instances: fome13, graph40-40, neos1, neso-2245, pds-100, pypsa-eur-240t5n, savsched1, steiner_t.

Measured total includes failed runs at their observed duration. PAR2 charges
each failure twice its solver time limit. Input I/O and backend warmup are
excluded; scaling, pilot work and the solver final KKT check are included.
