80 completed runs across 20 instances.

| Configuration | Verified | Measured total (s) | PAR2 (s) | Common solved total (s) |
|---|---:|---:|---:|---:|
| auto (guard + CR/Ruiz/PC) | 2/20 | 914.130 | 6513.622 | 33.622 |
| Ruiz+PC | 2/20 | 899.996 | 6512.252 | 32.252 |
| PC-only | 2/20 | 901.560 | 6511.725 | 31.725 |
| two-pilot selector | 2/20 | 1023.571 | 6527.341 | 47.341 |

Common verified subset: 2 instances: graph40-40, savsched1.

Measured total includes failed runs at their observed duration. PAR2 charges
each failure twice its solver time limit. Input I/O and backend warmup are
excluded; scaling, pilot work and the solver final KKT check are included.
