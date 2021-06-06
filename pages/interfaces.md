---
title: Interfaces
permalink: /interfaces/
---
### Language interfaces

There are HiGHS interfaces for C, C#, FORTRAN, Julia and Python in `HiGHS/src/interfaces`, with example driver files in `HiGHS/examples`. 

Further documentation will be published after the refactoring currently in progress is complete.

#### GAMS

Set custom options with `-D<option>=<value>` during the configuration step ( `cmake ..` ):

* `GAMS_ROOT` :

    path to GAMS system: enables building of GAMS interface

If build with GAMS interface, then HiGHS can be made available as solver
in GAMS by adding an entry for HiGHS to the file gmscmpun.txt in the GAMS
system folder (gmscmpnt.txt on Windows):

``` 
HIGHS 11 5 0001020304 1 0 2 LP RMIP
gmsgenus.run
gmsgenux.out
/path/to/libhighs.so his 1 1
```

#### OSI

* `OSI_ROOT` :

    path to COIN-OR/Osi build/install folder (OSI_ROOT/lib/pkg-config/osi.pc should exist)