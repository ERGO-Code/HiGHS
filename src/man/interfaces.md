### Language interfaces

There are HiGHS interfaces for C, C#, FORTRAN, Julia and Python in [`HiGHS/interfaces`](https://github.com/ERGO-Code/HiGHS/tree/master/src/interfaces), with example driver files in [`HiGHS/examples`](https://github.com/ERGO-Code/HiGHS/tree/master/examples). 

Julia
-----

- A Julia interface is available at https://github.com/jump-dev/HiGHS.jl.

Rust
----

- HiGHS can be used from rust through the [`highs` crate](https://crates.io/crates/highs). The rust linear programming modeler [**good_lp**](https://crates.io/crates/good_lp) supports HiGHS. 

Javascript
----------

HiGHS can be used from javascript directly inside a web browser thanks to [highs-js](https://github.com/lovasoa/highs-js). See the [demo](https://lovasoa.github.io/highs-js/) and the [npm package](https://www.npmjs.com/package/highs).

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
