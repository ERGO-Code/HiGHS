### Language interfaces

There are HiGHS interfaces for C, C#, FORTRAN, Julia and Python in [`HiGHS/interfaces`](https://github.com/ERGO-Code/HiGHS/tree/master/src/interfaces), with example driver files in [`HiGHS/examples`](https://github.com/ERGO-Code/HiGHS/tree/master/examples). 

Julia
-----

- A Julia interface is available at [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl).

Rust
----

- HiGHS can be used from rust through the [`highs` crate](https://crates.io/crates/highs). The rust linear programming modeler [**good_lp**](https://crates.io/crates/good_lp) supports HiGHS. 

Javascript
----------

HiGHS can be used from javascript directly inside a web browser thanks to [highs-js](https://github.com/lovasoa/highs-js). See the [demo](https://lovasoa.github.io/highs-js/) and the [npm package](https://www.npmjs.com/package/highs).

GAMS
----

- A GAMS interface is available at [GAMSlinks](https://github.com/coin-or/GAMSlinks/), including [pre-build libraries](https://github.com/coin-or/GAMSlinks/releases).

OSI
---

* `OSI_ROOT` :

    path to COIN-OR/Osi build/install folder (OSI_ROOT/lib/pkg-config/osi.pc should exist)
