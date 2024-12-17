## Build changes

## Code changes

HiGHS now handles multiple linear objectives by either blending using weights, or performing lexicographic optimization: see https://ergo-code.github.io/HiGHS/stable/guide/further/#guide-multi-objective-optimization

Fixed minor bug in bound checking in presolve

Fixed bug in `floor(HighsCDouble x)` and `ceil(HighsCDouble x)` when argument is small

Added some sanity checks to Highs::writeLocalModel to prevent segfaults if called directly by a user



