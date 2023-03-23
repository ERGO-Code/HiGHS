# [Further features](@id guide-further)

### [Hot start](@id hot-start)

It may be possible for HiGHS to start solving a model using data
obtained by solving a related model, or supplied by a user. Whether
this is possible depends on the the class of model being solved, the
solver to be used, and the modifications (if any) that have been to
the incumbent model since it was last solved.

#### LP

To run HiGHS from a user-defined solution or basis, this is passed to HiGHS
using the methods [setSolution](@ref Set-solution) or [setBasis](@ref Set-basis).
