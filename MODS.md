# Modifications for v1.7.1

* Python build update
   * Setuptools from v1.7.0 only worked on Python 3.12
   * Drop setuptools, use scikit-build with cmake

# Modifications for v1.7.0 

* Introduces exciting new first order LP solver (cuPDLP-c)
   * Currently only serial CPU
   * GPU WIP
* The HiGHS interior point solver now allows analytic centre calculations to be performed
* Python interface update
   * Drop meson, go back to setuptools and cmake

# Modifications for v1.3.0 since v1.2.2 

* Highs::setSolution can now be used to give a solution to the simplex solver (#775)
* Highs::addVar; Highs::addVars; Highs::deleteVars(Interval/set/mask) introduced for more natural modelling
* logHeader now written as first output, even when using libraries (#784)
* Highs::presolve and Highs::postsolve completed
* Highs::resetGlobalScheduler added to reset the global scheduler
* write_solution_style option modified: 
   * "Old" HiGHS raw solution style now indicated by value `kSolutionStyleOldRaw = -1`
   * Raw and pretty solution formats for Glpsol now indicated by values `kSolutionStyleGlpsolRaw = 2` and `kSolutionStyleGlpsolPretty = 3`
* Many minor fixes handling marginal LP fle behaviour
* Highs::crossover completed (#815)
* scaled_model_status_ removed from Highs (#814)
* Major revisions of CMake build system


# Planned modifications beyond v1.3.0

* Make use of HFactor in critical parts of IPX


