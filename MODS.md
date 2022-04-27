# Modifications for v1.3.0 since v1.2.2 

* Added (partial) Python wrapper `highspy`

* Highs::setSolution can now be used to give a solution to the simplex solver (#775)

* Highs::addVar; Highs::addVars; Highs::deleteVars(Interval/set/mask) introduced for more natural modelling

* logHeader now written as first output, even when using libraries (#784)

* Highs::presolve and Highs::postsolve completed

* Highs::resetGlobalScheduler added to reset the global scheduler

# Planned modifications for v1.3.0

* Highs::crossover completed (#815)

* scaled_model_status_ removed from Highs (#814)

# Planned modifications beyond v1.3.0

* Make use of HFactor in critical parts of IPX


