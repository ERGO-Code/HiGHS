include(CTest)

# define the set of feasible instances
set(successInstances
    "brandy\;0.1518509896488e4"
    )

set(failInstances
    )

# define settings
set(settings
    " "
    )

# define a macro to add tests
#
# add_instancetests takes an instance group and a status
# that the solver should report as arguments
macro(add_instancetests instances solutionstatus)
# loop over the instances
foreach(instance ${${instances}})
    # loop over all settings
    foreach(setting ${settings})
        # treat the instance as a tuple (list) of two values
        list(GET instance 0 name)
        list(GET instance 1 optval)
        # specify the instance and the settings load command
        set(inst "${PROJECT_SOURCE_DIR}/check/${name}.mps")
        # add a test and test properties
        add_test(NAME ${name}${setting} COMMAND $<TARGET_FILE:h2gmp> ${setting} -f ${inst})
        set_tests_properties (${name}${setting} PROPERTIES PASS_REGULAR_EXPRESSION "Solve plain: ${solutionstatus}")
    endforeach(setting)
endforeach(instance)
endmacro(add_instancetests)

# add tests for success and fail instances
add_instancetests( successInstances "OPTIMAL" )
add_instancetests( failInstances "NOT-OPT" )
