include(CTest)

# define the set of feasible instances
set(successInstances
    "adlittle\;         0.2254949631e6"
    "afiro\;           -0.4647531428e3"
    "agg\;             -0.3599176728e8"
    "beaconfd\;         0.3359248580e5"
    "blend\;           -0.3081214984e2"
    "bore3d\;           0.1373080394e4"
    "brandy\;           0.1518509896e4"
    "capri\;            0.2690012913e4"
    "etamacro\;        -0.7557152333e3"
    "finnis\;           0.1727910655e6"
    "grow7\;           -0.4778781181e8"
    "israel\;          -0.8966448218e6"
    "kb2\;             -0.1749900129e4"
    "lotfi\;           -0.2526470606e2"
    "recipe\;          -0.266616e3"
    "sc105\;           -0.5220206121e2"
    "sc205\;           -0.5220206121e2"
    "sc50a\;           -0.6457507705e2"
    "sc50b\;           -0.7e2"
    "scagr25\;         -0.1475343306e8"
    "scagr7\;          -0.2331389824e7"
    "scfxm1\;           0.1841675902e5"
    "scorpion\;         0.1878124822e4"
    "scrs8\;            0.9042969538e3"
    "scsd1\;            0.8666666674e1"
    "seba\;             0.157116e5"
    "share1b\;         -0.7658931857e5"
    )

set(infeasibleInstances
    "bgetam\;        infeasible"
    "box1\;          infeasible"
    "ex72a\;         infeasible"
    "forest6\;       infeasible"
    "galenet\;       infeasible"
    "gams10am\;      infeasible"
    "klein1\;        infeasible"
    "refinery\;      infeasible"
    "woodinfe\;      infeasible"
    )

set(unboundedInstances
     "gas11\;         unbounded"
    )

set(failInstances
    )

# define settings
set(settings
    " "
    "-m"
    "-pOn"
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
        set(inst "${PROJECT_SOURCE_DIR}/check/instances/${name}.mps")
        # add a test and test properties
        add_test(NAME ${name}${setting} COMMAND $<TARGET_FILE:highs> ${setting} -f ${inst})
        set_tests_properties (${name}${setting} PROPERTIES PASS_REGULAR_EXPRESSION "(Solve plain|multi|Presolve|Postsolve): ${solutionstatus}")
    endforeach(setting)
endforeach(instance)
endmacro(add_instancetests)

# add tests for success and fail instances
add_instancetests(successInstances "OPTIMAL")
add_instancetests(failInstances "NOT-OPT")
add_instancetests(infeasibleInstances "NOT-OPT")
add_instancetests(unboundedInstances "NOT-OPT")
