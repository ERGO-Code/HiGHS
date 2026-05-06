# Callbacks

The HiGHS callback allows user actions to be performed within HiGHS. There is
one generic callback method that can be defined by a user, with specific
callback scenarios communicated to the user via a parameter. Particular
callbacks must be activated (and can be deactivated) as described below. The
user callback can be given any name and, below, is called `userCallback`. Its
definition is

```cpp
void userCallback(const int callback_type,
                  const char* message,
                  const HighsCallbackDataOut* data_out,
                  HighsCallbackDataIn* data_in,
                  void* user_callback_data);
```

where

* `callback_type` is the type of callback
* `message` is a line of text associated with the callback.
* `data_out` communicates data from HiGHS to the user
* `data_in` communicates data from the user to HiGHS
* `user_callback_data` allows the user to pass data to the callback

The logging callback type is a cast of the relevant member of the C++ enum
`HighsCallbackType`, and is available in C as a constant.

The user's callback method is communicated to HiGHS via the method that in the HiGHS C++ class is
```cpp
HighsStatus setCallback(void (*userCallback)(const int, const char*, const HighsCallbackDataOut*,
                        HighsCallbackDataIn*, void*), void* user_callback_data);
```
and, in the HiGHS C API is
```cpp
HighsInt Highs_setCallback(
    void* highs,
    void (*userCallback)(const int, const char*,
                         const struct HighsCallbackDataOut*,
                         struct HighsCallbackDataIn*, void*),
    void* user_callback_data)
```
The current callback scenarios are set out below, and the particular callback is activated in C++ by calling

```cpp
HighsStatus startCallback(const int callback_type);
```
and, in C, by calling
```cpp
HighsInt Highs_startCallback(void* highs, const int callback_type);
```
, and de-activated in C++ by calling
```cpp
HighsStatus stopCallback(const int callback_type);
```
and, in C, by calling
```cpp
HighsInt Highs_stopCallback(void* highs, const int callback_type);
```

### User interrupt

For the non-logging callbacks below, if the `user_interrupt` member of the
`HighsCallbackDataIn` struct is set to a nonzero value, then the
corresponding solver will be interrupted, and HiGHS will return to the
user.

### Logging callback

The logging callback type is a cast of `kCallbackLogging` in the C++
enum `HighsCallbackType` and, in C, is the constant
`kHighsCallbackLogging`. The logging type is a cast of the particular
member of the enum class 'HighsLogType', and is available as a
constant in C. It is passed as the member `log_type` in the
`HighsCallbackDataOut` struct, and the message is passed as the `const
char*` parameter.

### Simplex interrupt callback

The simplex interrupt is called once every simplex iteration, and its
callback type is a cast of `kCallbackSimplexInterrupt` in the C++ enum
`HighsCallbackType`, and the `kHighsCallbackSimplexInterrupt` constant
in C. The simplex iteration count is passed as the
`simplex_iteration_count` member of the `HighsCallbackDataOut` struct.

### IPM interrupt callback

The IPM interrupt is called once every interior point iteration, and
its callback type is a cast of `kCallbackIpmInterrupt` in the C++ enum
`HighsCallbackType`, and the `kHighsCallbackIpmInterrupt` constant in
C. The IPM iteration count is passed as the `ipm_iteration_count`
member of the `HighsCallbackDataOut` struct.

### MIP improving solution callback

The MIP improving solution is called whenever the MIP solver
identifies an improving integer feasible solution, and its callback
type is a cast of `kCallbackMipImprovingSolution` in the C++ enum
`HighsCallbackType`, and the `kHighsCallbackMipImprovingSolution`
constant in C. A pointer to the improving solution is passed as the
`objective_function_value` and `mip_solution` members of the
`HighsCallbackDataOut` struct.

### MIP solution callback

The MIP solution is called whenever the MIP solver
identifies an integer feasible solution, and its callback
type is a cast of `kCallbackMipSolution` in the C++ enum
`HighsCallbackType`, and the `kHighsCallbackMipSolution`
constant in C. A pointer to the solution is passed as the
`objective_function_value` and `mip_solution` members of the
`HighsCallbackDataOut` struct.

### MIP logging callback

The MIP logging callback is called once every time MIP logging takes
place, and its callback type is a cast of `kCallbackMipLogging` in the
C++ enum `HighsCallbackType`, and `kHighsCallbackMipLogging` in C.

### MIP interrupt callback

The MIP interrupt callback is called when the MIP solver checks
whether computation limits (such as time, node, leaves, improving
solutions, and target objective) have been reached, and its callback
type is a cast of `kCallbackMipInterrupt` in the C++ enum
`HighsCallbackType`, and the `kHighsCallbackMipInterrupt` constant in
C. The simplex iteration count is passed as the
`simplex_iteration_count` member of the `HighsCallbackDataOut` struct.

### MIP user solution callback

The MIP user solution callback is called after setting up the MIP
solver, at five points when exploring the root node, and before each
dive in the branch-and-bound tree search. The origin of the call is
given by the value of the `external_solution_query_origin` member of
the `HighsCallbackDataOut` struct. The aim is to allow potential
feasible primal solutions generated externally to be passed to the
HiGHS MIP solver. Note that by introducing data to the HiGHS MIP
solver at the user's discretion, its behaviour will generally be
non-deterministic.

### MIP cut pool callback

The MIP cut pool callback is called after generating cuts at the root
node, and its callback type is a cast of `kCallbackMipGetCutPool` in
the C++ enum `HighsCallbackType`, and the
`kHighsCallbackMipGetCutPool` constant in C. The data supplied
consists of the compressed sparse column representation of the cutpool
constraint matrix, and the corresponding lower and upper bounds.

## Callback data output

The `HighsCallbackDataOut` struct supplies data to the user that is
relevant to the particular callback. The general data are

* `log_type`: An integer cast of the `HighsLogType` value, indicating the severity of the logging message--relevant to the logging callback.
* `running_time`: The excution time of HiGHS--relevant to the interrupt callbacks.
* `simplex_iteration_count`: The number of simplex iterations performed--relevant to the simplex interrupt callback.
* `ipm_iteration_count`: The number of IPM iterations performed--relevant to the IPM interrupt callback.
* `pdlp_iteration_count`: The number of PDLP iterations performed--relevant to the PDLP interrupt callback.

### MIP callback data

For each of the MIP callbacks, the following `HighsCallbackDataOut`
struct members will also have their value set

* `objective_function_value`: the objective function value of the best integer feasible solution found
* `mip_node_count`: the number of MIP nodes explored to date
* `mip_primal_bound`: the primal bound
* `mip_dual_bound`: the dual bound
* `mip_gap`: the (relative) difference between the primal and dual bounds
