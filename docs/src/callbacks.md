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

### MIP callback data

For each of the MIP callbacks, the following `HighsCallbackDataOut` struct members will have value set

* `running_time`: execution time of HiGHS
* `objective_function_value`: the objective function value of the best integer feasible solution found
* `mip_node_count`: the number of MIP nodes explored to date
* `mip_primal_bound`: the primal bound
* `mip_dual_bound`: the dual bound
* `mip_gap`: the (relative) difference between the primal and dual bounds



