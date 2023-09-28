# Callbacks

The HiGHS callback allows user actions to be performed within HiGHS. There is one generic callback method that can be defined by a user, with specific callback scenarios communicated to the user via a parameter. It can be given any name and, below, is called `userCallback`. Its definition is

```bash
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

The user's callback method is communicated to HiGHS via the method that in the HiGHS C++ class is
```bash
HighsStatus setCallback(void (*userCallback)(const int, const char*, const HighsCallbackDataOut*,
                        HighsCallbackDataIn*, void*), void* user_callback_data);
```
and, in the HiGHS C API is
```bash
HighsInt Highs_setCallback(
    void* highs,
    void (*userCallback)(const int, const char*,
                         const struct HighsCallbackDataOut*,
                         struct HighsCallbackDataIn*, void*),
    void* user_callback_data)
```
There current callback scenarios are set out below, and the particular callback is activated in C++ by calling

```bash
HighsStatus startCallback(const int callback_type);
```
and, in C, by calling
```bash
HighsInt Highs_startCallback(void* highs, const int callback_type);
```
, and de-activated in C++ by calling
```bash
HighsStatus stopCallback(const int callback_type);
```
and, in C, by calling
```bash
HighsInt Highs_stopCallback(void* highs, const int callback_type);
```

### Logging callback

The logging callback type is `kCallbackLogging` in the C++ enum
`HighsCallbackType` and, in C, `kHighsCallbackLogging`. The logging
type is passed as the member `log_type` in the `HighsCallbackDataOut`
struct, and the message is passed as the `const char*` parameter.

### Simplex interrupt callback

The simplex interrupt is called once every simplex iteration, and its
callback type is `kCallbackSimplexInterrupt` in the C++ enum
`HighsCallbackType`, and `kHighsCallbackSimplexInterrupt` in C. The
simplex iteration count is passed as the `simplex_iteration_count`
member of the `HighsCallbackDataOut` struct.

If the `user_interrupt` member of the `HighsCallbackDataIn` struct is
set to a nonzero value, then the simplex solver will be interrupted,
and HiGHS will return to the user.

### IPM interrupt callback

The IPM interrupt is called once every interior point
iteration, and its callback type is `kCallbackIpmInterrupt` in the C++
enum `HighsCallbackType`, and `kHighsCallbackIpmInterrupt` in C. The
IPM iteration count is passed as the `ipm_iteration_count` member of
the `HighsCallbackDataOut` struct.

If the `user_interrupt` member of the `HighsCallbackDataIn` struct is
set to a nonzero value, then the IPM solver will be interrupted, and
HiGHS will return to the user.

### MIP improving solution callback

The MIP improving solution is called whenever the MIP
solver identifies an improving integer feasible solution, and its
callback type is `kCallbackMipImprovingSolution` in the C++ enum
`HighsCallbackType`, and `kHighsCallbackMipImprovingSolution` in C. A
pointer to the improving solution is passed as the
`objective_function_value` and `mip_solution` members of the
`HighsCallbackDataOut` struct.


### MIP logging callback

The MIP logging callback is called once every time MIP logging takes
place, and its callback type is `kCallbackMipLogging` in the C++ enum
`HighsCallbackType`, and `kHighsCallbackMipLogging` in C.

### MIP interrupt callback

The MIP interrupt is called once every interior point
iteration, and its callback type is `kCallbackMipInterrupt` in the C++
enum `HighsCallbackType`, and `kHighsCallbackMipInterrupt` in C. The
simplex iteration count is passed as the `simplex_iteration_count`
member of the `HighsCallbackDataOut` struct.

If the `user_interrupt` member of the `HighsCallbackDataIn` struct is
set to a nonzero value, then the MIP solver will be interrupted,
and HiGHS will return to the user.


### MIP callback data

For each of the MIP callbacks, the `HighsCallbackDataOut` struct will have the following values set

* `running_time`: execution time of HiGHS
* `objective_function_value`: the objective function value of the best integer feasible solution found
* `num_nodes`: the number of MIP nodes explored to date
* `primal_bound`: the primal bound
* `dual_bound`: the dual bound
* `mip_gap`: the (relative) difference between tht primal and dual bounds

### User interrupt

For the non-logging callbacks, if the `user_interrupt` member of the
`HighsCallbackDataIn` struct is set to a nonzero value, then the
corresponding solver will be interrupted, and HiGHS will return to the
user.


