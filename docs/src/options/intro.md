# [Introduction](@id options-intro)

The options that control HiGHS are of four types: `boolean`, `integer`, `double`
and `string`. Their values can be specified:

 * via the command line when running the [Executable](@ref executable)
 * via method calls when running HiGHS in an application.

## Options file

When running the [Executable](@ref executable) via the command line, some options values
can be set explicitly in the command, and all options can be set by means of an
options file.

A sample options file, giving documentation of all the options is written to the
console by the command:

```shell
$ bin/highs --options_file=""
```

## Option methods

The following code illustrates how to access the HiGHS options in Python. To set the value of option `name`, call:

``` python
status = h.setOptionValue(name, value)
```

where the value passed can be an identifier of the appropriate type, or an
explicit value.

To get the value of option `name`, call:

``` python
[status, value] = h.getOptionValue(name)
```

To get the type of option `name`, call:

``` python
[status, type] = h.getOptionType(name)
```

Examples of calls to options methods are given in the [examples section](@ref example-py).

