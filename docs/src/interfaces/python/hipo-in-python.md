# [HiPO in Python](@id hipo-in-python)

The HiPO Interior Point Method (IPM) solver currently uses external dependencies to provide enhanced performance for linear and quadratic programming problems. The required dependencies are packaged in the `highspy-extras`. The packaged dependencies have licensing terms different from HiGHS, such as Apache 2.0. Other algorithms may also rely on `highspy-extras` in the future.

HiPO can enhance performance on many large problem instances. It is not very well suited for smaller or easier LPs.

## Installation

Install directly:

```bash
pip install highspy-extras
```

Or install via the highspy optional dependency:

```bash
pip install highspy[extras]
```

At present, the optional dependency installs support needed for HiPO. Both `highspy` and `highspy-extras` are available on PyPI and conda-forge.

## Usage

When `highspy-extras` is installed, HiGHS can use algorithms that depend on these external libraries. At present this primarily means the HiPO solver. Note that `highspy-extras` is automatically consumed by `highspy` and does not need to be imported manually. You can explicitly select HiPO:

```python
import highspy

# Create a HiGHS instance
h = highspy.Highs()

# Load your model
h.readModel("model.mps")

# Set solver to use HiPO
h.setOptionValue("solver", "hipo")

# Solve
h.run()
```

### BLAS threads

HiPO runs OpenBLAS in serial for better performance. The number of OpenBLAS threads is a process-wide setting, so after HiPO has run, any other library in the same Python process that is linked to the same OpenBLAS library also runs it in serial. (NumPy and SciPy wheels from PyPI bundle their own copy of OpenBLAS, so are not affected, but this may not hold for other installations such as conda.) If this is a problem, disable it with

```python
h.setOptionValue("hipo_serial_blas", False)
```

and set the number of OpenBLAS threads yourself, e.g. using the `OPENBLAS_NUM_THREADS` environment variable. See [BLAS threads](@ref parallel-blas-threads) for details.

For debugging library packaging issues, you can also query the ABI version reported directly by the shared library:

```python
import highspy_extras

print(highspy_extras.__version__)
print(highspy_extras.get_library_version())
```

## Local installation requirements

To install locally, you will need

- Python >= 3.9
- BLAS library (bundled or system)

To install locally without an existing OpenBLAS installation, run

```
python -m pip install ./highspy-extras
python -m pip install .
```
To install locally with an existing OpenBLAS installation, run
```
python -m pip install ./highspy-extras --config-settings=cmake.define.BUILD_OPENBLAS=OFF
python -m pip install .
```
If the OpenBLAS installation path is not in the default set, it could be provided with
```
    --config-settings=cmake.define.BLAS_LIBRARIES=/path/to/openblas/library
```


## Uninstall

To remove the HiPO support and go back to the MIT-licenced highspy, simply remove highspy-extras:

```bash
pip uninstall highspy-extras
```

## License

Apache 2.0 - see the license and `THIRD_PARTY_NOTICES` the [HiGHS repository](https://github.com/ERGO-Code/HiGHS) for details.
