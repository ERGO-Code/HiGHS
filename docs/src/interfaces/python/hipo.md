<a id="hipo-in-python"></a>

# HiPO in Python

The HiPO Interior Point Method (IPM) solver currently uses external
dependencies to provide enhanced performance for linear and quadratic
programming problems. The required dependencies are packaged in
`highspy-extras`. The packaged dependencies have licensing terms different
from HiGHS, such as Apache 2.0. Other algorithms may also rely on
`highspy-extras` in the future.

HiPO can enhance performance on many large problem instances. It is not very
well suited for smaller or easier LPs.

<a id="installation"></a>

## Installation

Install directly:

```shell
pip install highspy-extras
```

Or install via the `highspy` optional dependency:

```shell
pip install highspy[extras]
```

At present, the optional dependency installs support needed for HiPO. Both
`highspy` and `highspy-extras` are available on PyPI and conda-forge.

<a id="usage"></a>

## Usage

When `highspy-extras` is installed, HiGHS can use algorithms that depend on
these external libraries. At present this primarily means the HiPO solver.
Note that `highspy-extras` is automatically consumed by `highspy` and does
not need to be imported manually. You can explicitly select HiPO:

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

For debugging library packaging issues, you can also query the ABI version
reported directly by the shared library:

```python
import highspy_extras

print(highspy_extras.__version__)
print(highspy_extras.get_library_version())
```

<a id="local-installation-requirements"></a>

## Local Installation Requirements

To install locally, you will need:

* Python >= 3.9
* A BLAS library (bundled or system)

To install locally without an existing OpenBLAS installation, run:

```shell
python -m pip install ./highspy-extras
python -m pip install .
```

To install locally with an existing OpenBLAS installation, run:

```shell
python -m pip install ./highspy-extras --config-settings=cmake.define.BUILD_OPENBLAS=OFF
python -m pip install .
```

If the OpenBLAS installation path is not in the default search path, it can
be provided with:

```shell
--config-settings=cmake.define.BLAS_LIBRARIES=/path/to/openblas/library
```

<a id="uninstall"></a>

## Uninstall

To remove HiPO support and go back to the MIT-licensed `highspy`, remove
`highspy-extras`:

```shell
pip uninstall highspy-extras
```

<a id="license"></a>

## License

Apache 2.0 - see the license and `THIRD_PARTY_NOTICES` in the
[HiGHS repository](https://github.com/ERGO-Code/HiGHS) for details.
