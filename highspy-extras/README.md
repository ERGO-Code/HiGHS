# highspy-extras

Extension package for [highspy](https://pypi.org/project/highspy/) that enables access to external dependencies with licensing terms different from HiGHS, such as Apache 2.0.

The HiPO Interior Point Method (IPM) solver currently uses these external dependencies to provide enhanced performance for linear and quadratic programming problems. Other algorithms may also rely on `highspy-extras` in the future.

## Installation

Install directly:

```bash
pip install highspy-extras
```

Or install via the highspy optional dependency:

```bash
pip install highspy[extras]
```

At present, the optional dependency installs support needed for HiPO.

## Usage

When `highspy-extras` is installed, HiGHS can use algorithms that depend on these external libraries. At present this primarily means the HiPO solver. You can explicitly select HiPO:

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

For debugging library packaging issues, you can also query the ABI version reported directly by the shared library:

```python
import highspy-extras

print(highspy-extras.__version__)
print(highspy-extras.get_library_version())
```

## Requirements

- Python >= 3.8

- BLAS library
The highspy-extras wheels on PyPI and conda-forge ship with a bundled OpenBLAS.

When installing locally from source, the requirement is:
  - OpenBLAS on Windows and Linux:
    - Install with vcpkg
    - Install with apt on Linux
    - Set BUILD_OPENBLAS=ON in cmake, to fetch and build OpenBLAS as a subproject
  - Apple Accelerate on MacOS:
    - MacOS comes with the Accelerate framework pre-installed

## License

Apache 2.0 - see the license and `THIRD_PARTY_NOTICES` in the [HiGHS repository](https://github.com/ERGO-Code/HiGHS) for details.
