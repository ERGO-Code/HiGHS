# highspy-extras

Extension package for [highspy](https://pypi.org/project/highspy/) that enables access to external dependencies with licensing terms different from HiGHS, such as Apache 2.0.

The HiPO Interior Point Method (IPM) solver currently uses these external dependencies to provide enhanced performance for linear programming problems. Other algorithms may also rely on `highspy-extras` in the future.

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
import highspy_extras

print(highspy_extras.__version__)
print(highspy_extras.get_library_version())
```

## Requirements

- Python >= 3.8
- BLAS library (bundled or system)

## License

Apache - see the [HiGHS repository](https://github.com/ERGO-Code/HiGHS) for details.
