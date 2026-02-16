# highspy-hipo

HiPO IPM solver extension for [highspy](https://pypi.org/project/highspy/).

HiPO is an advanced Interior Point Method (IPM) solver that provides enhanced performance for linear programming problems.

## Installation

Install directly:

```bash
pip install highspy-hipo
```

Or install via the highspy optional dependency:

```bash
pip install highspy[hipo]
```

## Usage

When `highspy-hipo` is installed, HiGHS will automatically detect and use the HiPO solver when appropriate. You can explicitly select the HiPO solver:

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

## Requirements

- Python >= 3.8
- highspy >= 1.12.0
- BLAS library (bundled or system)

## License

Apache - see the [HiGHS repository](https://github.com/ERGO-Code/HiGHS) for details.
