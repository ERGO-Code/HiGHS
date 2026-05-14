# [HiPO in Python](@id hipo-in-python)

The HiPO Interior Point Method (IPM) solver currently uses external dependencies to provide enhanced performance for linear and quadratic programming problems. The required dependencies are packaged in the `highspy-extras`. The packaged dependencies have licensing terms different from HiGHS, such as Apache 2.0. Other algorithms may also rely on `highspy_extras` in the future.

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

When `highspy_extras` is installed, HiGHS can use algorithms that depend on these external libraries. At present this primarily means the HiPO solver. You can explicitly select HiPO:

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

## Local installation requirements

To install locally, you will need

- Python >= 3.8
- BLAS library (bundled or system)

## Uninstall

To remove the HiPO support and go back to the MIT-licenced highspy, simply remove highspy-extras:

```bash
pip uninstall highspy-extras
```

## License

Apache 2.0 - see the license and `THIRD_PARTY_NOTICES` the [HiGHS repository](https://github.com/ERGO-Code/HiGHS) for details.
