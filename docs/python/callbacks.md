# Callbacks

`highspy` exposes callback helpers for receiving solver events and, for MIP
models, interacting with incumbent solutions, logging, interrupts, and cut
pool events.

The high-level wrapper stores callback event objects on the
{py:attr}`highspy.Highs.callbacks` collection and provides convenience events
such as `cbMipGetCutPool`.

```python
import highspy
import numpy as np

h = highspy.Highs()
x = h.addBinaries(2, 3)
h.addConstrs(x[:, j].sum() == 1 for j in range(3))

def print_cuts(event):
    for cut in event.cuts:
        print(cut)

h.cbMipGetCutPool += print_cuts
h.minimize((np.ones((2, 3)) * x).sum())
```

## Callback API

```{eval-rst}
.. autosummary::

   highspy.HighsCallback
   highspy.HighsCallbackEvent
   highspy._core.cb.HighsCallbackInput
   highspy._core.cb.HighsCallbackOutput
   highspy._core.cb.HighsCallbackType
```
