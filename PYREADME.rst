pyHiGHS
=======

Cython wrappers over `HiGHS <https://github.com/ERGO-Code/HiGHS>`_.  This project is a fork of the main repo: `fork <https://github.com/mckib2/HiGHS/>`_.

The package is available on pypi as `scikit-highs <https://pypi.org/project/scikit-highs/>`_ and successfully builds on Ubuntu 18, gcc 7.5.0, Python 3.6.9.

It has dependencies on numpy and Cython that require these to be installed before installation of scikit-highs.  Once numpy and Cython are installed, install using pip like this:

.. code:: bash

    pip install scikit-highs

Example usage from script:

.. code:: python

    from pyHiGHS import highs_wrapper
    from pyHiGHS import linprog_mps

Examples can be run like this once it has been pip-installed:

.. code:: bash

    python -m pyHiGHS.examples.linprog_interface
    python -m pyHiGHS.examples.solve_mps

It is anticipated that scikit-highs can only be installed on Linux systems (including WSL) as the HiGHS project is currently not building under Windows (see `issue #270 <https://github.com/ERGO-Code/HiGHS/issues/270>`_).  It will take a second to compile when installing, I did not create any wheels, etc, so always compiles from source.
