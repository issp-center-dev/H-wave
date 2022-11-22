# Test suite for UHFr and UHFk

This test suite runs UHFr and UHFk solvers with various interaction types
and checks the results with the reference data.
It is prepared to work with python unittest and ctest framework. 

## How to run

### Using python unittest

```bash
python3 -m unittest
```

or run tests using `pytest` at the top directory.
(Usually you need to install `pytest` package.)

```bash
pytest -v
```

These command finds a series of tests from `test_uhf.py`
and other scripts in `tests` directory.
You may add test scripts in the `unittest` manner to perform various checks. 

### Using ctest

To run the tests prepared in the cmake/ctest framework, type in as follows:

```bash
mkdir build && cd build
cmake ..
make test
```

The test suite uses CMakeList.txt, tests/runtest.sh, and tests/compfile.py.
