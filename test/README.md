# Test suite for UHFr and UHFk

This test suite runs UHFr and UHFk solvers with various interaction types
and checks the results with the reference data.
It is prepared to work with ctest framework. 

## How to run

```bash
mkdir build && cd build
cmake ..
make test
```