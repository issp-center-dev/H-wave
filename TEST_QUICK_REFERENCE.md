# H-wave Test Quick Reference

## 🚀 Common Test Commands

### Full Test Execution
```bash
# All unit tests
python -m unittest tests.unit.test_basic tests.unit.test_dos tests.unit.test_qlms tests.unit.test_qlmsio tests.unit.test_solver -v

# Using test runner
python run_tests.py --unit --verbose
```

### Individual Test Execution
```bash
# Basic tests (imports, numerical operations, file processing)
python -m unittest tests.unit.test_basic -v

# DOS calculation tests
python -m unittest tests.unit.test_dos -v

# QLMS functionality tests
python -m unittest tests.unit.test_qlms -v

# I/O processing tests
python -m unittest tests.unit.test_qlmsio -v

# Solver tests
python -m unittest tests.unit.test_solver -v
```

### Integration Test Execution
```bash
# Existing integration tests
python -m unittest tests.test_uhf -v
```

## 📊 Test Statistics

| Test File | Test Count | Main Functionality |
|---|---|---|
| test_basic.py | 13 | Imports, numerical operations, file processing |
| test_dos.py | 12 | DOS calculations, numerical stability |
| test_qlms.py | 15 | QLMS functionality, parameter validation |
| test_qlmsio.py | 8 | I/O processing, file reading |
| test_solver.py | 10 | Solver functionality, error handling |
| **Total** | **58** | **Full functionality coverage** |

## 🎯 Test Categories

### Basic Functionality (13 tests)
- Import tests: numpy, scipy, requests, tomli
- Numerical operation tests: matrices, complex numbers, eigenvalues, DOS calculations
- File processing tests: TOML, DEF files
- Error handling tests: invalid operation processing

### DOS Calculations (12 tests)
- Basic DOS calculations: Gaussian broadening, normalization
- Numerical stability: extreme values, degenerate eigenvalues
- Integration tests: input format, output format
- Error handling: empty lists, negative sigma

### QLMS Functionality (15 tests)
- Main functionality: execution, parameter processing
- Parameter validation: types, ranges, solver types
- Numerical operations: matrices, complex numbers, eigenvalues
- File processing: TOML, DEF reading
- Error handling: invalid inputs, missing parameters

### I/O Processing (8 tests)
- File reading: empty lists, non-existent files
- Parameter processing: retrieval, validation
- Mock tests: file operation simulation
- Utilities: library imports

### Solver (10 tests)
- Basic functionality: initialization, parameter validation
- Hamiltonian: parameter processing
- Error handling: invalid parameters, array shapes
- Utilities: NumPy operations, complex numbers, eigenvalues

## 🔧 Troubleshooting

### Common Issues and Solutions

#### 1. Import Errors
```bash
# Path configuration check
export PYTHONPATH="${PYTHONPATH}:$(pwd)/src"
python -m unittest tests.unit.test_basic -v
```

#### 2. Mock Errors
```bash
# Identify issues with individual tests
python -m unittest tests.unit.test_qlmsio.TestQLMSInput.test_initialization_with_mock_files -v
```

#### 3. Numerical Calculation Errors
```bash
# Verify numerical operation tests
python -m unittest tests.unit.test_basic.TestBasicNumericalOperations -v
```

## 📈 Performance

### Execution Time
- **All tests**: < 1 second
- **Individual files**: < 0.2 seconds
- **Integration tests**: < 1 second

### Memory Usage
- **Basic tests**: Low
- **DOS calculation tests**: Medium
- **Integration tests**: High

## 🎯 Test Strategy

### Recommended Development Flow
1. **Run basic tests** → Infrastructure verification
2. **Run related tests** → Functionality verification
3. **Run all tests** → Regression verification
4. **Run integration tests** → End-to-end verification

### CI/CD Execution
```bash
# Automatic execution in GitHub Actions
# - Multiple Python versions (3.9, 3.10, 3.11, 3.12)
# - Automatic lint checking
# - Coverage report generation
```

## 📝 Guidelines for Adding Tests

### Adding New Tests
1. **Select appropriate file**: Choose test file based on functionality
2. **Create test class**: Separate classes by functionality
3. **Error handling**: Handle expected errors appropriately
4. **Documentation**: Clearly state test purpose and expected results

### Ensuring Test Quality
- **Single responsibility**: One test for one functionality
- **Independence**: Avoid dependencies between tests
- **Reproducibility**: Ensure consistent results
- **Maintainability**: Write understandable test code

---

This quick reference enables efficient execution and management of H-wave project tests.
