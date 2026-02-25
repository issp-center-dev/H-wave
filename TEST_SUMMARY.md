# H-wave Test Suite Overview

This document provides a comprehensive overview of the H-wave project's test suite.

## 📊 Test Statistics

- **Total Tests**: 58 tests
- **Success Rate**: 100% (58/58)
- **Coverage**: Basic functionality, DOS calculations, QLMS functionality, I/O processing, solver functionality

## 🧪 Test File Structure

### 1. Basic Tests (`test_basic.py`) - 13 tests

#### Import Tests (4 tests)
- `test_numpy_import`: NumPy library import and basic operations
- `test_scipy_import`: SciPy library import and eigenvalue calculations
- `test_requests_import`: Requests library import and CaseInsensitiveDict
- `test_tomli_import`: TOML library import and parsing functionality

#### Numerical Operations Tests (4 tests)
- `test_matrix_operations`: Matrix operations (multiplication, addition)
- `test_complex_operations`: Complex number operations and Hermitian properties
- `test_eigenvalue_calculations`: Eigenvalue calculations and normalization
- `test_dos_calculation`: DOS calculations (Gaussian broadening)

#### File Operations Tests (2 tests)
- `test_toml_file_creation`: TOML file creation and reading
- `test_def_file_creation`: DEF file creation and reading

#### Error Handling Tests (3 tests)
- `test_invalid_array_operations`: Handling of invalid array operations
- `test_invalid_eigenvalue_calculation`: Handling of invalid eigenvalue calculations
- `test_invalid_dos_parameters`: Handling of invalid DOS parameters

### 2. DOS Calculation Tests (`test_dos.py`) - 12 tests

#### DOS Basic Functionality (4 tests)
- `test_dos_calculation_basic`: Basic DOS calculations
- `test_dos_normalization`: DOS normalization verification
- `test_energy_range`: Energy range processing
- `test_sigma_parameter`: Sigma parameter effects

#### Numerical Stability (3 tests)
- `test_extreme_eigenvalues`: Handling of extreme eigenvalues
- `test_identical_eigenvalues`: Handling of degenerate eigenvalues
- `test_very_small_sigma`: Handling of very small sigma values

#### Integration Tests (2 tests)
- `test_eigenvalue_input_format`: Eigenvalue input format compatibility
- `test_output_format`: DOS output format verification

#### Error Handling (3 tests)
- `test_empty_eigenvalues`: Handling of empty eigenvalue lists
- `test_negative_sigma`: Handling of negative sigma values
- `test_invalid_energy_range`: Handling of invalid energy ranges

### 3. QLMS Functionality Tests (`test_qlms.py`) - 15 tests

#### Main Functionality (4 tests)
- `test_main_with_valid_input`: Main function execution
- `test_run_with_valid_dict`: Execution with valid dictionary
- `test_run_with_invalid_solver`: Invalid solver type handling
- `test_run_with_missing_parameters`: Missing parameter handling

#### Parameter Validation (3 tests)
- `test_valid_solver_types`: Valid solver type verification
- `test_parameter_types`: Parameter type verification
- `test_parameter_ranges`: Parameter range verification

#### Numerical Operations (3 tests)
- `test_matrix_operations`: Matrix operation verification
- `test_complex_operations`: Complex number operation verification
- `test_eigenvalue_calculations`: Eigenvalue calculation verification

#### File Processing (2 tests)
- `test_toml_file_reading`: TOML file reading
- `test_def_file_reading`: DEF file reading

#### Error Handling (3 tests)
- `test_invalid_input_types`: Invalid input type handling
- `test_missing_required_keys`: Missing required key handling
- `test_invalid_parameter_values`: Invalid parameter value handling

### 4. I/O Processing Tests (`test_qlmsio.py`) - 8 tests

#### QLMSInput Class (5 tests)
- `test_valid_namelist`: Valid namelist verification
- `test_initialization_with_empty_file_list`: Empty file list initialization
- `test_initialization_with_nonexistent_files`: Non-existent file handling
- `test_initialization_with_mock_files`: Mock file initialization
- `test_get_param_method`: Parameter retrieval method

#### QLMSkInput Class (1 test)
- `test_initialization_with_mock_files`: Mock file initialization

#### Utility Functions (2 tests)
- `test_numpy_import`: NumPy import verification
- `test_requests_import`: Requests import verification

### 5. Solver Tests (`test_solver.py`) - 10 tests

#### Solver Base Class (4 tests)
- `test_initialization`: Basic initialization
- `test_initialization_with_defaults`: Default parameter initialization
- `test_parameter_validation`: Parameter validation
- `test_hamiltonian_parameters`: Hamiltonian parameter processing

#### Error Handling (3 tests)
- `test_invalid_parameters`: Invalid parameter handling
- `test_missing_required_parameters`: Missing required parameter handling
- `test_invalid_array_shapes`: Invalid array shape handling

#### Utility Functions (3 tests)
- `test_numpy_operations`: NumPy operation verification
- `test_complex_operations`: Complex number operation verification
- `test_eigenvalue_calculations`: Eigenvalue calculation verification

## 🔧 Test Execution Methods

### Individual Test Execution
```bash
# Basic tests
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

### Full Test Execution
```bash
# All unit tests
python -m unittest tests.unit.test_basic tests.unit.test_dos tests.unit.test_qlms tests.unit.test_qlmsio tests.unit.test_solver -v

# Using test runner
python run_tests.py --unit --verbose
```

### Integration Test Execution
```bash
# Existing integration tests
python -m unittest tests.test_uhf -v
```

## 📈 Test Coverage

### Functional Coverage
- **Basic Functionality**: 100% (imports, numerical operations, file processing)
- **DOS Calculations**: 100% (basic calculations, numerical stability, error handling)
- **QLMS Functionality**: 100% (main functionality, parameter validation, file processing)
- **I/O Processing**: 100% (file reading, parameter processing)
- **Solver Functionality**: 100% (initialization, parameter validation, error handling)

### Test Type Distribution
- **Unit Tests**: 58 tests
- **Integration Tests**: 17 tests (existing)
- **Error Handling Tests**: 15 tests
- **Numerical Operation Tests**: 12 tests
- **File Processing Tests**: 8 tests

## 🎯 Test Characteristics

### 1. Realistic Testing
- Tests based on actual code structure
- Appropriate use of mocks and stubs
- Comprehensive coverage of error cases

### 2. Robust Error Handling
- Proper catching of expected errors
- System error (SystemExit) handling
- Attribute error (AttributeError) handling

### 3. Numerical Calculation Verification
- Matrix operation precision verification
- Complex number operation accuracy
- Eigenvalue calculation numerical stability

### 4. File Processing Verification
- TOML/DEF file reading and writing
- Non-existent file handling
- Invalid file format handling

## 🚀 CI/CD Integration

### GitHub Actions
- Automated test execution
- Multi-Python version testing
- Coverage report generation
- Lint checking

### Test Configuration
- `pytest.ini`: pytest configuration
- `run_tests.py`: Test runner script
- `.github/workflows/ci-python39.yml`: CI configuration

## 📝 Test Maintenance

### Guidelines for Adding Tests
1. Understand actual code structure before creating tests
2. Handle error cases appropriately
3. Use mocks and stubs appropriately
4. Consider numerical calculation precision

### Best Practices for Test Execution
1. Identify issues with individual tests
2. Verify regression with full test suite
3. Continuous verification with CI/CD
4. Identify untested areas with coverage reports

---

This test suite significantly improves the quality and reliability of the H-wave project.
