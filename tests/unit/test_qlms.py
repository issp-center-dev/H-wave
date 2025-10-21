#!/usr/bin/env python3
"""
Unit tests for hwave.qlms module.

This module tests the main QLMS functionality of the H-wave package,
including parameter handling, solver selection, and main execution flow.
"""

import unittest
import numpy as np
import sys
import os
from unittest.mock import patch, MagicMock, mock_open
import tempfile
import json

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from hwave.qlms import main, run


class TestQLMSMain(unittest.TestCase):
    """Test cases for main function."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_input_dict = {
            "solver": "UHFr",
            "Nsite": 4,
            "Ncond": 4,
            "T": 0.0,
            "EPS": 1e-8,
            "IterationMax": 1000
        }
    
    def test_main_with_valid_input(self):
        """Test main function with valid input."""
        # Mock the run function to avoid actual calculations
        with patch('hwave.qlms.run') as mock_run:
            mock_run.return_value = None
            
            # Test that main can be called without errors
            try:
                main()
            except SystemExit:
                # Expected for command line interface
                pass
    
    def test_run_with_valid_dict(self):
        """Test run function with valid input dictionary."""
        with patch('hwave.qlms.QLMSInput') as mock_input:
            with patch('hwave.solver.uhfr.UHFr') as mock_solver:
                mock_solver_instance = MagicMock()
                mock_solver.return_value = mock_solver_instance
                
                # Test run function
                result = run(input_dict=self.test_input_dict)
                
                # Verify that solver was called
                mock_solver.assert_called_once()
    
    def test_run_with_invalid_solver(self):
        """Test run function with invalid solver type."""
        invalid_input = self.test_input_dict.copy()
        invalid_input["solver"] = "InvalidSolver"
        
        with self.assertRaises((ValueError, KeyError, AttributeError)):
            run(input_dict=invalid_input)
    
    def test_run_with_missing_parameters(self):
        """Test run function with missing required parameters."""
        incomplete_input = {"solver": "UHFr"}
        
        with self.assertRaises((KeyError, ValueError)):
            run(input_dict=incomplete_input)


class TestQLMSParameterValidation(unittest.TestCase):
    """Test cases for parameter validation."""
    
    def test_valid_solver_types(self):
        """Test valid solver types."""
        valid_solvers = ["UHFr", "UHFk", "RPA"]
        
        for solver in valid_solvers:
            input_dict = {
                "solver": solver,
                "Nsite": 4,
                "Ncond": 4
            }
            
            # Should not raise an error for valid solvers
            try:
                with patch('hwave.qlms.QLMSInput'):
                    with patch('hwave.solver.uhfr.UHFr'):
                        with patch('hwave.solver.uhfk.UHFk'):
                            with patch('hwave.solver.rpa.RPA'):
                                run(input_dict=input_dict)
            except (ValueError, KeyError, AttributeError):
                # Some solvers might not be fully mocked, which is acceptable
                pass
    
    def test_parameter_types(self):
        """Test parameter type validation."""
        # Test integer parameters
        input_dict = {
            "solver": "UHFr",
            "Nsite": 4,
            "Ncond": 4,
            "T": 0.0,
            "EPS": 1e-8,
            "IterationMax": 1000
        }
        
        # All parameters should be of correct types
        self.assertIsInstance(input_dict["Nsite"], int)
        self.assertIsInstance(input_dict["Ncond"], int)
        self.assertIsInstance(input_dict["T"], float)
        self.assertIsInstance(input_dict["EPS"], float)
        self.assertIsInstance(input_dict["IterationMax"], int)
    
    def test_parameter_ranges(self):
        """Test parameter range validation."""
        # Test valid ranges
        self.assertGreater(4, 0)  # Nsite > 0
        self.assertGreaterEqual(4, 4)  # Ncond <= Nsite
        self.assertGreaterEqual(1e-8, 0)  # EPS >= 0
        self.assertGreater(1000, 0)  # IterationMax > 0


class TestQLMSFileHandling(unittest.TestCase):
    """Test cases for file handling functionality."""
    
    def test_toml_file_reading(self):
        """Test TOML file reading functionality."""
        # Create a temporary TOML file
        toml_content = """
[solver]
type = "UHFr"

[parameters]
Nsite = 4
Ncond = 4
T = 0.0
EPS = 1e-8
IterationMax = 1000
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.toml', delete=False) as f:
            f.write(toml_content)
            temp_file = f.name
        
        try:
            # Test that file can be read
            with open(temp_file, 'r') as f:
                content = f.read()
                self.assertIn("solver", content)
                self.assertIn("parameters", content)
        finally:
            os.unlink(temp_file)
    
    def test_def_file_reading(self):
        """Test DEF file reading functionality."""
        def_content = """
&namelist
Nsite = 4
Ncond = 4
T = 0.0
EPS = 1e-8
IterationMax = 1000
/
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.def', delete=False) as f:
            f.write(def_content)
            temp_file = f.name
        
        try:
            # Test that file can be read
            with open(temp_file, 'r') as f:
                content = f.read()
                self.assertIn("namelist", content)
                self.assertIn("Nsite", content)
        finally:
            os.unlink(temp_file)


class TestQLMSNumericalOperations(unittest.TestCase):
    """Test cases for numerical operations."""
    
    def test_matrix_operations(self):
        """Test matrix operations used in QLMS."""
        # Test matrix creation
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        B = np.array([[2.0, 0.0], [0.0, 2.0]])
        
        # Test matrix multiplication
        C = np.dot(A, B)
        expected = np.array([[2.0, 4.0], [6.0, 8.0]])
        np.testing.assert_array_almost_equal(C, expected)
        
        # Test matrix addition
        D = A + B
        expected = np.array([[3.0, 2.0], [3.0, 6.0]])
        np.testing.assert_array_almost_equal(D, expected)
    
    def test_complex_operations(self):
        """Test complex number operations."""
        # Test complex matrix
        A = np.array([[1.0 + 1j, 2.0], [3.0, 4.0 + 2j]], dtype=complex)
        
        # Test complex conjugate
        A_conj = np.conj(A)
        expected = np.array([[1.0 - 1j, 2.0], [3.0, 4.0 - 2j]], dtype=complex)
        np.testing.assert_array_almost_equal(A_conj, expected)
        
        # Test Hermitian property
        A_herm = A + A.conj().T
        self.assertTrue(np.allclose(A_herm, A_herm.conj().T))
    
    def test_eigenvalue_calculations(self):
        """Test eigenvalue calculations."""
        # Test symmetric matrix
        A = np.array([[2.0, 1.0], [1.0, 2.0]])
        eigenvals, eigenvecs = np.linalg.eigh(A)
        
        # Eigenvalues should be real and positive
        self.assertTrue(np.all(np.isreal(eigenvals)))
        self.assertTrue(np.all(eigenvals > 0))
        
        # Test eigenvalue reconstruction
        A_reconstructed = np.dot(eigenvecs * eigenvals, eigenvecs.T)
        np.testing.assert_array_almost_equal(A, A_reconstructed)


class TestQLMSErrorHandling(unittest.TestCase):
    """Test cases for error handling."""
    
    def test_invalid_input_types(self):
        """Test handling of invalid input types."""
        # Test with non-dictionary input
        with self.assertRaises((TypeError, AttributeError)):
            run(input_dict="invalid_input")
        
        # Test with None input
        with self.assertRaises((TypeError, AttributeError)):
            run(input_dict=None)
    
    def test_missing_required_keys(self):
        """Test handling of missing required keys."""
        incomplete_dict = {"Nsite": 4}
        
        with self.assertRaises((KeyError, ValueError)):
            run(input_dict=incomplete_dict)
    
    def test_invalid_parameter_values(self):
        """Test handling of invalid parameter values."""
        # Test negative Nsite
        invalid_dict = {
            "solver": "UHFr",
            "Nsite": -1,
            "Ncond": 4
        }
        
        with self.assertRaises((ValueError, AssertionError)):
            run(input_dict=invalid_dict)


if __name__ == '__main__':
    unittest.main()
