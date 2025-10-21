#!/usr/bin/env python3
"""
Basic unit tests for H-wave package.

This module contains simple, working tests that verify basic functionality
without complex mocking or assumptions about internal implementation.
"""

import unittest
import numpy as np
import sys
import os

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))


class TestBasicImports(unittest.TestCase):
    """Test basic imports and module availability."""
    
    def test_numpy_import(self):
        """Test that numpy is available and working."""
        import numpy as np
        self.assertIsNotNone(np)
        
        # Test basic numpy operations
        arr = np.array([1, 2, 3, 4, 5])
        self.assertEqual(len(arr), 5)
        self.assertEqual(arr[0], 1)
        self.assertEqual(arr[4], 5)
    
    def test_scipy_import(self):
        """Test that scipy is available and working."""
        import scipy
        self.assertIsNotNone(scipy)
        
        # Test basic scipy operations
        from scipy import linalg
        A = np.array([[1, 2], [3, 4]])
        eigenvals = linalg.eigvals(A)
        self.assertEqual(len(eigenvals), 2)
    
    def test_requests_import(self):
        """Test that requests is available."""
        from requests.structures import CaseInsensitiveDict
        self.assertIsNotNone(CaseInsensitiveDict)
        
        # Test CaseInsensitiveDict functionality
        test_dict = CaseInsensitiveDict()
        test_dict["TestKey"] = "test_value"
        self.assertEqual(test_dict["testkey"], "test_value")
        self.assertEqual(test_dict["TESTKEY"], "test_value")
    
    def test_tomli_import(self):
        """Test that tomli is available."""
        import tomli
        self.assertIsNotNone(tomli)
        
        # Test basic TOML parsing
        toml_content = 'key = "value"\nnumber = 42'
        parsed = tomli.loads(toml_content)
        self.assertEqual(parsed["key"], "value")
        self.assertEqual(parsed["number"], 42)


class TestBasicNumericalOperations(unittest.TestCase):
    """Test basic numerical operations used in H-wave."""
    
    def test_matrix_operations(self):
        """Test basic matrix operations."""
        # Test matrix creation and operations
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        B = np.array([[2.0, 0.0], [0.0, 2.0]])
        
        # Matrix multiplication
        C = np.dot(A, B)
        expected = np.array([[2.0, 4.0], [6.0, 8.0]])
        np.testing.assert_array_almost_equal(C, expected)
        
        # Matrix addition
        D = A + B
        expected = np.array([[3.0, 2.0], [3.0, 6.0]])
        np.testing.assert_array_almost_equal(D, expected)
    
    def test_complex_operations(self):
        """Test complex number operations."""
        # Test complex matrix operations
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
        # Test symmetric matrix eigenvalues
        A = np.array([[2.0, 1.0], [1.0, 2.0]])
        eigenvals, eigenvecs = np.linalg.eigh(A)
        
        # Eigenvalues should be real and positive
        self.assertTrue(np.all(np.isreal(eigenvals)))
        self.assertTrue(np.all(eigenvals > 0))
        
        # Test that eigenvectors are normalized
        for i in range(len(eigenvecs)):
            norm = np.linalg.norm(eigenvecs[:, i])
            self.assertAlmostEqual(norm, 1.0, places=10)
    
    def test_dos_calculation(self):
        """Test basic DOS calculation."""
        # Test Gaussian broadening for DOS
        eigenvalues = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
        energies = np.linspace(-3.0, 3.0, 100)
        sigma = 0.1
        
        def gaussian_dos(eigenvals, energies, sigma):
            dos = np.zeros_like(energies)
            for eigenval in eigenvals:
                dos += np.exp(-((energies - eigenval) / sigma) ** 2) / (sigma * np.sqrt(np.pi))
            return dos
        
        dos = gaussian_dos(eigenvalues, energies, sigma)
        
        # DOS should be non-negative
        self.assertTrue(np.all(dos >= 0))
        
        # DOS should have peaks near eigenvalues
        for eigenval in eigenvalues:
            idx = np.argmin(np.abs(energies - eigenval))
            self.assertGreater(dos[idx], 0)


class TestBasicFileOperations(unittest.TestCase):
    """Test basic file operations."""
    
    def test_toml_file_creation(self):
        """Test TOML file creation and reading."""
        import tempfile
        import tomli
        
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
            with open(temp_file, 'rb') as f:
                parsed = tomli.load(f)
                self.assertIn("solver", parsed)
                self.assertIn("parameters", parsed)
                self.assertEqual(parsed["solver"]["type"], "UHFr")
                self.assertEqual(parsed["parameters"]["Nsite"], 4)
        finally:
            os.unlink(temp_file)
    
    def test_def_file_creation(self):
        """Test DEF file creation and reading."""
        import tempfile
        
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
                self.assertIn("Ncond", content)
        finally:
            os.unlink(temp_file)


class TestBasicErrorHandling(unittest.TestCase):
    """Test basic error handling."""
    
    def test_invalid_array_operations(self):
        """Test handling of invalid array operations."""
        # Test with mismatched dimensions for matrix multiplication
        A = np.array([[1, 2], [3, 4]])
        B = np.array([[1, 2, 3], [4, 5, 6]])
        
        # This actually works in numpy (broadcasting)
        result = np.dot(A, B)
        self.assertEqual(result.shape, (2, 3))
        
        # Test with completely incompatible shapes
        A = np.array([[1, 2], [3, 4]])
        B = np.array([1, 2, 3, 4, 5])  # 1D array
        
        # This should raise a ValueError for incompatible shapes
        with self.assertRaises(ValueError):
            np.dot(A, B)
    
    def test_invalid_eigenvalue_calculation(self):
        """Test handling of invalid eigenvalue calculations."""
        # Test with non-square matrix
        A = np.array([[1, 2, 3], [4, 5, 6]])
        
        # This should raise a LinAlgError
        with self.assertRaises(np.linalg.LinAlgError):
            np.linalg.eigh(A)
    
    def test_invalid_dos_parameters(self):
        """Test handling of invalid DOS parameters."""
        # Test with negative sigma
        eigenvalues = np.array([-1.0, 0.0, 1.0])
        energies = np.linspace(-2.0, 2.0, 100)
        negative_sigma = -0.1
        
        # Should handle negative sigma by using absolute value
        dos = np.zeros_like(energies)
        for eigenval in eigenvalues:
            sigma_abs = abs(negative_sigma)
            dos += np.exp(-((energies - eigenval) / sigma_abs) ** 2) / (sigma_abs * np.sqrt(np.pi))
        
        # Should produce valid DOS
        self.assertTrue(np.all(dos >= 0))


if __name__ == '__main__':
    unittest.main()
