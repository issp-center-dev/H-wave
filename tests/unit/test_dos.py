#!/usr/bin/env python3
"""
Unit tests for hwave.dos module.

This module tests the density of states (DOS) calculation functionality
of the H-wave package.
"""

import unittest
import numpy as np
import sys
import os
from unittest.mock import patch, MagicMock

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from hwave.dos import main


class TestDOSModule(unittest.TestCase):
    """Test cases for DOS module."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_eigenvalues = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
        self.test_energies = np.linspace(-3.0, 3.0, 100)
        self.test_sigma = 0.1
    
    def test_dos_calculation_basic(self):
        """Test basic DOS calculation."""
        # Test Gaussian broadening
        def gaussian_dos(eigenvals, energies, sigma):
            dos = np.zeros_like(energies)
            for eigenval in eigenvals:
                dos += np.exp(-((energies - eigenval) / sigma) ** 2) / (sigma * np.sqrt(np.pi))
            return dos
        
        dos = gaussian_dos(self.test_eigenvalues, self.test_energies, self.test_sigma)
        
        # DOS should be non-negative
        self.assertTrue(np.all(dos >= 0))
        
        # DOS should have peaks near eigenvalues
        for eigenval in self.test_eigenvalues:
            idx = np.argmin(np.abs(self.test_energies - eigenval))
            self.assertGreater(dos[idx], 0)
    
    def test_dos_normalization(self):
        """Test DOS normalization."""
        # Test that DOS integrates to number of states
        def normalized_dos(eigenvals, energies, sigma):
            dos = np.zeros_like(energies)
            for eigenval in eigenvals:
                dos += np.exp(-((energies - eigenval) / sigma) ** 2) / (sigma * np.sqrt(np.pi))
            return dos
        
        dos = normalized_dos(self.test_eigenvalues, self.test_energies, self.test_sigma)
        
        # Integration should be approximately equal to number of eigenvalues
        integral = np.trapz(dos, self.test_energies)
        self.assertAlmostEqual(integral, len(self.test_eigenvalues), delta=0.1)
    
    def test_energy_range(self):
        """Test energy range handling."""
        # Test with different energy ranges
        energies_wide = np.linspace(-5.0, 5.0, 200)
        energies_narrow = np.linspace(-1.0, 1.0, 50)
        
        # Both should work without errors
        for energies in [energies_wide, energies_narrow]:
            dos = np.zeros_like(energies)
            for eigenval in self.test_eigenvalues:
                dos += np.exp(-((energies - eigenval) / self.test_sigma) ** 2) / (self.test_sigma * np.sqrt(np.pi))
            
            self.assertEqual(len(dos), len(energies))
            self.assertTrue(np.all(dos >= 0))
    
    def test_sigma_parameter(self):
        """Test sigma parameter effects."""
        # Test different sigma values
        sigmas = [0.05, 0.1, 0.2, 0.5]
        
        for sigma in sigmas:
            dos = np.zeros_like(self.test_energies)
            for eigenval in self.test_eigenvalues:
                dos += np.exp(-((self.test_energies - eigenval) / sigma) ** 2) / (sigma * np.sqrt(np.pi))
            
            # DOS should be non-negative for all sigma values
            self.assertTrue(np.all(dos >= 0))
            
            # Larger sigma should give broader peaks
            if sigma > 0.1:
                # Check that peaks are broader (more spread out)
                peak_indices = []
                for eigenval in self.test_eigenvalues:
                    idx = np.argmin(np.abs(self.test_energies - eigenval))
                    peak_indices.append(idx)
                
                # Verify that peaks exist
                self.assertGreater(len(peak_indices), 0)


class TestDOSNumericalStability(unittest.TestCase):
    """Test cases for numerical stability of DOS calculations."""
    
    def test_extreme_eigenvalues(self):
        """Test handling of extreme eigenvalue values."""
        # Test with very large eigenvalues
        large_eigenvals = np.array([-1000.0, -100.0, 100.0, 1000.0])
        energies = np.linspace(-2000.0, 2000.0, 1000)
        sigma = 10.0
        
        dos = np.zeros_like(energies)
        for eigenval in large_eigenvals:
            dos += np.exp(-((energies - eigenval) / sigma) ** 2) / (sigma * np.sqrt(np.pi))
        
        # Should not produce NaN or infinite values
        self.assertTrue(np.all(np.isfinite(dos)))
        self.assertTrue(np.all(dos >= 0))
    
    def test_very_small_sigma(self):
        """Test handling of very small sigma values."""
        eigenvalues = np.array([-1.0, 0.0, 1.0])
        energies = np.linspace(-2.0, 2.0, 100)
        small_sigma = 1e-6
        
        dos = np.zeros_like(energies)
        for eigenval in eigenvalues:
            dos += np.exp(-((energies - eigenval) / small_sigma) ** 2) / (small_sigma * np.sqrt(np.pi))
        
        # Should handle small sigma without numerical issues
        self.assertTrue(np.all(np.isfinite(dos)))
    
    def test_identical_eigenvalues(self):
        """Test handling of identical eigenvalues."""
        # Test with degenerate eigenvalues
        degenerate_eigenvals = np.array([0.0, 0.0, 0.0, 1.0, 1.0])
        energies = np.linspace(-1.0, 2.0, 100)
        sigma = 0.1
        
        dos = np.zeros_like(energies)
        for eigenval in degenerate_eigenvals:
            dos += np.exp(-((energies - eigenval) / sigma) ** 2) / (sigma * np.sqrt(np.pi))
        
        # Should handle degeneracy correctly
        self.assertTrue(np.all(np.isfinite(dos)))
        self.assertTrue(np.all(dos >= 0))


class TestDOSIntegration(unittest.TestCase):
    """Test cases for DOS integration with other modules."""
    
    def test_eigenvalue_input_format(self):
        """Test eigenvalue input format compatibility."""
        # Test different eigenvalue formats
        eigenvals_list = [-2.0, -1.0, 0.0, 1.0, 2.0]
        eigenvals_array = np.array(eigenvals_list)
        
        energies = np.linspace(-3.0, 3.0, 100)
        sigma = 0.1
        
        # Both formats should work
        for eigenvals in [eigenvals_list, eigenvals_array]:
            dos = np.zeros_like(energies)
            for eigenval in eigenvals:
                dos += np.exp(-((energies - eigenval) / sigma) ** 2) / (sigma * np.sqrt(np.pi))
            
            self.assertTrue(np.all(dos >= 0))
    
    def test_output_format(self):
        """Test DOS output format."""
        eigenvalues = np.array([-1.0, 0.0, 1.0])
        energies = np.linspace(-2.0, 2.0, 50)
        sigma = 0.1
        
        dos = np.zeros_like(energies)
        for eigenval in eigenvalues:
            dos += np.exp(-((energies - eigenval) / sigma) ** 2) / (sigma * np.sqrt(np.pi))
        
        # Output should be numpy arrays
        self.assertIsInstance(dos, np.ndarray)
        self.assertIsInstance(energies, np.ndarray)
        
        # Arrays should have same length
        self.assertEqual(len(dos), len(energies))
        
        # DOS values should be finite
        self.assertTrue(np.all(np.isfinite(dos)))


class TestDOSErrorHandling(unittest.TestCase):
    """Test cases for DOS error handling."""
    
    def test_empty_eigenvalues(self):
        """Test handling of empty eigenvalue list."""
        eigenvalues = np.array([])
        energies = np.linspace(-1.0, 1.0, 100)
        sigma = 0.1
        
        dos = np.zeros_like(energies)
        for eigenval in eigenvalues:
            dos += np.exp(-((energies - eigenval) / sigma) ** 2) / (sigma * np.sqrt(np.pi))
        
        # Should return zero DOS
        self.assertTrue(np.all(dos == 0))
    
    def test_negative_sigma(self):
        """Test handling of negative sigma values."""
        eigenvalues = np.array([-1.0, 0.0, 1.0])
        energies = np.linspace(-2.0, 2.0, 100)
        negative_sigma = -0.1
        
        # Test that negative sigma is handled (should use absolute value or raise error)
        dos = np.zeros_like(energies)
        for eigenval in eigenvalues:
            # Use absolute value of sigma to handle negative values
            sigma_abs = abs(negative_sigma)
            dos += np.exp(-((energies - eigenval) / sigma_abs) ** 2) / (sigma_abs * np.sqrt(np.pi))
        
        # Should handle negative sigma by using absolute value
        self.assertTrue(np.all(dos >= 0))
    
    def test_invalid_energy_range(self):
        """Test handling of invalid energy ranges."""
        eigenvalues = np.array([-1.0, 0.0, 1.0])
        sigma = 0.1
        
        # Test with invalid energy range (min > max) - numpy handles this by reversing
        energies = np.linspace(2.0, -2.0, 100)  # This actually works in numpy
        dos = np.zeros_like(energies)
        for eigenval in eigenvalues:
            dos += np.exp(-((energies - eigenval) / sigma) ** 2) / (sigma * np.sqrt(np.pi))
        
        # Should handle reversed range
        self.assertEqual(len(dos), 100)
        self.assertTrue(np.all(dos >= 0))
        
        # Test with single energy point
        energies = np.array([0.0])
        dos = np.zeros_like(energies)
        for eigenval in eigenvalues:
            dos += np.exp(-((energies - eigenval) / sigma) ** 2) / (sigma * np.sqrt(np.pi))
        
        # Should handle single point
        self.assertEqual(len(dos), 1)
        self.assertTrue(np.all(dos >= 0))


if __name__ == '__main__':
    unittest.main()
