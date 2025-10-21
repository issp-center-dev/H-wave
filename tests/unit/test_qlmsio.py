#!/usr/bin/env python3
"""
Unit tests for hwave.qlmsio module.

This module tests the input/output functionality of the H-wave package,
including file reading, parameter validation, and data structure handling.
"""

import unittest
import tempfile
import os
import numpy as np
from unittest.mock import patch, mock_open
import sys
import io

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from hwave.qlmsio.read_input import QLMSInput
from hwave.qlmsio.read_input_k import QLMSkInput
from hwave.qlmsio.wan90 import Wannier90Input


class TestQLMSInput(unittest.TestCase):
    """Test cases for QLMSInput class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.valid_namelist = [
            "trans", "coulombinter", "coulombintra", "pairhop", 
            "hund", "exchange", "ising", "pairlift", "interall", 
            "initial", "onebodyg"
        ]
    
    def test_valid_namelist(self):
        """Test that valid_namelist contains expected items."""
        expected_items = [
            "trans", "coulombinter", "coulombintra", "pairhop", 
            "hund", "exchange", "ising", "pairlift", "interall", 
            "initial", "onebodyg"
        ]
        for item in expected_items:
            self.assertIn(item, QLMSInput.valid_namelist)
    
    def test_initialization_with_empty_file_list(self):
        """Test initialization with empty file list."""
        with self.assertRaises((FileNotFoundError, IndexError)):
            QLMSInput([])
    
    def test_initialization_with_nonexistent_files(self):
        """Test initialization with non-existent files."""
        with self.assertRaises(FileNotFoundError):
            QLMSInput(["nonexistent.def"])
    
    @patch('builtins.open', new_callable=mock_open)
    @patch('os.path.exists', return_value=True)
    def test_initialization_with_mock_files(self, mock_exists, mock_file):
        """Test initialization with mocked files."""
        # Mock file content
        mock_file.return_value.read.return_value = "1 2 3\n4 5 6\n"
        
        file_list = ["test.def"]
        qlms_input = QLMSInput(file_list)
        
        self.assertEqual(qlms_input.file_names, file_list)
        self.assertIsInstance(qlms_input.ham_param, dict)
        self.assertIsInstance(qlms_input.green, dict)
    
    def test_get_param_method(self):
        """Test get_param method."""
        with patch('builtins.open', new_callable=mock_open):
            with patch('os.path.exists', return_value=True):
                qlms_input = QLMSInput(["test.def"])
                
                # Test getting existing parameter
                qlms_input.ham_param["TestParam"] = "test_value"
                self.assertEqual(qlms_input.get_param("TestParam"), "test_value")
                
                # Test getting non-existing parameter
                self.assertIsNone(qlms_input.get_param("NonExistentParam"))


class TestQLMSkInput(unittest.TestCase):
    """Test cases for QLMSkInput class."""
    
    def setUp(self):
        """Set up test fixtures."""
        pass
    
    @patch('builtins.open', new_callable=mock_open)
    @patch('os.path.exists', return_value=True)
    def test_initialization_with_mock_files(self, mock_exists, mock_file):
        """Test initialization with mocked files."""
        # Mock file content
        mock_file.return_value.read.return_value = "1 2 3\n4 5 6\n"
        
        file_list = ["test.def"]
        qlmsk_input = QLMSkInput(file_list)
        
        self.assertEqual(qlmsk_input.file_names, file_list)
        self.assertIsInstance(qlmsk_input.ham_param, dict)
        self.assertIsInstance(qlmsk_input.green, dict)


class TestWannier90Input(unittest.TestCase):
    """Test cases for Wannier90Input class."""
    
    def setUp(self):
        """Set up test fixtures."""
        pass
    
    def test_initialization_with_empty_file_list(self):
        """Test initialization with empty file list."""
        with self.assertRaises((FileNotFoundError, IndexError)):
            Wannier90Input([])
    
    @patch('builtins.open', new_callable=mock_open)
    @patch('os.path.exists', return_value=True)
    def test_initialization_with_mock_files(self, mock_exists, mock_file):
        """Test initialization with mocked files."""
        # Mock file content
        mock_file.return_value.read.return_value = "1 2 3\n4 5 6\n"
        
        file_list = ["test.def"]
        wan90_input = Wannier90Input(file_list)
        
        self.assertEqual(wan90_input.file_names, file_list)
        self.assertIsInstance(wan90_input.ham_param, dict)


class TestUtilityFunctions(unittest.TestCase):
    """Test cases for utility functions."""
    
    def test_numpy_import(self):
        """Test that numpy is properly imported."""
        import numpy as np
        self.assertIsNotNone(np)
        self.assertTrue(hasattr(np, 'array'))
        self.assertTrue(hasattr(np, 'zeros'))
        self.assertTrue(hasattr(np, 'ones'))
    
    def test_requests_import(self):
        """Test that requests is properly imported."""
        from requests.structures import CaseInsensitiveDict
        self.assertIsNotNone(CaseInsensitiveDict)
        
        # Test CaseInsensitiveDict functionality
        test_dict = CaseInsensitiveDict()
        test_dict["TestKey"] = "test_value"
        self.assertEqual(test_dict["testkey"], "test_value")
        self.assertEqual(test_dict["TESTKEY"], "test_value")


if __name__ == '__main__':
    unittest.main()
