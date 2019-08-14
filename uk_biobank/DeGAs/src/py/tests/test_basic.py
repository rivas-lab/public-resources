# -*- coding: utf-8 -*-

import unittest
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import rivas_decomposition_py

class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test(self):
        assert True

if __name__ == '__main__':
    unittest.main()

