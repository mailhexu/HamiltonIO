#!/usr/bin/env python3
"""
Unit tests for the WignerData class.

This module contains comprehensive tests for the WignerData dataclass
used for handling Wigner-Seitz data from EPW calculations.
"""

import unittest
import tempfile
import os
import numpy as np
from wigner import WignerData


class TestWignerData(unittest.TestCase):
    """Unit tests for WignerData class."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_file = "test/up/wigner.fmt"
        self.assertTrue(os.path.exists(self.test_file), 
                      f"Test file {self.test_file} not found")

    def test_file_reading(self):
        """Test reading wigner data from file."""
        data = WignerData.from_file(self.test_file)
        
        # Check basic properties
        self.assertEqual(data.dims, 1)
        self.assertEqual(data.dims2, 1)
        self.assertEqual(data.nrr_k, 27)
        self.assertEqual(data.nrr_q, 27)
        self.assertEqual(data.nrr_g, 27)

    def test_array_shapes(self):
        """Test that arrays have correct Pythonic shapes."""
        data = WignerData.from_file(self.test_file)
        
        # Test electron data shapes
        self.assertEqual(data.irvec_k.shape, (27, 3))
        self.assertEqual(data.ndegen_k.shape, (27, 1, 1))
        self.assertEqual(data.wslen_k.shape, (27,))
        
        # Test phonon data shapes
        self.assertEqual(data.irvec_q.shape, (27, 3))
        self.assertEqual(data.ndegen_q.shape, (27, 1, 1))
        self.assertEqual(data.wslen_q.shape, (27,))
        
        # Test electron-phonon data shapes
        self.assertEqual(data.irvec_g.shape, (27, 3))
        self.assertEqual(data.ndegen_g.shape, (27, 1, 1))
        self.assertEqual(data.wslen_g.shape, (27,))

    def test_array_types(self):
        """Test that arrays have correct data types."""
        data = WignerData.from_file(self.test_file)
        
        # Integer arrays
        self.assertTrue(np.issubdtype(data.irvec_k.dtype, np.integer))
        self.assertTrue(np.issubdtype(data.ndegen_k.dtype, np.integer))
        self.assertTrue(np.issubdtype(data.irvec_q.dtype, np.integer))
        self.assertTrue(np.issubdtype(data.ndegen_q.dtype, np.integer))
        self.assertTrue(np.issubdtype(data.irvec_g.dtype, np.integer))
        self.assertTrue(np.issubdtype(data.ndegen_g.dtype, np.integer))
        
        # Float arrays
        self.assertTrue(np.issubdtype(data.wslen_k.dtype, np.floating))
        self.assertTrue(np.issubdtype(data.wslen_q.dtype, np.floating))
        self.assertTrue(np.issubdtype(data.wslen_g.dtype, np.floating))

    def test_sample_data_values(self):
        """Test specific data values from the file."""
        data = WignerData.from_file(self.test_file)
        
        # Test first k-vector
        np.testing.assert_array_equal(data.irvec_k[0], [-1, -1, -1])
        self.assertAlmostEqual(data.wslen_k[0], 1.732051, places=5)
        
        # Test first q-vector
        np.testing.assert_array_equal(data.irvec_q[0], [-1, -1, -1])
        self.assertAlmostEqual(data.wslen_q[0], 1.732051, places=5)
        
        # Test first g-vector
        np.testing.assert_array_equal(data.irvec_g[0], [-1, -1, -1])
        self.assertAlmostEqual(data.wslen_g[0], 1.732051, places=5)

    def test_round_trip_file_io(self):
        """Test reading and writing data preserves all information."""
        # Read original data
        data1 = WignerData.from_file(self.test_file)
        
        # Write to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fmt', delete=False) as tmp_file:
            temp_filename = tmp_file.name
        
        try:
            data1.to_file(temp_filename)
            
            # Read back the data
            data2 = WignerData.from_file(temp_filename)
            
            # Compare all fields
            self.assertEqual(data1.dims, data2.dims)
            self.assertEqual(data1.dims2, data2.dims2)
            self.assertEqual(data1.nrr_k, data2.nrr_k)
            self.assertEqual(data1.nrr_q, data2.nrr_q)
            self.assertEqual(data1.nrr_g, data2.nrr_g)
            
            # Compare arrays
            np.testing.assert_array_equal(data1.irvec_k, data2.irvec_k)
            np.testing.assert_array_equal(data1.ndegen_k, data2.ndegen_k)
            np.testing.assert_allclose(data1.wslen_k, data2.wslen_k)
            
            np.testing.assert_array_equal(data1.irvec_q, data2.irvec_q)
            np.testing.assert_array_equal(data1.ndegen_q, data2.ndegen_q)
            np.testing.assert_allclose(data1.wslen_q, data2.wslen_q)
            
            np.testing.assert_array_equal(data1.irvec_g, data2.irvec_g)
            np.testing.assert_array_equal(data1.ndegen_g, data2.ndegen_g)
            np.testing.assert_allclose(data1.wslen_g, data2.wslen_g)
            
        finally:
            # Clean up temporary file
            if os.path.exists(temp_filename):
                os.unlink(temp_filename)

    def test_summary_method(self):
        """Test the summary method returns expected format."""
        data = WignerData.from_file(self.test_file)
        summary = data.summary()
        
        self.assertIn("WignerData Summary:", summary)
        self.assertIn("1 Wannier functions", summary)
        self.assertIn("1 atoms", summary)
        self.assertIn("27 WS vectors", summary)

    def test_utility_methods(self):
        """Test utility methods work correctly."""
        data = WignerData.from_file(self.test_file)
        
        # Test find_origin_indices (this method still exists)
        k_orig, q_orig, g_orig = data.find_origin_indices()
        self.assertGreaterEqual(k_orig, 0)
        self.assertGreaterEqual(q_orig, 0)
        self.assertGreaterEqual(g_orig, 0)
        
        # Verify the origin vectors are actually [0, 0, 0]
        np.testing.assert_array_equal(data.irvec_k[k_orig], [0, 0, 0])
        np.testing.assert_array_equal(data.irvec_q[q_orig], [0, 0, 0])
        np.testing.assert_array_equal(data.irvec_g[g_orig], [0, 0, 0])

    def test_shape_validation(self):
        """Test that shape validation works correctly."""
        # Create data with correct shapes
        dims, dims2 = 2, 3
        nrr_k, nrr_q, nrr_g = 5, 7, 9
        
        data = WignerData(
            dims=dims, dims2=dims2,
            nrr_k=nrr_k, 
            irvec_k=np.zeros((nrr_k, 3), dtype=int),
            ndegen_k=np.ones((nrr_k, dims, dims), dtype=int),  # Use ones for positive degeneracies
            wslen_k=np.zeros(nrr_k),
            nrr_q=nrr_q,
            irvec_q=np.zeros((nrr_q, 3), dtype=int),
            ndegen_q=np.ones((nrr_q, dims2, dims2), dtype=int),
            wslen_q=np.zeros(nrr_q),
            nrr_g=nrr_g,
            irvec_g=np.zeros((nrr_g, 3), dtype=int),
            ndegen_g=np.ones((nrr_g, dims, dims2), dtype=int),
            wslen_g=np.zeros(nrr_g)
        )
        
        # Should not raise any assertion errors
        self.assertEqual(data.dims, dims)
        self.assertEqual(data.dims2, dims2)

    def test_shape_validation_failure(self):
        """Test that shape validation catches incorrect shapes."""
        dims, dims2 = 2, 3
        nrr_k = 5
        
        # Test with wrong irvec_k shape
        with self.assertRaises(AssertionError):
            WignerData(
                dims=dims, dims2=dims2,
                nrr_k=nrr_k,
                irvec_k=np.zeros((3, nrr_k), dtype=int),  # Wrong shape (Fortran style)
                ndegen_k=np.ones((nrr_k, dims, dims), dtype=int),
                wslen_k=np.zeros(nrr_k),
                nrr_q=1,
                irvec_q=np.zeros((1, 3), dtype=int),
                ndegen_q=np.ones((1, dims2, dims2), dtype=int),
                wslen_q=np.zeros(1),
                nrr_g=1,
                irvec_g=np.zeros((1, 3), dtype=int),
                ndegen_g=np.ones((1, dims, dims2), dtype=int),
                wslen_g=np.zeros(1)
            )

    def test_down_spin_data(self):
        """Test reading down-spin data if available."""
        down_file = "test/down/wigner.down.fmt"
        if os.path.exists(down_file):
            data_down = WignerData.from_file(down_file)
            data_up = WignerData.from_file(self.test_file)
            
            # Should have same structure
            self.assertEqual(data_up.dims, data_down.dims)
            self.assertEqual(data_up.dims2, data_down.dims2)
            self.assertEqual(data_up.nrr_k, data_down.nrr_k)
            self.assertEqual(data_up.nrr_q, data_down.nrr_q)
            self.assertEqual(data_up.nrr_g, data_down.nrr_g)


if __name__ == "__main__":
    unittest.main(verbosity=2)