#!/usr/bin/python3
import unittest
import sys
sys.path.append("../")
import FortranDriver
import numpy as np

class DriverSmokeTest(unittest.TestCase):
    """ Set of "smoke tests" to check the driver for the Fortran backend runs at all.
    
    These tests are somewhat unsophisticated and are not strictly unit tests, 
    but provide value because if they fail then there's no point in continuing
    to test."""
    def setUp(self):
        """ Initialise the Fortran backend."""
        self.md = FortranDriver.MDInterface()
        self.md.setup()
    
    def test_md_short(self):
        """ Test a single step of the MD routine. """
        x_start, y_start, z_start = self.md.x.copy(), self.md.y.copy(), self.md.z.copy()
        px_start, py_start, pz_start = self.md.px.copy(), self.md.py.copy(), self.md.pz.copy()
        temp_start = self.md.temp

        self.md.run(1)

        x_end, y_end, z_end = self.md.x.copy(), self.md.y.copy(), self.md.z.copy()
        px_end, py_end, pz_end = self.md.px.copy(), self.md.py.copy(), self.md.pz.copy()
        # Particle position should be different to when we start
        self.assertFalse(np.all(x_start == x_end))
        self.assertFalse(np.all(y_start == y_end))
        self.assertFalse(np.all(z_start == z_end))

        # Particles should not be outside the box bounds (0 < x < md.box_bounds[0])
        # (x,y,z) >= 0
        self.assertTrue(np.all(self.md.x >= 0))
        self.assertTrue(np.all(self.md.y >= 0))
        self.assertTrue(np.all(self.md.z >= 0))
        # (x,y,z) < (xmax,ymax,zmax)
        self.assertTrue(np.all(self.md.x < self.md.box_bounds[0]))
        self.assertTrue(np.all(self.md.y < self.md.box_bounds[0]))
        self.assertTrue(np.all(self.md.z < self.md.box_bounds[0]))
        # Basic thermodynamic properties
        self.assertGreater(self.md.temp, 0.0)
        self.assertAlmostEqual(self.md.temp, temp_start)

        # Conservation of total momentum. This should hold because we initialise the momentum distribution
        # to the desired temperature, which we don't change during this test. Using assertAlmostEqual
        # to allow for finite floating point precision
        self.assertAlmostEqual(np.sum(px_start), np.sum(px_end))
        self.assertAlmostEqual(np.sum(py_start), np.sum(py_end))
        self.assertAlmostEqual(np.sum(pz_start), np.sum(pz_end))

    def test_md_long(self):
        """ Test many steps of the MD routine"""
        x_start, y_start, z_start = self.md.x.copy(), self.md.y.copy(), self.md.z.copy()
        px_start, py_start, pz_start = self.md.px.copy(), self.md.py.copy(), self.md.pz.copy()
        temp_start = self.md.temp

        self.md.run(1000)

        x_end, y_end, z_end = self.md.x.copy(), self.md.y.copy(), self.md.z.copy()
        px_end, py_end, pz_end = self.md.px.copy(), self.md.py.copy(), self.md.pz.copy()
        # Particle position should be different to when we start
        self.assertFalse(np.all(x_start == x_end))
        self.assertFalse(np.all(y_start == y_end))
        self.assertFalse(np.all(z_start == z_end))

        # Particles should not be outside the box bounds (0 < x < md.box_bounds[0])
        # (x,y,z) >= 0
        self.assertTrue(np.all(self.md.x >= 0))
        self.assertTrue(np.all(self.md.y >= 0))
        self.assertTrue(np.all(self.md.z >= 0))
        # (x,y,z) < (xmax,ymax,zmax)
        self.assertTrue(np.all(self.md.x < self.md.box_bounds[0]))
        self.assertTrue(np.all(self.md.y < self.md.box_bounds[0]))
        self.assertTrue(np.all(self.md.z < self.md.box_bounds[0]))

        # Basic thermodynamic properties. Temperature can never go below zero, and shouldn't change too much
        # over the course of this run.
        self.assertGreater(self.md.temp, 0.0)
        self.assertAlmostEqual(self.md.temp, temp_start)

        # Conservation of total momentum. This should hold because we initialise the momentum distribution
        # to the desired temperature, which we don't change during this test. Using assertAlmostEqual
        # to allow for finite floating point precision
        self.assertAlmostEqual(np.sum(px_start), np.sum(px_end))
        self.assertAlmostEqual(np.sum(py_start), np.sum(py_end))
        self.assertAlmostEqual(np.sum(pz_start), np.sum(pz_end))

        # Ensure mean-squared displacement is nonzero
        msdx, msdy, msdz = self.md.compute_msd()
        self.assertGreater(msdx, 0.0)
        self.assertGreater(msdy, 0.0)
        self.assertGreater(msdz, 0.0)

if __name__ == '__main__':
    unittest.main()
