#!/usr/bin/python3

""" Test specific implementations of surface objects """

import copy
import unittest
import numpy as np
import plato_pylib.shared.ucell_class as UCell

import gen_basis_helpers.shared.surfaces as tCode

class TestHcp001Surface(unittest.TestCase):

	def setUp(self):
		baseLattAngles = [90.0, 90.0, 120.0]
		baseLengths = [3.21,3.21,5.91]
		baseFractCoords = None
		self.baseUnitCell = UCell.UnitCell(lattParams=baseLengths, lattAngles=baseLattAngles, fractCoords = baseFractCoords)
		self.vacuumRegion = 4.0
		self.nLayers = 1
		self.createSurfaceObj()

	def createSurfaceObj(self):
		self.testObj = tCode.Hcp0001Surface(self.baseUnitCell, self.nLayers, self.vacuumRegion)

	def testWrongAnglesLeadsToValueError(self):
		lattAngles = self.baseUnitCell.lattAngles
		lattAngles["gamma"] = 110.0 #Not actually hcp if this is the angle, so should break
		self.baseUnitCell.lattAngles = lattAngles
		with self.assertRaises(ValueError):
			self.createSurfaceObj()

	def testExpectedSurfaceAreaGiven(self):
		lattParams = self.baseUnitCell.lattParams
		lattParams["b"] = 2*lattParams["b"] #As it would be in a supercell
		self.baseUnitCell.lattParams = lattParams
		self.createSurfaceObj()
		expSurfaceArea = 17.847224726270476
		actSurfaceArea = self.testObj.surfaceArea
		self.assertAlmostEqual(expSurfaceArea, actSurfaceArea)

	def testExpCartCoordsGiven(self):
		""" For hcp0001 the z-coords of each atom should be shifted by 0.5* vacuum length. This is easy to check
		for a single layer """
		self.baseUnitCell.fractCoords = [ [0.0,0.0,0.0,"Mg"], [0.33,0.66,0.5,"Mg"] ]
		startCartCoords = copy.deepcopy(self.baseUnitCell.cartCoords)
		self.createSurfaceObj()

		#Form the expected co-ordinates
		expCoords = np.array( copy.deepcopy([x[:3] for x in startCartCoords]) )
		expCoords[0,2] += 0.5*self.vacuumRegion
		expCoords[1,2] += 0.5*self.vacuumRegion
		actCoords = np.array( [x[:3] for x in self.testObj.unitCell.cartCoords] )

		self.assertTrue( np.allclose(expCoords,actCoords) )

	def testExpNumbAtomsInMultiLayerCell(self):
		""" Test we get the expected number of atoms when we apply multiple layers """
		self.baseUnitCell.fractCoords = [ [0.0,0.0,0.0,"Mg"], [0.33,0.66,0.5,"Mg"] ]
		self.nLayers = 4
		self.createSurfaceObj()

		#Test the initially constructed
		expectedNumbAtoms = 8
		actNumbAtoms = len( self.testObj.unitCell.cartCoords )
		self.assertEqual(expectedNumbAtoms, actNumbAtoms)

		#Test it after i've changed nLayers post construction
		self.testObj._nLayers = 3
		expectedNumbAtoms = 6
		actNumbAtoms = len( self.testObj.unitCell.cartCoords )
		self.assertEqual(expectedNumbAtoms, actNumbAtoms)

if __name__ == '__main__':
	unittest.main()
