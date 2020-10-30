
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.adsorption.surface_detectors as tCode

class TestSurfaceDetectorBasedOnElement(unittest.TestCase):

	def setUp(self):
		self.lattParamsA = [10,10,10]
		self.lattAnglesA = [90,90,90]
		self.cartCoordsA = [ [5,5,5,"X"],
		                     [7,7,7,"Y"],
		                     [3,3,3,"X"],
		                     [2,2,2,"Z"] ]
		self.elesToInclude = ["X"]
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"lattParams":self.lattParamsA, "lattAngles":self.lattAnglesA}
		self.testCellA = uCellHelp.UnitCell(**kwargDict)
		self.testCellA.cartCoords = self.cartCoordsA
		self.detectorA = tCode.DetectSurfaceBasedOnElementsPresent( self.elesToInclude )

	def testExpectedCoordsBasedOnElement(self):
		expCoords = [ [5,5,5,"X"], [3,3,3,"X"] ]
		actCoords = self.detectorA(self.testCellA)
		self.checkCoordsEqual( expCoords, actCoords )

	def testExpectedCoordsFromTwoElements(self):
		self.elesToInclude = ["X","Y"]
		self.createTestObjs()
		expCoords = [ [5,5,5,"X"], [7,7,7,"Y"], [3,3,3,"X"] ]
		actCoords = self.detectorA(self.testCellA)
		self.checkCoordsEqual( expCoords, actCoords )

	def checkCoordsEqual(self, coordsA, coordsB):
		for cA, cB in it.zip_longest(coordsA,coordsB):
			self.assertEqual( len(cA), len(cB) )
			[self.assertAlmostEqual(e,a) for e,a in zip(cA[:3],cB[:3])]
			self.assertEqual( cA[-1], cB[-1] )
