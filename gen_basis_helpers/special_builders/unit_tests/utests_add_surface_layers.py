
import copy
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.special_builders.add_surface_layers as tCode

class TestGetExtraSurfacePlaneOfAtomsByCopyAndTransAdjacentPlane(unittest.TestCase):

	def setUp(self):
		self.lattParamsA = [10,10,10]
		self.lattAnglesA = [90,90,90]

		self.surfCoordsA = [  [2,2,2,"A"],
		                      [3,3,2,"A"],
		                      [5,5,3,"B"],
		                      [6,6,3,"B"],
		                      [7,7,4,"C"],
		                      [8,8,4,"C"] ]
		self.top = True
		self.surfEles = None
		self.distTol = 1e-1
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.cellA.cartCoords = self.surfCoordsA

	def _runTestFunct(self):
		return tCode.getExtraSurfacePlaneOfAtomsByCopyAndTransAdjacentPlane(self.cellA)

	def testExpectedSimpleCaseA(self):
		expExtraCoords = [ [5,5,5,"B"],
		                   [6,6,5,"B"] ]
		actExtraCoords = self._runTestFunct()
		self._checkExpAndActCoordsEqual(expExtraCoords,actExtraCoords)

	def _checkExpAndActCoordsEqual(self, expCoords, actCoords):
		expCell, actCell = copy.deepcopy(self.cellA), copy.deepcopy(self.cellA)
		expCell.cartCoords = expCoords
		actCell.cartCoords = actCoords
		self.assertEqual(expCell,actCell)



