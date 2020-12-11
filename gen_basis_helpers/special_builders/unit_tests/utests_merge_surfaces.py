
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.special_builders.merge_surfaces as tCode

class TestGetSurfaceAndAdsCellsMergedStandard(unittest.TestCase):

	def setUp(self):
		self.lattParamsSurf , self.lattAnglesSurf = [4,4,4], [90,90,90]
		self.lattParamsAbove, self.lattAnglesAbove = [4,4,8], [90,90,90]
		self.lattParamsBelow, self.lattAnglesBelow = [4,4,2], [90,90,90]

		self.cartCoordsSurf =  [[2,2,2,"surf"], [2,2,3,"surf"] ]
		self.cartCoordsAbove = [[2,2,4,"sAbove"]]
		self.cartCoordsBelow = [[2,2,1,"sBelow"]]
		self.deltaVac = 0
		self.createTestObjs()

	def createTestObjs(self):
		self._createCells()

	def _createCells(self):
		self.cellSurf  = uCellHelp.UnitCell(lattParams=self.lattParamsSurf , lattAngles=self.lattAnglesSurf )
		self.cellAbove = uCellHelp.UnitCell(lattParams=self.lattParamsAbove, lattAngles=self.lattAnglesAbove)
		self.cellBelow = uCellHelp.UnitCell(lattParams=self.lattParamsBelow, lattAngles=self.lattAnglesBelow)

		self.cellSurf.cartCoords  = self.cartCoordsSurf
		self.cellAbove.cartCoords = self.cartCoordsAbove
		self.cellBelow.cartCoords = self.cartCoordsBelow 

	def _runTestFunct(self):
		return tCode.getGeomForSurfaceMergedWithCells(self.cellSurf, self.cellAbove, geomBelowSurface=self.cellBelow, deltaVac=self.deltaVac)

	def testRaisesForIncompatibleCellAbove_diffAngles(self):
		self.lattAnglesAbove[-1] = self.lattAnglesSurf[-1] + 2
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self._runTestFunct()

	def testRaisesForIncompatibleCellBelow_diffBParam(self):
		self.lattParamsBelow[1] = self.lattParamsSurf[1]+0.5
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self._runTestFunct()

	def testExpecetedForSimpleCubicAboveOnly(self):
		self.cellBelow = None
		expOutCell = self._getExpectedCellForAboveOnly()
		actOutCell = self._runTestFunct()
		self.assertEqual(expOutCell,actOutCell)

	def testExpecetedForSimpleCubicAboveOnly_cartAboveNotCentred(self):
		self.cartCoordsAbove = [[2,2,7,"sAbove"]]
		self.createTestObjs()
		self.cellBelow = None
		expOutCell = self._getExpectedCellForAboveOnly()
		actOutCell = self._runTestFunct()
		self.assertEqual(expOutCell,actOutCell)

	def testExpecetedForSimpleCubicAboveAndBelow(self):
		expOutCell = self._getExpectedCellForAboveAndBelow()
		actOutCell = self._runTestFunct()
		self.assertEqual(expOutCell, actOutCell)

	def testExpectedForSimpleCubicAboveBelowAndDeltaVac(self):
		self.deltaVac = 2
		expOutCell = self._getExpectedCellForAboveAndBelow()
		actOutCell = self._runTestFunct()
		self.assertEqual(expOutCell, actOutCell)

	def _getExpectedCellForAboveOnly(self):
		heightSurface = max([x[-2] for x in self.cartCoordsSurf]) - min ([x[-2] for x in self.cartCoordsSurf])
		heightCellAbove = self.lattParamsAbove[-1]
		totalHeight = heightSurface + heightCellAbove
		expLattParams = [4,4,heightSurface+heightCellAbove]
		#NOTE: By default we put the original surface at the centre of the cell
		topSurfPosZ = (totalHeight/2)+0.5
		expCartCoords = [ [2, 2, topSurfPosZ-1, "surf"],
		                  [2, 2, topSurfPosZ, "surf"],
		                  [2, 2, topSurfPosZ + 4,  "sAbove"] ]
		outCell = uCellHelp.UnitCell(lattParams=expLattParams, lattAngles=self.lattAnglesSurf)
		outCell.cartCoords = expCartCoords

		return outCell

	def _getExpectedCellForAboveAndBelow(self):
		heightSurface = max([x[-2] for x in self.cartCoordsSurf]) - min ([x[-2] for x in self.cartCoordsSurf])
		heightCellAbove = self.lattParamsAbove[-1]
		heightCellBelow = self.lattParamsBelow[-1]
		totalHeight = heightSurface + heightCellAbove + heightCellBelow + self.deltaVac
		expLattParams = [4,4,totalHeight]
		topSurfPosZ = (totalHeight/2) + heightSurface/2
		expCartCoords = [ [2, 2, topSurfPosZ-heightSurface, "surf"],
		                  [2, 2, topSurfPosZ              , "surf"],
		                  [2, 2, topSurfPosZ + 4          , "sAbove"],
		                  [2, 2, topSurfPosZ-heightSurface-1, "sBelow"] ]
		outCell = uCellHelp.UnitCell(lattParams=expLattParams, lattAngles=self.lattAnglesSurf)
		outCell.cartCoords = expCartCoords
		return outCell
	


