

import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp


import gen_basis_helpers.adsorption.adsorbate_detectors as adsDetectHelp
import gen_basis_helpers.lammps_interface.pure_water_map_objs as tCode

class TestGetBondsFromInpGeom(unittest.TestCase):

	def setUp(self):
		self.minOHDist = 0.1
		self.maxOHDist = 2.0
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = _loadSimpleCubicCellWithWaterDimerA()
		kwargs = {"minBondLength":self.minOHDist, "maxBondLength":self.maxOHDist}
		adsDetector = adsDetectHelp.DetectH2OAdsorbatesFromInpGeomStandard(**kwargs)
		self.testObjA = tCode.GetPureWaterBondsFromInpGeomCentralCellOnly(adsDetector)

	def testExpectedForSimpleDimer(self):
		#NOTE: order is ID, type, atom1, atom2
		expBonds = [ [1, 1, 1, 2],
		             [2, 1, 1, 3],
		             [3, 1, 4, 5],
		             [4, 1, 4, 6] ]
		actBonds = self.testObjA(self.cellA)
		self.assertEqual(expBonds, actBonds)

	def testRaisesForSimpleDimerWhenBondsGreaterThanFour(self):
		self.maxOHDist = 12
		self.createTestObjs()
		with self.assertRaises(AssertionError):
			self.testObjA(self.cellA)

class TestGetAnglesFromInpGeom(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = _loadSimpleCubicCellWithWaterDimerA()
		adsDetector = _loadStandardWaterAdsDetectorA()
		self.testObjA = tCode.GetPureWaterAnglesFromInpGeomCentralCellOnly(adsDetector)

	def testExpectedForSimpleDimer(self):
		#ID, type, atomA, atomB, atomC (NOTE: for h-o-h bond O has to be the atomB integer)
		expAngles = [ [1, 1, 2, 1, 3],
		              [2, 1, 5, 4, 6] ]
		actAngles = self.testObjA(self.cellA)
		self.assertEqual(expAngles,actAngles)


class TestGetMoleculeIDListFromInpGeom(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = _loadSimpleCubicCellWithWaterDimerA()
		adsDetector = _loadStandardWaterAdsDetectorA()
		self.testObjA = tCode.GetPureWaterMoleculeIndicesFromInpGeomCentralCellOnly(adsDetector)

	def testExpectedForSimpleDimer(self):
		expVals = [1,1,1,2,2,2]
		actVals = self.testObjA(self.cellA)
		self.assertEqual(expVals, actVals)

	def testExpectedWithCartCoordsReordered(self):
		#Mess with the co-ordinate order
		cartCoords = self.cellA.cartCoords
		cartCoords[1], cartCoords[4] = cartCoords[4], cartCoords[1]
		cartCoords[2], cartCoords[5] = cartCoords[5], cartCoords[2]
		self.cellA.cartCoords = cartCoords

		#Test we get expected values
		expVals = [1,2,2,2,1,1]
		actVals = self.testObjA(self.cellA)
		self.assertEqual(expVals, actVals)


def _loadSimpleCubicCellWithWaterDimerA():
	lattParams = [10,10,10]
	lattAngles = [90,90,90]
	cartCoordsA = [ [1  , 1  , 1  , "O"],
	                [1  , 1.5, 1.5, "H"],
	                [1  , 0.5, 1.5, "H"],
	                [6  , 1  , 1  , "O"],
	                [6  , 1.5, 1.5, "H"],
	                [6  , 0.5, 1.5, "H"] ]
	outCell = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	outCell.cartCoords = cartCoordsA
	return outCell


def _loadStandardWaterAdsDetectorA(minOHDist=None, maxOHDist=None):
	minOHDist = 0.1 if minOHDist is None else minOHDist
	maxOHDist = 2.0 if maxOHDist is None else maxOHDist
	kwargs = {"minBondLength":minOHDist, "maxBondLength":maxOHDist}
	return adsDetectHelp.DetectH2OAdsorbatesFromInpGeomStandard(**kwargs)




