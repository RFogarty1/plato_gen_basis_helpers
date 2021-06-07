
import copy
import math
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.ads_sites_impl as adsSiteImplHelp
import gen_basis_helpers.analyse_md.ads_sites_utils as tCode


class TestAssignAdsIndicesToAdsorptionSites(unittest.TestCase):

	def setUp(self):
		self.topIndices = [1,2]
		self.adsIndices = [0]
		self.maxHozDist = 4 #If somethings further than this its considered to be "None"
		self.siteName = "site_a"
		self.maxTotDist = None
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = loadTestCellA()
		self.adsObjs = [adsSiteImplHelp.TopStandard(x, siteName=self.siteName) for x in self.topIndices]

	def _runTestFunct(self):
		args = [self.cellA, self.adsObjs, self.adsIndices]
		kwargs = {"maxHozDist":self.maxHozDist, "maxTotDist":self.maxTotDist}
		return tCode.assignAdsIndicesToAdsorptionSites(*args, **kwargs)

	def testExpectedCaseA(self):
		expDict = {"None":list(), "site_a":[0]}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)


	def testExpectedCaseB(self):
		self.maxHozDist = 2 #sqrt(8) is the hoz dist for the None site
		self.topIndices = [0]
		self.adsIndices = [1,2]
		self.createTestObjs()
		expDict = {"None":[2], "site_a":[1]}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)

	def testExpectedTooFarFromHozDist(self):
		self.maxHozDist = 1
		expDict = {"None":[0], "site_a":list()}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)

	def testExpectedWithMaxTotalDist_closeEnoughToBoth(self):
		self.maxTotDist = 10
		expDict = {"None":list(), "site_a":[0]}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)

	def testExpectedWithMaxTotalDist_tooFarFromBoth(self):
		self.maxTotDist = math.sqrt(2) #Closest is sqrt(3)
		expDict = {"None":[0], "site_a":list()}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)

	def testExpectedWhenPbcsMatter(self):
		cartCoords = [ [9,9,9,"X"],
		               [1,1,1,"Y"],
		               [6,6,5,"Z"] ]
		self.cellA.cartCoords = cartCoords
		self.maxHozDist = math.sqrt(8) + 0.1 #Index 1 should work here, if PBCs accounted for
		expDict = {"None":list(), "site_a":[0]}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)


class TestAddAdsSitesStandard(unittest.TestCase):

	def setUp(self):
		self.topIndices = [1,2] 
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = loadTestCellA()
		self.topSites = [adsSiteImplHelp.TopStandard(idx,siteName="top") for idx in self.topIndices]
		self.testObjA = tCode.AddAdsSitesToGeomsStandard(self.topSites)

	def _runTestFunct(self):
		args = self.cellA
		return self.testObjA.addToUnitCell(self.cellA)

	def testAddToUnitCellSimple(self):
		expCoords = [ [4,4,6,"X"],
		              [5,5,7,"Y"],
		              [6,6,8,"Z"],
		              [5,5,7,"top"],
		              [6,6,8,"top"] ]
		self._runTestFunct()
		actCoords = self.cellA.cartCoords
		testVal = areExpAndActCoordsEqualUsingTemplateCell(expCoords, actCoords, self.cellA)
		self.assertTrue(testVal)


def areExpAndActCoordsEqualUsingTemplateCell(expCoords, actCoords, templCell):
	expCell, actCell = copy.deepcopy(templCell), copy.deepcopy(templCell)
	expCell.cartCoords = expCoords
	actCell.cartCoords = actCoords
	if expCell == actCell:
		return True
	else:
		return False

def loadTestCellA():
	lattParams, lattAngles = [10,10,10], [90,90,90]
	cartCoords = [ [4,4,6,"X"],
	               [5,5,7,"Y"],
	               [6,6,8,"Z"] ]

	outCell = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	outCell.cartCoords = cartCoords
	return outCell


