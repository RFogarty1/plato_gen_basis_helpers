

import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.ads_sites_impl as adsImplHelp

import gen_basis_helpers.analyse_md.ads_sites_hcp0001 as tCode

class TestGetBridgeSiteObjs(unittest.TestCase):

	def setUp(self):
		self.surfIndices = [1,4,16]
		self.maxDist = None #To limit problems when using subset of the surface indices
		self.distTol = 1 #Nearest nebs are those within this dist of the central atom
		self.siteNames = "bridge_test"
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = loadTestCellA()

	def _runTestFunct(self):
		args = [self.cellA, self.surfIndices]
		kwargs = {"distTol":self.distTol, "maxDist":self.maxDist, "siteNames":self.siteNames}
		return tCode.getHcp0001BridgeAdsSites(*args, **kwargs)

	#Super limited for simplicity
	def testExpectedCaseA(self):
		expIdxPairs = [ [1,4], [4,16] ]
		expObjs = [adsImplHelp.BridgeStandard(idxPair, siteName=self.siteNames) for idxPair in expIdxPairs]
		actObjs = self._runTestFunct()
		self.checkExpectedAndActualObjsMatch(expObjs,actObjs)

	def testNoneReturnedWithMaxDistSetSmall(self):
		self.maxDist = 5
		self.createTestObjs()
		expObjs = list()
		actObjs = self._runTestFunct()
		self.checkExpectedAndActualObjsMatch(expObjs, actObjs)

	def testExpected_largeDistTol(self):
		self.distTol = 20
		expIdxPairs = [ [1,4], [1,16], [4,16] ]
		expObjs = [adsImplHelp.BridgeStandard(idxPair, siteName=self.siteNames) for idxPair in expIdxPairs]
		actObjs = self._runTestFunct()
		self.checkExpectedAndActualObjsMatch(expObjs, actObjs)

#This is the distance matrix
#array([[ 0.        ,  6.06603193, 10.50664412],
#      [ 6.06603193,  0.        ,  6.06601334],
#       [10.50664412,  6.06601334,  0.        ]])

	def testExpected_largeDistTolAndShortMaxDist(self):
		self.distTol = 20
		self.maxDist = 7
		expIdxPairs = [ [1,4], [4,16] ]
		expObjs = [adsImplHelp.BridgeStandard(idxPair, siteName=self.siteNames) for idxPair in expIdxPairs]
		actObjs = self._runTestFunct()
		self.checkExpectedAndActualObjsMatch(expObjs, actObjs)


	def checkExpectedAndActualObjsMatch(self, expObjs, actObjs):
		sortedExpected = sorted(expObjs, key=lambda x:x.atomIndices)
		sortedActual = sorted(actObjs, key=lambda x:x.atomIndices)
		self.assertEqual(sortedExpected, sortedActual)


def loadTestCellA():
	cellVects = [ [18.19806, 0, 0],
	              [-9.099031, 15.75998, 0],
	              [0, 0, 54.76426] ]

	fractCoords = [
	               [0.000943, -0.000800, 0.590575, "Mg"],
	               [0.112937, 0.220712, 0.680002, "Mg"],
	               [0.334276, -0.000800, 0.590575, "Mg"],
	               [0.667610, -0.000800, 0.590575, "Mg"],
	               [0.446271, 0.220712, 0.680002, "Mg"],
	               [0.779604, 0.220712, 0.680002, "Mg"],
	               [0.000943, 0.332533, 0.590575, "Mg"],
	               [0.000943, 0.665867, 0.590575, "Mg"],
	               [0.112937, 0.554045, 0.680002, "Mg"],
	               [0.112937, 0.887379, 0.680002, "Mg"],
	               [0.334276, 0.332533, 0.590575, "Mg"],
	               [0.334276, 0.665867, 0.590575, "Mg"],
	               [0.667610, 0.332533, 0.590575, "Mg"],
	               [0.667610, 0.665867, 0.590575, "Mg"],
	               [0.446271, 0.554045, 0.680002, "Mg"],
	               [0.446271, 0.887379, 0.680002, "Mg"],
	               [0.779604, 0.554045, 0.680002, "Mg"],
	               [0.779604, 0.887379, 0.680002, "Mg"],
	               [0.890829, 0.776179, 0.713129, "O"],
	               [0.890775, 0.776191, 0.747008, "H"],
	               [0.224162, 0.776179, 0.713129, "O"],
	               [0.224109, 0.776191, 0.747008, "H"],
	               [0.557495, 0.776179, 0.713129, "O"],
	               [0.557442, 0.776191, 0.747008, "H"],
	               [0.890829, 0.109512, 0.713129, "O"],
	               [0.890775, 0.109525, 0.747008, "H"],
	               [0.890828, 0.442845, 0.713129, "O"],
	               [0.890775, 0.442858, 0.747008, "H"],
	               [0.224162, 0.109512, 0.713129, "O"],
	               [0.224109, 0.109525, 0.747008, "H"],
	               [0.224162, 0.442845, 0.713129, "O"],
	               [0.224109, 0.442858, 0.747008, "H"],
	               [0.557495, 0.109512, 0.713129, "O"],
	               [0.557442, 0.109525, 0.747008, "H"],
	               [0.557495, 0.442846, 0.713129, "O"],
	               [0.557442, 0.442858, 0.747008, "H"],
	              ]

	outCell = uCellHelp.UnitCell.fromLattVects( cellVects )
	outCell.fractCoords = fractCoords
	return outCell



