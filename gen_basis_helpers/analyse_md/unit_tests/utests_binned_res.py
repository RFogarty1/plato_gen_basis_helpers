

import copy
import unittest
import unittest.mock as mock

import gen_basis_helpers.analyse_md.binned_res as tCode

class TestBinnedResultsObj(unittest.TestCase):

	def setUp(self):
		self.binCentres = [0.5, 1.5, 2.5]
		self.binEdges = [0, 1, 2, 3]
		self.valsA = [2, 3, 4]
		self.valsB = [5, 6, 7]
		self.createTestObjs()

	def createTestObjs(self):
		self.binVals = {"valsA":self.valsA, "valsB":self.valsB}
		args = [self.binCentres, self.binEdges, self.binVals]
		self.testObjA = tCode.BinnedResultsStandard(*args)

	def testEqualCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA,objB)

	def testUnequalCompareUnequal_diffLenBins(self):
		objA = copy.deepcopy(self.testObjA)
		self.binCentres.append(1)
		self.binEdges.append(2)
		self.valsA.append(3)
		self.valsB.append(4)
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)
	
	def testUnequalCompareUnequal_diffBinCentreVals(self):
		objA = copy.deepcopy(self.testObjA)
		self.binCentres[-1] += 1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalCompareUnequal_diffBinValKeys(self):
		objA = copy.deepcopy(self.testObjA)
		self.testObjA.binVals["new_key"] = self.valsB
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalCompareUnequal_diffLenValsB(self):
		objA = copy.deepcopy(self.testObjA)
		self.valsA.append(2)
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalCompareUnequal_diffValsB(self):
		objA = copy.deepcopy(self.testObjA)
		self.valsA[-1] += 1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testFromConstantWidth(self):
		binWidth = 1
		expObj = self.testObjA
		actObj = tCode.BinnedResultsStandard.fromConstantBinWidth(self.binCentres, binWidth, self.binVals)
		self.assertEqual(expObj,actObj)

	def testFromBinEdges(self):
		binEdges = copy.deepcopy(self.testObjA.binEdges)
		expObj = copy.deepcopy(self.testObjA)
		actObj = tCode.BinnedResultsStandard.fromBinEdges(binEdges, binVals=self.binVals)
		self.assertEqual(expObj,actObj)


class TestBinResultsForTwoDimData(unittest.TestCase):

	def setUp(self):
		self.dataA = [ [1,2], [2,4], [3,9], [4,16] ]
		self.binEdges = [ 1, 3.01, 5 ]
		self.binVals = dict()

		self.dimZeroLabel, self.dimOneLabel = "dim_zero_label", "dim_one_label"
		self.raiseIfValsOutsideBins = False
		self.createTestObjs()

	def createTestObjs(self):
		self.binResA = tCode.BinnedResultsStandard.fromBinEdges(self.binEdges, binVals=self.binVals)

	def _runTestFunct(self):
		args = [self.dataA, self.binResA]
		kwargs = {"dimZeroLabel":self.dimZeroLabel, "dimOneLabel":self.dimOneLabel,
		          "raiseIfValsOutsideBins":self.raiseIfValsOutsideBins}
		return tCode.binResultsFromTwoDimDataSimple(*args, **kwargs)
	
	def testExpectedValsA(self):
		expBinVals = {self.dimZeroLabel: [ [1,2,3], [4] ], self.dimOneLabel: [ [2,4,9], [16] ] }
		self._runTestFunct()
		actBinVals = self.binResA.binVals
		self.assertEqual(expBinVals, actBinVals)

	def testAppendsWhenBinValsAlreadyContainsVals(self):
		self.binVals = {self.dimZeroLabel:[ [1], [2] ], self.dimOneLabel: [ [4], [5] ]}
		self.createTestObjs()
		expBinVals = {self.dimZeroLabel: [ [1,1,2,3], [2,4] ], self.dimOneLabel: [ [4,2,4,9], [5,16] ] }
		self._runTestFunct()
		actBinVals = self.binResA.binVals
		self.assertEqual(expBinVals, actBinVals)

	def testRaisesWhenValueNotInAnyBin(self):
		self.dataA.append( [20,5] )
		self.raiseIfValsOutsideBins = True
		with self.assertRaises(ValueError):
			self._runTestFunct()

class TestCreateEmptyBinsStandard(unittest.TestCase):

	def setUp(self):
		self.inpVals = [1,2,4,5]
		self.binWidth = 1

	def _runTestFunct(self):
		return tCode.getEmptyBinResultsForValsStandard(self.inpVals, self.binWidth)

	def testExpectedCaseA(self):
		expEdges = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
		expObj = tCode.BinnedResultsStandard.fromBinEdges(expEdges)
		actObj = self._runTestFunct()
		self.assertEqual(expObj, actObj)

	def testExpectedCaseB(self):
		self.inpVals = [1,2,4,5.1]
		expEdges = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
		expObj = tCode.BinnedResultsStandard.fromBinEdges(expEdges)
		actObj = self._runTestFunct()
		self.assertEqual(expObj, actObj)


