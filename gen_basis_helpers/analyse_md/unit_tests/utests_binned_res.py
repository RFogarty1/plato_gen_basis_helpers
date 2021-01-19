

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

