
import os
import copy
import itertools as it
import unittest
import unittest.mock as mock

import numpy as np

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

	def testEqualCompareEqual_iterOfIterBinVals(self):
		self.valsA = [ [x] for x in self.valsA ]
		self.createTestObjs()
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

	def testUnequalCompareUnequal_diffLenIterOfIters(self):
		self.valsA = [ [x] for x in self.valsA ]
		self.createTestObjs()
		objA = copy.deepcopy(self.testObjA)
		self.valsA[0].append(2)
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)

	def testUnequalCompreUnequal_diffValsInIterOfIters(self):
		self.valsA = [ [x] for x in self.valsA ]
		self.createTestObjs()
		objA = copy.deepcopy(self.testObjA)
		self.valsA[0][0] += 2
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)		

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


class TestAddProbabilityDensitiesToSimpleBin(unittest.TestCase):

	def setUp(self):
		self.binEdges = [1,2,4]
		self.counts = [4,6]
		self.countKey = "counts"
		self.probabilityKey = "pdf"
		self.createTestObjs()

	def createTestObjs(self):
		self.testObj = tCode.BinnedResultsStandard.fromBinEdges(self.binEdges)
		self.testObj.binVals["counts"] = self.counts

	def _runTestFunct(self):
		currKwargs =  {"countKey":self.countKey, "probabilityKey":self.probabilityKey}
		tCode.addProbabilityDensitiesToOneDimBinObj(self.testObj, **currKwargs)

	def testExpectedCaseA(self):
		expVals = [ 4/10, 3/10 ] #Dont need to add to 1 since their DENSITIES
		self._runTestFunct()
		actVals = self.testObj.binVals["pdf"]
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]


class TestBinCountResultsForOneDimData(unittest.TestCase):

	def setUp(self):
		self.dataA = [ 1, 2, 5, 2, 2]
		self.binEdges = [0, 3, 6]
		self.countLabel = "counts"
		self.raiseIfValsOutsideBins = True
		self.initIfNeeded = True
		self.binVals = None
		self.createTestObjs()

	def createTestObjs(self):
		self.binResObjA = tCode.BinnedResultsStandard.fromBinEdges(self.binEdges,binVals=self.binVals)

	def _runTestFunct(self):
		args = [self.dataA, self.binResObjA]
		currKwargs = {"countKey":self.countLabel, "raiseIfValsOutsideBins":self.raiseIfValsOutsideBins,
		              "initIfNeeded":self.initIfNeeded}
		return tCode.binCountsFromOneDimDataSimple(*args,**currKwargs)

	def testExpectedResultsA_uninitialised(self):
		expResObj = copy.deepcopy(self.binResObjA)
		expResObj.binVals = {self.countLabel: [4,1] }
		self._runTestFunct()
		actResObj = self.binResObjA
		self.assertEqual(expResObj, actResObj)

	def testExpectedResultsA_initialised(self):
		self.binVals = {self.countLabel:[10,12]}
		self.createTestObjs()
		expResObj = copy.deepcopy(self.binResObjA)
		expResObj.binVals = {self.countLabel:[14,13]}
		self._runTestFunct()
		actResObj = self.binResObjA
		self.assertEqual(expResObj, actResObj)

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

	def testExpectedWhenInpDataUnsorted(self):
		self.dataA = [ [3,9], [1,2], [4,16], [2,4] ]
		expBinVals = {self.dimZeroLabel: [ [1,2,3], [4] ], self.dimOneLabel: [ [2,4,9], [16] ] }
		self._runTestFunct()
		actBinVals = self.binResA.binVals
		self.assertEqual(expBinVals, actBinVals)


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

#Mostly covered by "TestCreateEmptyBinsStandard" really
class TestCreateBinsFromMinMaxAndWidthStandard(unittest.TestCase):

	def setUp(self):
		self.minVal = 4
		self.maxVal = 7
		self.width = 1
		self.extremesAtCentre = True

	def _runTestFunct(self):
		args = [self.minVal, self.maxVal, self.width]
		kwargs = {"extremesAtCentre":self.extremesAtCentre}
		return tCode.getEmptyBinResultsFromMinMaxAndWidthStandard(*args,**kwargs)

	def testExpectedCaseA(self):
		expEdges = [3.5, 4.5, 5.5, 6.5, 7.5]
		expObj = tCode.BinnedResultsStandard.fromBinEdges(expEdges)
		actObj = self._runTestFunct()
		self.assertEqual(expObj, actObj)

	def testExpected_extremesAtEdges(self):
		self.extremesAtCentre = False
		expEdges = [4,5,6,7]
		expObj = tCode.BinnedResultsStandard.fromBinEdges(expEdges)
		actObj = self._runTestFunct()
		self.assertEqual(expObj, actObj)

class TestConvertBinListsIntoAverages(unittest.TestCase):

	def setUp(self):
		self.keys = ["xVals", "yVals"]
		self.inpKeys = copy.deepcopy(self.keys)
		self.binEdges = [1,2,3] #Means two bins
		self.binValsKeyA = [ [5,7,12], [9,11,16] ] 
		self.binValsKeyB = [ [3,5,10], [8,10,12] ]

		self.createTestObjs()

	def createTestObjs(self):
		binVals = {self.keys[0]: self.binValsKeyA,
		           self.keys[1]: self.binValsKeyB}
		self.binObjA = tCode.BinnedResultsStandard.fromBinEdges(self.binEdges, binVals=binVals)

	def _runTestFunct(self):
		return tCode.averageItersInEachBin(self.binObjA, keys=self.inpKeys)

	def testExpectedValsA(self):
		expBinVals = {"xVals": [8, 12], "yVals":[6, 10]}
		expBinObj = tCode.BinnedResultsStandard.fromBinEdges(self.binEdges, binVals=expBinVals)
		actBinObj = self.binObjA
		self._runTestFunct()
		self.assertEqual(expBinObj, actBinObj)

	def testExpected_noKeysPassed(self):
		self.inpKeys = None
		expBinVals = {"xVals": [8, 12], "yVals":[6, 10]}
		expBinObj = tCode.BinnedResultsStandard.fromBinEdges(self.binEdges, binVals=expBinVals)
		actBinObj = self.binObjA
		self._runTestFunct()
		self.assertEqual(expBinObj, actBinObj)

	def testExpected_singleKeyPassed(self):
		self.inpKeys = ["yVals"]
		expBinVals = {"yVals":[6, 10], "xVals":self.binValsKeyA}
		expBinObj = tCode.BinnedResultsStandard.fromBinEdges(self.binEdges, binVals=expBinVals)
		actBinObj = self.binObjA
		self._runTestFunct()
		self.assertEqual(expBinObj, actBinObj)

	def testNoEffectIfBinValsNotIterOfIters(self):
		self.binValsKeyA = [1,2]
		self.binValsKeyB = [3,4]
		self.createTestObjs()
		expBin = copy.deepcopy(self.binObjA)
		self._runTestFunct()
		actBin = self.binObjA
		self.assertEqual(expBin, actBin)


#Mostly covered by average vals function; since they use the same moments-calculating backend function
class TestGetHigherMomentsForPdf(unittest.TestCase):

	def setUp(self):
		self.betweenVals = None
		self.normaliseBySum = False

		#Pdf vals come from a Gaussian with a mean of 4 and a variance of 0.8
		#Since we sample it so sparsely the calculated mean/variance wont be so neat but.....
		#Note: Mean comes out at about 4.02
		self.binEdges = [0, 2, 3, 5, 8] #Centres are 1,2.5,4,6.5
		self.pdfVals = [0.001608639066848, 0.109304604496336, 0.446031029038193, 0.008972268309668]


	def testExpectedVariance(self):
		expVal = 0.447706606977464 #Just calculated in excel
		currKwargs = {"betweenVals":self.betweenVals, "normaliseBySum":self.normaliseBySum}
		actVal = tCode.getVarianceForPdfValsSimple(self.binEdges, self.pdfVals, **currKwargs)
		self.assertAlmostEqual(expVal, actVal)

	def testExpectedSkew_noRange(self):
		expVal = -0.061496122223159 #Just calculated in excel, would be zero if sampling was better but...
		currKwargs = {"betweenVals":self.betweenVals, "normaliseBySum":self.normaliseBySum}
		actVal = tCode.getSkewForPdfValsSimple(self.binEdges, self.pdfVals, **currKwargs)
		self.assertAlmostEqual(expVal, actVal)

class TestGetAverageValsForPdf(unittest.TestCase):

	def setUp(self):
		self.betweenVals = None
		self.normaliseSum = False
		self.centreInBetweenVals = False

		#Edges are [0,2,3,5,8]
		self.binEdges =   [0, 2, 3, 5, 8]
		self.pdfVals = [0.4, 0.2, 0.3, 0.6] #Dont care how physical these are


	def _runTestFunct(self):
		args = [self.binEdges, self.pdfVals]
		currKwargs = {"betweenVals":self.betweenVals,"normaliseBySum":self.normaliseSum, "useCentreForBetweenVals":self.centreInBetweenVals}
		return tCode.getAverageForPdfValsSimple(*args, **currKwargs)

	def testExpectedNoRange(self):
		expVal = 15.4
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testRaisesIfEdgesOutOfOrder(self):
		self.binEdges = [2,3,1,5,3.1]
		with self.assertRaises(AssertionError):
			self._runTestFunct()

	def testExpectedBetweenVals(self):
		""" Note both edges of each bin need to be in range """
		self.betweenVals = [0.5, 5.5]
		expVal = 2.9
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testExpectedAddToTotal_binCentresAndEdgesStraddle(self):
		""" Test that we get the total sum from segments of betweenVals when using bin-centres to decide if bin is between vals """
		self.centreInBetweenVals = True
		self.betweenValsA, self.betweenValsB = [0,2.5], [2.5,10]
		expVals = [0.4*2, (2.5*0.2) + (4*0.3*2) + (6.5*0.6*3)] #This is what we'd get for no betweenVals limits
	
		#Run funct	
		actVals = []
		for betweenVals in [self.betweenValsA, self.betweenValsB]:
			self.betweenVals = betweenVals
			actVals.append( self._runTestFunct() )

		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	def testExpectedNormaliseWithBetweenValsSet(self):
		""" Make sure normalise factor uses only relevant pdf when betweenVals is set """
		self.betweenVals = [0.5,5.5]
		self.normaliseSum = True
		expVal = 2.9/0.8 #Note we use divide by the total probability rather than sum of densities 
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

class TestGetEdgesFromCentres_fixedWidth(unittest.TestCase):

	def setUp(self):
		self.centres = [1,3,5]

	def _runTestFunct(self):
		return tCode.getBinEdgesFromCentresFixedWidthAssumed(self.centres)

	def testExpectedCaseA(self):
		expVals = [0,2,4,6]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals, actVals)]

	def testRaisesWhenCentresOutOfOrder(self):
		self.centres = [x for x in reversed(self.centres)]
		with self.assertRaises(AssertionError):
			self._runTestFunct()

	def testRaisesWhenEqualWidthsImpossible(self):
		self.centres = [1,3,10]
		with self.assertRaises(AssertionError):
			self._runTestFunct()


#Same backend used as TestGetAverageValsForPdf; hence why we limit number of tests here
class TestIntegrateOverPdfSimple(unittest.TestCase):

	def setUp(self):
		self.betweenVals = None
		self.normByFullSum  = False
		self.binEdges =   [0, 2, 3, 5, 8]
		self.pdfVals = [0.4, 0.2, 0.3, 0.6] #Dont care how physical these are

	def _runTestFunct(self):
		args = [self.binEdges, self.pdfVals]
		currKwargs = {"betweenVals":self.betweenVals, "normByFullSum":self.normByFullSum}
		return tCode.getIntegralOverPdfSimple(*args, **currKwargs)

	def testExpectedCaseA(self):
		expVal = 3.4
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testExpectedCaseB_betweenVals(self):
		self.betweenVals = [0.5,5.5]
		expVal = 1.6
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal,actVal)

	def testExpected_betweenVals_normByTotal(self):
		self.normByFullSum = True
		self.betweenVals = [0.5,5.5]
		expVal = 1.6/3.4
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)


class TestIntegrateOverRdf(unittest.TestCase):

	def setUp(self):
		self.betweenVals = None
		self.prefactor = 5 #N_tot / V_tot generally gonna be the formula
		self.binEdges =   [0, 2, 3, 5, 8]
		self.rdfVals = [0.4, 0.2, 0.3, 0.6] #Dont care how physical these are

		#NOTE: Expected Surface areas are:
#		12.5663706143592, 78.5398163397448, 201.061929829747, 530.929158456675

	def _runTestFunct(self):
		currArgs = [self.binEdges, self.rdfVals]
		currKwargs = {"betweenVals":self.betweenVals, "prefactor":self.prefactor}
		return tCode.getWeightedIntegralOverRdf(*currArgs, **currKwargs)

	def testExpectedCaseA(self):
		expVal = 5510.3535143965
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testExpected_betweenVals(self):
		self.betweenVals = [-1,6]
		expVal = 731.991088286422
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

class TestIntegralOverCircularRdf(unittest.TestCase):

	def setUp(self):
		self.betweenVals = None
		self.prefactor = 5 #N_tot / V_tot generally gonna be the formula
		self.binEdges =   [0, 2, 3, 5, 8]
		self.rdfVals = [0.4, 0.2, 0.3, 0.6] #Dont care how physical these are

		#NOTE: Expected circumferences are:
		#6.28318530717959, 15.707963267949, 25.1327412287183,40.8407044966673

	def _runTestFunct(self):
		currArgs = [self.binEdges, self.rdfVals]
		currKwargs = {"betweenVals":self.betweenVals, "prefactor":self.prefactor}
		return tCode.getWeightedIntegralOverCircularRdf(*currArgs, **currKwargs)

	def testExpectedCaseA(self):
		expVal = 483.805268652828
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)





