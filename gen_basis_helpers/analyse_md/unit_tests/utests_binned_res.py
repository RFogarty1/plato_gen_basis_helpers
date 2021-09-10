
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


#Needs slightly different setup to the other tests for this class
class TestNDimensionalBinObj_equality(unittest.TestCase):

	def setUp(self):
		self.edges = [ [1,2,3],
		               [6,5,4,3] ]
		self.binVals = None
		self.createTestObjs()

	def createTestObjs(self):
		self.testObj = tCode.NDimensionalBinnedResults(self.edges, binVals=self.binVals)
		self.testObj.initialiseCountsMatrix()

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObj)
		self.createTestObjs()
		objB = self.testObj
		self.assertEqual(objA, objB)

	def testUnequalCompareUnequal_extraDimensionInOne(self):
		objA = copy.deepcopy(self.testObj)
		self.edges.append( [5,6] )
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testUnequalCompareUnequal_diffEdgeValues(self):
		objA = copy.deepcopy(self.testObj)
		self.edges[1][0] += 2
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testUnequalCompareUnequal_diffKeysInBinVals(self):
		#binVals should really only ever be np arrays; but shouldnt matter here
		objA = copy.deepcopy(self.testObj)
		self.binVals = {"any_key":"any_val"}
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testUnequalCompareUnequal_diffValsInCounts(self):
		objA = copy.deepcopy(self.testObj)
		objB = copy.deepcopy(self.testObj)

		objA.initialiseCountsMatrix(), objB.initialiseCountsMatrix()
		objA.binVals["counts"][0][0] += 1
		self.assertNotEqual(objA, objB)


class TestNDimensionalBinObj(unittest.TestCase):

	def setUp(self):
		self.edgesA = [1,2,3] 
		self.edgesB = [6,5,4,3]
		self.createTestObjs()

	def createTestObjs(self):
		edges = [ self.edgesA, self.edgesB ]
		self.testObj = tCode.NDimensionalBinnedResults(edges)
		self.testObj.initialiseCountsMatrix()

	def testBinEdgesArray_2dim(self):
		expArray = [ [None,None,None], [None,None,None] ]
		expArray[0][0], expArray[0][1]  = [ [1,2], [6,5] ], [ [1,2], [5,4] ]
		expArray[0][2], expArray[1][0]  = [ [1,2], [4,3] ], [ [2,3], [6,5] ]
		expArray[1][1], expArray[1][2]  = [ [2,3], [5,4] ], [ [2,3], [4,3] ]

		expArray = np.array(expArray)
		actArray = self.testObj.binEdgesArray

		self.assertTrue( np.allclose(expArray,actArray) )

	def testBinEdgesArray_3dim(self):
		self.edgesC = [7,8,9]
		edges = [self.edgesA, self.edgesB, self.edgesC]
		self.testObj = tCode.NDimensionalBinnedResults(edges)

		#2,3,2 are number of bins along first 3 dimensions.
		# Then we need 3 slots to hold EACH of the bin edge pairs
		# THEN we need another 2 slots to hold upper and lower values of the bin
		expArray = np.zeros( [2,3,2,len(edges),2] )

		expArray[0][0][0] = [ np.array([1,2]), np.array([6,5]), np.array([7,8]) ] 
		expArray[0][0][1] = [ np.array([1,2]), np.array([6,5]), np.array([8,9]) ]
		expArray[0][1][0] = [ np.array([1,2]), np.array([5,4]), np.array([7,8]) ]
		expArray[0][1][1] = [ np.array([1,2]), np.array([5,4]), np.array([8,9]) ]
		expArray[0][2][0] = [ np.array([1,2]), np.array([4,3]), np.array([7,8]) ]
		expArray[0][2][1] = [ np.array([1,2]), np.array([4,3]), np.array([8,9]) ]

		expArray[1][0][0] = [ np.array([2,3]), np.array([6,5]), np.array([7,8]) ]
		expArray[1][0][1] = [ np.array([2,3]), np.array([6,5]), np.array([8,9]) ]
		expArray[1][1][0] = [ np.array([2,3]), np.array([5,4]), np.array([7,8]) ]
		expArray[1][1][1] = [ np.array([2,3]), np.array([5,4]), np.array([8,9]) ]
		expArray[1][2][0] = [ np.array([2,3]), np.array([4,3]), np.array([7,8]) ]
		expArray[1][2][1] = [ np.array([2,3]), np.array([4,3]), np.array([8,9]) ]

		actArray = self.testObj.binEdgesArray

		self.assertTrue( np.allclose(expArray, actArray) )


	def testExpectedBinCentresArray_3dim(self):
		self.edgesC = [7,8,10] #Putting in a diff width one
		edges = [self.edgesA, self.edgesB, self.edgesC]
		self.testObj = tCode.NDimensionalBinnedResults(edges)

		expArray = np.zeros( [2,3,2,len(edges)] )

		#Only [X][X][1] are different from theyd be if using the matrix above to calculate
		expArray[0][0][0] = [ 1.5, 5.5, 7.5 ] 
		expArray[0][0][1] = [ 1.5, 5.5, 9 ]
		expArray[0][1][0] = [ 1.5, 4.5, 7.5 ]
		expArray[0][1][1] = [ 1.5, 4.5, 9 ]
		expArray[0][2][0] = [ 1.5, 3.5, 7.5 ]
		expArray[0][2][1] = [ 1.5, 3.5, 9 ]

		expArray[1][0][0] = [ 2.5, 5.5, 7.5 ]
		expArray[1][0][1] = [ 2.5, 5.5, 9 ]
		expArray[1][1][0] = [ 2.5, 4.5, 7.5 ]
		expArray[1][1][1] = [ 2.5, 4.5, 9 ]
		expArray[1][2][0] = [ 2.5, 3.5, 7.5 ]
		expArray[1][2][1] = [ 2.5, 3.5, 9 ]

		actArray = self.testObj.binCentresArray
		self.assertTrue( np.allclose(expArray,actArray) )

	def testAddValsToBinCounts(self):

		self.edgesA = [1,2,3] 
		self.edgesB = [6,5,4,3]

		valsToBin = [ [2.5, 3.5], #[1,2]
		              [1.2, 5.5], #[0,0]
		              [1.5,5.5], #[0,0] 
		              [1.5,5], # [0,0] Edge case
		              [2, 4.5] ] #Edge case; but its [1,1]

		expCountsMatrix = copy.deepcopy(self.testObj.binVals["counts"])
		expCountsMatrix[0][0] = 3
		expCountsMatrix[1][1] = 1
		expCountsMatrix[1][2] = 1

		self.testObj.addBinValuesToCounts( valsToBin )

		actCountsMatrix = self.testObj.binVals["counts"]

		self.assertTrue( np.allclose(expCountsMatrix,actCountsMatrix) )

	def testToAndFromDictConsistent(self):
		self.testObj.binVals["counts"][0][0] += 2
		outDict = self.testObj.toDict()
		expObj = self.testObj
		actObj = tCode.NDimensionalBinnedResults.fromDict(outDict)
		self.assertEqual(expObj, actObj)

class TestDumpAndReadNDimBinnedResultsJson(unittest.TestCase):

	def setUp(self):
		self.edgesA = [1,2,3] 
		self.edgesB = [6,5,4,3]
		self.tempFileName = "_tempDumpBinResFile.json"
		self.createTestObjs()

	def tearDown(self):
		os.remove(self.tempFileName)

	def createTestObjs(self):
		edges = [ self.edgesA, self.edgesB ]
		self.testObj = tCode.NDimensionalBinnedResults(edges)
		self.testObj.initialiseCountsMatrix()
		self.outIter = [self.testObj]

	def _dumpFile(self):
		tCode.dumpIterOfNDimBinnedResultsToJson(self.outIter, self.tempFileName)

	def _readFromFile(self):
		return tCode.readIterNDimensionalBinnedResFromJson(self.tempFileName)

	def testConsistentA(self):
		objA = copy.deepcopy(self.testObj)
		self.testObj.binVals["counts"][0][1] += 2
		objB = self.testObj

		expIter = [objA, objB]
		self.outIter = [objA, objB]

		self._dumpFile()
		actIter = self._readFromFile()
		self.assertEqual(expIter, actIter)


class TestGetLowerDimBinObj(unittest.TestCase):

	def setUp(self):
		#Define bins
		self.edgesA = [1,2,3]
		self.edgesB = [4,5,6]
		self.edgesC = [7,8,9]

		#Run options
		self.keepDims = [0,2]
		self.useIdxOther = [1]

		self.createTestObjs()

	def createTestObjs(self):
		self.comboEdges = [self.edgesA, self.edgesB, self.edgesC]
		self.binObjA = tCode.NDimensionalBinnedResults(self.comboEdges)

		#Put various values there
		#Use a weird key to make sure we dont rely on it being defualt val
		self.binObjA.initialiseCountsMatrix(countKey="fake_counts") 
		self.binObjA.binVals["fake_counts"][0][0][0] = 2
		self.binObjA.binVals["fake_counts"][0][1][1] = 3
		self.binObjA.binVals["fake_counts"][0][1][0] = 5
		self.binObjA.binVals["fake_counts"][1][0][1] = 7

	def _runTestFunct(self):
		args = [self.binObjA, self.keepDims, self.useIdxOther]
		return tCode.getLowerDimNDimBinObj_takeSingleBinFromOthers(*args)

	def _loadExpectedBinsCaseA(self):
		outObj = tCode.NDimensionalBinnedResults( [self.edgesA, self.edgesC] )
		outObj.initialiseCountsMatrix(countKey="fake_counts")
		outObj.binVals["fake_counts"][0][0] = 5
		outObj.binVals["fake_counts"][0][1] = 3
		return outObj

	def testExpectedCaseA_3dim_to_2dim(self):
		expBins = self._loadExpectedBinsCaseA()
		actBins = self._runTestFunct()
		self.assertEqual(expBins, actBins)

	def testExpectedCaseB_3dim_to_1dim(self):
		#Setup
		self.keepDims = [0]
		self.useIdxOther = [0,0]

		#Figure out what we expect
		expBinEdges = [self.edgesA]
		expBinCounts = np.array( [[2,0]] )
		expBinObj = tCode.NDimensionalBinnedResults( expBinEdges, binVals={"fake_counts":expBinCounts} )

		#Run and test
		actBinObj = self._runTestFunct()
		self.assertEqual(expBinObj, actBinObj)


	def testExpectedCaseSameBins(self):
		self.keepDims = [0,1,2]
		self.useIdxOther = list()
		expBins = copy.deepcopy(self.binObjA)
		actBins = self._runTestFunct()
		self.assertEqual(expBins, actBins)




class TestGetLowerDimObj_integrationMethod(unittest.TestCase):

	def setUp(self):
		#Define bins
		self.edgesA = [1,2,3]
		self.edgesB = [4,5,6]
		self.edgesC = [7,8,9]

		#Run options
		self.keepDims = [0,2]

		self.createTestObjs()

	def createTestObjs(self):
		self.comboEdges = [self.edgesA, self.edgesB, self.edgesC]
		self.binObjA = tCode.NDimensionalBinnedResults(self.comboEdges)

		#Put various values there
		#Use a weird key to make sure we dont rely on it being defualt val
		self.binObjA.initialiseCountsMatrix(countKey="fake_counts") 
		self.binObjA.binVals["fake_counts"][0][0][0] = 1
		self.binObjA.binVals["fake_counts"][0][0][1] = 2
		self.binObjA.binVals["fake_counts"][0][1][0] = 3
		self.binObjA.binVals["fake_counts"][0][1][1] = 4

		self.binObjA.binVals["fake_counts"][1][0][0] = 5
		self.binObjA.binVals["fake_counts"][1][0][1] = 6
		self.binObjA.binVals["fake_counts"][1][1][0] = 7
		self.binObjA.binVals["fake_counts"][1][1][1] = 8

	def _runTestFunct(self):
		return tCode.getLowerDimNDimBinObj_integrationMethod(self.binObjA, self.keepDims)

	def _loadExpectedBinsCaseA_3d_to_2d(self):
		outEdges = [self.edgesA, self.edgesC]
		outBinObj = tCode.NDimensionalBinnedResults(outEdges)

		#
		outBinObj.initialiseCountsMatrix(countKey="fake_counts")

		outBinObj.binVals["fake_counts"][0][0] = 1+3
		outBinObj.binVals["fake_counts"][0][1] = 2+4
		outBinObj.binVals["fake_counts"][1][0] = 5+7
		outBinObj.binVals["fake_counts"][1][1] = 6+8

		return outBinObj

	def _loadExpectedBinsCasesB_3d_to_1d(self):
		outEdges = [self.edgesB]
		outBinObj = tCode.NDimensionalBinnedResults(outEdges)

		outBinObj.initialiseCountsMatrix(countKey="fake_counts")
		outBinObj.binVals["fake_counts"][0] = 1+2+5+6
		outBinObj.binVals["fake_counts"][1] = 3+4+7+8

		return outBinObj	

	def testExpectedCaseA_3d_to_2d(self):
		expBinObj = self._loadExpectedBinsCaseA_3d_to_2d()
		actBinObj = self._runTestFunct()
		self.assertEqual(expBinObj, actBinObj)

	def testExpectedCaseB_3d_to_1d(self):
		self.keepDims = [1]
		expBinObj = self._loadExpectedBinsCasesB_3d_to_1d()
		actBinObj = self._runTestFunct()
		self.assertEqual(expBinObj, actBinObj)


#In reality I'll likely only actually do this on 1-dim objects but....
class TestAddRdfToBinObj(unittest.TestCase):

	def setUp(self):
		self.edgesA = [0,1,2,4]
		self.edgesB = [0,1,2]

		self.edgesTot = [self.edgesA, self.edgesB]

		self.volumeA = 40
		self.volumeB = 50

		self.numbFromA = 1
		self.numbFromB = 4
		self.numbToA = 20
		self.numbToB = 30
		self.volumes = [self.volumeA, self.volumeB]

		self.normCounts = [  [4,8],
		                     [3,6],
		                     [5,3] ]

		self.numbAtomsTo = [self.numbToA, self.numbToB]
		self.numbAtomsFrom = [self.numbFromA, self.numbFromB]

		self.createTestObjs()

	def createTestObjs(self):
		self.binObjA = tCode.NDimensionalBinnedResults(self.edgesTot, binVals={"normalised_counts":np.array(self.normCounts)})


	def _runTestFunct(self):
		tCode.addRdfValsToNDimBins(self.binObjA, self.numbAtomsFrom, self.numbAtomsTo, volumes=self.volumes)

	def _loadExpObj_2d(self):
		outObj = tCode.NDimensionalBinnedResults([self.edgesA,self.edgesB], binVals={"normalised_counts":self.normCounts})
		outObj.binVals["rdf"] = np.array( [ [48.6341681483221, 7.6553783196433],
		                                    [4.05284734569351, 0.637948193303608],
		                                    [0.450316371743724, 0.07088313258929] ] )
		outObj.binVals["rdf"] *= (1/(self.numbFromA*self.numbFromB))
		return outObj

	def testExpectedCase_1d(self):
		#Create expected
		expObj = tCode.NDimensionalBinnedResults([self.edgesA], binVals={"normalised_counts":np.array([ 12, 9,8 ]) })
		expObj.binVals["rdf"] = [ 7.63943726841098, 0.636619772367581, 0.070735530263065 ]

		#Create actual + run
		self.normCounts = [ 12, 9,8 ]
		self.volumes = self.volumeA
		self.edgesTot = [self.edgesA]
		self.numbAtomsTo = [self.numbToA]
		self.numbAtomsFrom = [self.numbFromA]
		self.createTestObjs()
		self._runTestFunct()
		actObj = self.binObjA

		self.assertEqual(expObj,actObj)

	def testExpectedCase_2d(self):
		expObj = self._loadExpObj_2d()
		self._runTestFunct()
		actObj = self.binObjA
		self.assertEqual(expObj, actObj)

	def testExpectedCase_2d_defaultVolumes(self):
		self.volumes = None
		expObj = tCode.NDimensionalBinnedResults([self.edgesA,self.edgesB], binVals={"normalised_counts":self.normCounts})
		expObj.binVals["rdf"] = np.array( [ [38.88, 6.12],
		                                       [3.24, 0.51],
		                                       [0.36, 0.056666666666667] ] )
		expObj.binVals["rdf"] *= (1/(self.numbFromA*self.numbFromB))
		self._runTestFunct()
		actObj = self.binObjA
		self.assertEqual(expObj, actObj)

	def testExpectedCase_2d_singleVol(self):
		self.volumes = 25
		expObj = tCode.NDimensionalBinnedResults([self.edgesA,self.edgesB], binVals={"normalised_counts":self.normCounts})
		expObj.binVals["rdf"] = np.array( [ [15.1981775463507, 2.39230572488853],
		                                    [1.26651479552922, 0.199358810407378],
		                                    [0.140723866169914, 0.022150978934153] ] )
		expObj.binVals["rdf"] *= (1/(self.numbFromA*self.numbFromB))
		self._runTestFunct()
		actObj = self.binObjA
		self.assertEqual(expObj, actObj)

class TestAddCircularRdfToBinObj(unittest.TestCase):

	def setUp(self):
		self.edgesA = [0,1,2,4]
		self.edgesB = [0,1,2]

		self.edgesTot = [self.edgesA, self.edgesB]

		self.areaA = 40
		self.areaB = 50

		self.numbFromA = 1
		self.numbFromB = 4
		self.numbToA = 20
		self.numbToB = 30
		self.areas = [self.areaA, self.areaB]

		self.normCounts = [  [4,8],
		                     [3,6],
		                     [5,3] ]

		self.numbAtomsTo = [self.numbToA, self.numbToB]
		self.numbAtomsFrom = [self.numbFromA, self.numbFromB]

		self.createTestObjs()

	def createTestObjs(self):
		self.binObjA = tCode.NDimensionalBinnedResults(self.edgesTot, binVals={"normalised_counts":np.array(self.normCounts)})

	def _runTestFunct(self):
		tCode.addCircularRdfToNDimBins(self.binObjA, self.numbAtomsFrom, self.numbAtomsTo, areas=self.areas)

	def _loadExpObj_2d(self):
		outObj = tCode.NDimensionalBinnedResults([self.edgesA,self.edgesB], binVals={"normalised_counts":self.normCounts})
		outObj.binVals["circular_rdf"] = np.array( [ [48.6341681483221, 22.9661349589299],
		                                             [12.1585420370805, 5.74153373973247],
		                                             [2.70189823046234, 1.27589638660722] ] )
		outObj.binVals["circular_rdf"] *= (1/(self.numbFromA*self.numbFromB))
		return outObj

	def testExpectedCase_2d(self):
		expObj = self._loadExpObj_2d()
		self._runTestFunct()
		actObj = self.binObjA
		self.assertEqual(expObj, actObj)



class TestAddProbabilitiesToBinObj(unittest.TestCase):

	def setUp(self):
		self.edgesA = [ 0, 1, 3, 5 ] #Non-uniform; to make things a bit annoying
		self.edgesB = [ 2, 4, 6 ]

		#Sort out the counts
		self.counts = np.zeros((3,2))
		self.counts[0][0] = 4
		self.counts[0][1] = 7
		self.counts[1][0] = 2
		self.counts[1][1] = 5
		self.counts[2][0] = 3
		self.counts[2][1] = 0

		#Args for the function
		self.countKey = "fake_counts_key"
		self.outKey = "prob_distrib_key"

		self.createTestObjs()

	def createTestObjs(self):
		binVals = {self.countKey: self.counts}
		edges = [self.edgesA, self.edgesB]
		self.testObjA = tCode.NDimensionalBinnedResults(edges, binVals=binVals) 

	def _runTestFunct(self):
		args = [self.testObjA]
		kwargs = {"countKey":self.countKey, "outKey":self.outKey}
		return tCode.addProbabilityDensitiesToNDimBinsSimple(*args, **kwargs)

	def _loadExpectedResultsA(self):
		expBinObj = copy.deepcopy(self.testObjA)
		expBinObj.initialiseCountsMatrix(countKey=self.outKey)

		outMatrix = expBinObj.binVals[self.outKey]
		totalArea = 5*4
		totalCounts = 4 + 7 + 2 +5 + 3 + 0
		outMatrix[0][0] = (4/totalCounts) * (1/2)  # * (2/totalArea)
		outMatrix[0][1] = (7/totalCounts) * (1/2) #* (2/totalArea)
		outMatrix[1][0] = (2/totalCounts) * (1/4) #* (4/totalArea)
		outMatrix[1][1] = (5/totalCounts) * (1/4) #* (4/totalArea)
		outMatrix[2][0] = (3/totalCounts) * (1/4) #* (4/totalArea)
		outMatrix[2][1] = 0

		return expBinObj

	def testExpectedResultsA(self):
		expBinObjA = self._loadExpectedResultsA()
		self._runTestFunct()
		actBinObjA = self.testObjA
		self.assertEqual(expBinObjA, actBinObjA)


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

		#Edges are [0,2,3,5,8]
		self.binEdges =   [0, 2, 3, 5, 8]
		self.pdfVals = [0.4, 0.2, 0.3, 0.6] #Dont care how physical these are


	def _runTestFunct(self):
		args = [self.binEdges, self.pdfVals]
		currKwargs = {"betweenVals":self.betweenVals,"normaliseBySum":self.normaliseSum}
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
		expVal = 0.8
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal,actVal)

	def testExpected_betweenVals_normByTotal(self):
		self.normByFullSum = True
		self.betweenVals = [0.5,5.5]
		expVal = 0.8/3.4
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





