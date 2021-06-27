

import copy
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

