
import copy
import os
import unittest

import numpy as np

import gen_basis_helpers.analyse_md.binned_res as tCode


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


class TestGetLowerDimObj_weightedAverageMethod(unittest.TestCase):

	def setUp(self):
		#define bins
		self.edgesA = [10,11,12] # e.g. planar distances
		self.edgesB = [0,1,2,3]  # e.g. number of hydrogen bonds

		#Run options
		self.keepDims = [0]

		self.createTestObjs()

	def createTestObjs(self):
		self.comboEdges = [self.edgesA, self.edgesB]
		self.binObjA = tCode.NDimensionalBinnedResults(self.comboEdges)

		#Put various values there
		#Use a weird key to make sure we dont rely on it being the default value
		self.binObjA.initialiseCountsMatrix(countKey="fake_counts")
		self.binObjA.binVals["fake_counts"][0][0] = 2
		self.binObjA.binVals["fake_counts"][0][1] = 5
		self.binObjA.binVals["fake_counts"][0][2] = 1
		self.binObjA.binVals["fake_counts"][1][0] = 8
		self.binObjA.binVals["fake_counts"][1][1] = 1

	def _runTestFunct(self):
		return tCode.getLowerDimBinObj_weightedAverageMethod(self.binObjA, self.keepDims, useCountKey="fake_counts")

	def testRaisesForMoreThanTwoDims(self):
		""" Not sure what to do for >2 dims at the moment (likely it'd be pointless anyway though) """
		self.binObjA = tCode.NDimensionalBinnedResults( self.comboEdges + [[6,7]] )
		with self.assertRaises(NotImplementedError):
			self._runTestFunct()

	def testExpectedValsA(self):
		#Create the expected output bin object
		expBinObj = tCode.NDimensionalBinnedResults([self.edgesA])
		expBinObj.initialiseCountsMatrix(countKey="fake_counts")
		expBinObj.binVals["fake_counts"][0] = self.binObjA.binVals["fake_counts"][0][0] + self.binObjA.binVals["fake_counts"][0][1] + self.binObjA.binVals["fake_counts"][0][2]
		expBinObj.binVals["fake_counts"][1] = self.binObjA.binVals["fake_counts"][1][0] + self.binObjA.binVals["fake_counts"][1][1]
		expBinObj.binVals["weighted_average"] = np.array( [sum([ (2/8)*0.5, (5/8)*1.5, (1/8)*2.5 ]),
		                                                  sum([ (8/9)*0.5, (1/9)*1.5 ])] )   

		actBinObj = self._runTestFunct()
		self.assertEqual(expBinObj,actBinObj)


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

		#Using Vol(bin) = Vol(upper)-Vol(lower) this should lead to bin volumes of:
		#binVolsA = [4.18879020478639, 29.3215314335047, 234.572251468038]
		#binVolsB = [4.18879020478639, 29.3215314335047]




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

#		outObj.binVals["rdf"] = np.array( [ [113.98633159763, 32.5675233136086],
#		                                    [12.2128212426032, 3.48937749788663],
#		                                    [2.54433775887567, 0.218086093617915] ] )

		outObj.binVals["rdf"] *= (1/(self.numbFromA*self.numbFromB))
		return outObj

	def testExpectedCase_1d(self):
		#Create expected
		expObj = tCode.NDimensionalBinnedResults([self.edgesA], binVals={"normalised_counts":np.array([ 12, 9,8 ]) })
		expObj.binVals["rdf"] = [ 5.72957795130823, 0.613883351925882, 0.068209261325098 ]

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

	@unittest.skip("Couldnt figure out expected answer again...")
	def testExpectedCase_2d(self):
		expObj = self._loadExpObj_2d()
		self._runTestFunct()
		actObj = self.binObjA
		self.assertEqual(expObj, actObj)

	@unittest.skip("Couldnt figure out expected answer again...")
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

	@unittest.skip("Couldnt figure out expected answer again...")
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


class TestGetNDimbinsSummedSimple(unittest.TestCase):

	def setUp(self):
		self.edgesA = [0,1,2]
		self.edgesB = [3,4,5,6]
		self.edgesC = [7,8,9,10] #Will be useful for checking errors raises properly

		#Run options
		self.countKey = "counts"
		self.normCountKey = "normalised_counts"

		#Setup counts arrays and add a few
		self.countsA = np.zeros( (2,3) )
		self.countsB = np.zeros( (2,3) )
	
		self.countsA[0][1] = 4 
		self.countsA[1][2] = 6
		self.countsA[1][1] = 2

		self.countsB[0][0] = 12
		self.countsB[1][1] = 14
	
		#Define nSteps for each bin so we can get normalised counts in a sensible way
		self.nStepsA, self.nStepsB = 1,2

		self.createTestObjs()

	def createTestObjs(self):
		self.binValsA = {self.countKey:self.countsA, self.normCountKey: self.countsA/self.nStepsA}
		self.binValsB = {self.countKey:self.countsB, self.normCountKey: self.countsB/self.nStepsB}

		self.binObjA = tCode.NDimensionalBinnedResults([self.edgesA,self.edgesB], binVals=self.binValsA)
		self.binObjB = tCode.NDimensionalBinnedResults([self.edgesA,self.edgesB], binVals=self.binValsB)

		self.binObjsToAdd = [self.binObjA, self.binObjB]

	def _runTestFunct(self):
		currKwargs = {"countKey":self.countKey,"normCountKey":self.normCountKey}
		return tCode.getSummedNDimensionalBinnedResSimple(self.binObjsToAdd, **currKwargs)

	def _loadExpCountsA(self):
		return np.array(  [ [12,4,0], [0,16,6] ] )

	def _loadExpNormCountsA(self):
		totalSteps = self.nStepsA + self.nStepsB
		expNormCounts = np.array( [ [12/totalSteps, 4/totalSteps, 0],
		                             [ 0, (14+2)/totalSteps, 6/totalSteps ] ] )
		return expNormCounts

	def testExpectedCaseA(self):
		#Manually get the exp counts and norm-counts matrices
		expCounts = self._loadExpCountsA()
		expNormCounts = self._loadExpNormCountsA()
		expBinVals = {self.countKey: expCounts, self.normCountKey: expNormCounts}

		#Get the expected bin obj
		expBinObj = tCode.NDimensionalBinnedResults( [self.edgesA,self.edgesB], binVals=expBinVals)
		actBinObj = self._runTestFunct()

		self.assertEqual(expBinObj,actBinObj)

	def testExpected_normCountsMissing(self):
		self.binObjA.binVals.pop(self.normCountKey)
		self.binObjB.binVals.pop(self.normCountKey)

		#Load expected bin values
		expCounts = self._loadExpCountsA()
		expBinVals = {self.countKey:expCounts}

		#
		expBinObj = tCode.NDimensionalBinnedResults( [self.edgesA,self.edgesB],binVals=expBinVals )
		actBinObj = self._runTestFunct()

		self.assertEqual(expBinObj,actBinObj)

	def testExpectedCountsMissing(self):
		self.binObjA.binVals.pop(self.countKey)
		self.binObjB.binVals.pop(self.countKey)

		with self.assertRaises(KeyError):
			self._runTestFunct()

	def testRaisesWhenEdgesDontMatch(self):
		self.binObjB = tCode.NDimensionalBinnedResults([self.edgesA,self.edgesC], binVals=self.binValsB)
		self.binObjsToAdd = [self.binObjA, self.binObjB]

		with self.assertRaises(ValueError):
			self._runTestFunct()

	def testExpected_diffCountKey(self):
		self.countKey = "new_count_key"
		self.createTestObjs()

		#Manually get the exp counts and norm-counts matrices
		expCounts = self._loadExpCountsA()
		expNormCounts = self._loadExpNormCountsA()
		expBinVals = {self.countKey: expCounts, self.normCountKey: expNormCounts}

		#Get the expected bin obj
		expBinObj = tCode.NDimensionalBinnedResults( [self.edgesA,self.edgesB], binVals=expBinVals)
		actBinObj = self._runTestFunct()

		self.assertTrue( "new_count_key" in expBinObj.binVals.keys() )
		self.assertEqual(expBinObj,actBinObj)

	def testExpected_threeBinsAdded(self):
		self.binObjsToAdd = [self.binObjA, self.binObjB, self.binObjA]

		#Generate expected
		expCounts = np.array(  [ [12,8,0],[0,18,12] ] )
		totalSteps = (2*self.nStepsA) + self.nStepsB
		expNormCounts = np.array( [ [12/totalSteps, 8/totalSteps,0],
		                            [0, 18/totalSteps, 12/totalSteps] ] )
		expBinVals = {self.countKey: expCounts, self.normCountKey: expNormCounts}

		#Get the expected bin object
		expBinObj = tCode.NDimensionalBinnedResults( [self.edgesA,self.edgesB], binVals=expBinVals)
		actBinObj = self._runTestFunct()

		self.assertEqual(expBinObj,actBinObj)


