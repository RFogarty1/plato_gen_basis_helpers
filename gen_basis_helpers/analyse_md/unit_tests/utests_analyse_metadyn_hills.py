
import copy
import itertools as it
import unittest
import unittest.mock as mock

import numpy as np

import gen_basis_helpers.analyse_md.analyse_metadyn_hills as tCode


class TestGetPotValForEachStepOverTime(unittest.TestCase):

	def setUp(self):
		#Create the info object
		self.infoObj = mock.Mock()

		#Options for the run function
		self.timeRange = None
		self.minTimeDiff = 1e-3
		self.timeTol = 1e-3
		self.oneDimOutput = True
		self.inpValsA = [ [1,2], [3,4] ]

		#Stuff for mocking
		self.times = [1,2]
		self.expOutValsA = [3,4]
		self.expOutValsB = [5,6]

		self.createTestObjs()

	def createTestObjs(self):
		self.expDataA = [ [inpVal,outVal] for inpVal,outVal in it.zip_longest(self.inpValsA, self.expOutValsA) ]
		self.expDataB = [ [inpVal,outVal] for inpVal,outVal in it.zip_longest(self.inpValsA, self.expOutValsB) ]

		#
		self.expPotObjs = [mock.Mock(), mock.Mock()]
		self.expPotObjs[0].evalFunctAtVals.side_effect = lambda *args,**kwargs: self.expOutValsA
		self.expPotObjs[1].evalFunctAtVals.side_effect = lambda *args,**kwargs: self.expOutValsB

	def _runTestFunct(self):
		args = [self.infoObj, self.inpValsA]
		currKwargs = {"timeRange":self.timeRange, "timeTol":self.timeTol,
		              "minTimeDiff":self.minTimeDiff, "oneDimOutput":self.oneDimOutput}
		return tCode.getTimeVsPotAtValsForEachHillAdded(*args, **currKwargs)

	@mock.patch("gen_basis_helpers.analyse_md.analyse_metadyn_hills._getTimeVsPotObjsForEachTimeHillAdded")
	def testExpectedValsA_twoDimCase(self, mockGetTimeVsPotObjects):
		#Setup test objs
		mockGetTimeVsPotObjects.side_effect = lambda *args,**kwargs: [ [t,potObj] for t,potObj in it.zip_longest(self.times, self.expPotObjs) ]

		expOutput = [ [self.times[0], self.expDataA],
		              [self.times[1], self.expDataB] ]

		actOutput = self._runTestFunct()

		mockGetTimeVsPotObjects.assert_called_with(self.infoObj, timeRange=self.timeRange, timeTol=self.timeTol, minTimeDiff=self.minTimeDiff)
		self.expPotObjs[0].evalFunctAtVals.assert_called_with(self.inpValsA)
		self.expPotObjs[1].evalFunctAtVals.assert_called_with(self.inpValsA)
		self.assertEqual(expOutput, actOutput)

	@mock.patch("gen_basis_helpers.analyse_md.analyse_metadyn_hills._getTimeVsPotObjsForEachTimeHillAdded")
	def testExpectedValsA_oneDimCase(self, mockGetTimeVsPotObjects):
		#Swap to 1-dimensional
		self.inpValsA = [ [2], [3] ]
		self.createTestObjs()

		#Setup test objs
		mockGetTimeVsPotObjects.side_effect = lambda *args,**kwargs: [ [t,potObj] for t,potObj in it.zip_longest(self.times, self.expPotObjs) ]

		expDataA = [ [x[0],y] for x,y in self.expDataA ]
		expDataB = [ [x[0],y] for x,y in self.expDataB ]
		expOutput = [ [self.times[0], expDataA], [self.times[1], expDataB] ]

		actOutput = self._runTestFunct()

		mockGetTimeVsPotObjects.assert_called_with(self.infoObj, timeRange=self.timeRange, timeTol=self.timeTol, minTimeDiff=self.minTimeDiff)
		self.assertEqual(expOutput, actOutput)

class TestGetPotObjForEachStepOverTime(unittest.TestCase):

	def setUp(self):
		#Create the info object
		self.times = [2,3,1]
		self.positions = [ [3,4], [5,6], [7,8] ]
		self.scales =  [ [4,5], [6,7], [8,9] ]
		self.heights = [ [2,3], [4,5], [6,7] ]
		self.sortTimes = False

		#Options for the run function
		self.timeRange = None
		self.minTimeDiff = 1e-3
		self.timeTol = 1e-3
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"times":self.times, "positions":self.positions, "scales":self.scales,
		              "heights":self.heights, "sortTimes":self.sortTimes}
		self.testObjA = tCode.MetadynHillsInfo(**currKwargs)

	def _runTestFunct(self):
		currKwargs = {"timeRange":self.timeRange, "timeTol":self.timeTol, "minTimeDiff":self.minTimeDiff}
		return tCode._getTimeVsPotObjsForEachTimeHillAdded(self.testObjA, **currKwargs)

	def _loadMultiDimHillsInTimeOrder(self):
		hillA = tCode.MultiDimGaussHill.fromIters(heights=[6,7], scales=[8,9], positions=[7,8])
		hillB = tCode.MultiDimGaussHill.fromIters(heights=[2,3], scales=[4,5], positions=[3,4])
		hillC = tCode.MultiDimGaussHill.fromIters(heights=[4,5], scales=[6,7], positions=[5,6])
		return [hillA, hillB, hillC]

	def testExpectedFunctValA_fullRange(self):
		#Figure out the expected potential objects
		hillA, hillB, hillC = self._loadMultiDimHillsInTimeOrder()

		groupA = tCode.GroupedMultiDimGaussHills([hillA])
		groupB = tCode.GroupedMultiDimGaussHills([hillA,hillB])
		groupC = tCode.GroupedMultiDimGaussHills([hillA,hillB,hillC])

		expTimesVsHills = [ [1,groupA], [2,groupB], [3,groupC] ]
		actTimesVsHills = self._runTestFunct()
		self.assertEqual(expTimesVsHills, actTimesVsHills)

	def testExpected_minTimeSet(self):
		hillA, hillB, hillC = self._loadMultiDimHillsInTimeOrder()
		self.timeRange = [1.4, 20]
		groupA = tCode.GroupedMultiDimGaussHills([hillB])
		groupB = tCode.GroupedMultiDimGaussHills([hillB, hillC])

		expTimesVsHills = [ [2,groupA], [3,groupB] ]
		actTimesVsHills = self._runTestFunct()
		self.assertEqual(expTimesVsHills, actTimesVsHills)

	def testExpected_maxTimeSet(self):
		hillA, hillB, hillC = self._loadMultiDimHillsInTimeOrder()
		self.timeRange = [0,2.2]
		groupA = tCode.GroupedMultiDimGaussHills([hillA])
		groupB = tCode.GroupedMultiDimGaussHills([hillA,hillB])

		expTimesVsHills = [ [1,groupA], [2,groupB] ]
		actTimesVsHills = self._runTestFunct()
		self.assertEqual(expTimesVsHills, actTimesVsHills)

	def testExpected_ridicLargeMinTimeDiff(self):
		self.minTimeDiff = 10
		hillA, hillB, hillC = self._loadMultiDimHillsInTimeOrder()
		groupA = tCode.GroupedMultiDimGaussHills([hillA,hillB,hillC])

		expTimesVsHills = [ [1,groupA] ]
		actTimesVsHills = self._runTestFunct()
		self.assertEqual(expTimesVsHills, actTimesVsHills)


class TestEvalPotOverTimeRangeHillsInfoObj(unittest.TestCase):

	def setUp(self):
		#Create the info object
		self.times = [2,3,1]
		self.positions = [ [3,4], [5,6], [7,8] ]
		self.scales =  [ [4,5], [6,7], [8,9] ]
		self.heights = [ [2,3], [4,5], [6,7] ]
		self.sortTimes = False

		#options for running the function
		self.timeRange = None
		self.timeTol = 1e-3
		self.evalAtVals = [ [1,1], [3,4] ]
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"times":self.times, "positions":self.positions, "scales":self.scales,
		              "heights":self.heights, "sortTimes":self.sortTimes}
		self.testObjA = tCode.MetadynHillsInfo(**currKwargs)

	def _runTestFunct(self):
		args = self.testObjA, self.evalAtVals
		kwargs = {"timeRange":self.timeRange, "timeTol":self.timeTol}
		return tCode.evalPotAddedOverTimeRangeForHillsInfoObj(*args, **kwargs)

	#Results generated from excel
	def testExpectedFullRangeA(self):
		actVals = self._runTestFunct()
		expVals = [7.30237550948453, 68.1534473578693]
		expVals = [40.2600020863019, 57.7416360491726]
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expVals, actVals)]

	def testExpectedSmallerTimeRange(self):
		self.timeRange = [0.8,2.2]
		expVals = [27.8511758765677, 39.5791113014523]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expVals, actVals)]

class TestGetCombinedMetaDynHillsClass(unittest.TestCase):

	def setUp(self):
		#First obj
		self.timesA = [ 1 ]
		self.positionsA = [ [3,4] ]
		self.scalesA = [ [4,5] ]
		self.heightsA = [ [6,7] ]
		self.sortTimesA = False

		#Second obj
		self.timesB = [ 2 ]
		self.positionsB = [ [6,4] ]
		self.scalesB = [ [7,3] ]
		self.heightsB = [ [8,2] ]
		self.sortTimesB = False

		#Anything needed for running test funct can go here
		self.copyData = False
		self.createTestObjs()

	def createTestObjs(self):
		kwargsA = {"times":self.timesA, "positions":self.positionsA, "scales":self.scalesA,
		           "heights":self.heightsA, "sortTimes":self.sortTimesA}
		kwargsB = {"times":self.timesB, "positions":self.positionsB, "scales":self.scalesB,
		           "heights":self.heightsB, "sortTimes":self.sortTimesB}

		self.testObjA = tCode.MetadynHillsInfo(**kwargsA)
		self.testObjB = tCode.MetadynHillsInfo(**kwargsB)

	def _runTestFunct(self):
		return tCode.getMergedMetadynHillsInfoInstance([self.testObjA, self.testObjB], copyData=self.copyData)

	def _loadExpObjCombinedAB(self):
		kwargs = {"times": self.timesA+self.timesB, "positions":self.positionsA+self.positionsB,
		          "scales": self.scalesA+self.scalesB, "heights":self.heightsA+self.heightsB}
		return tCode.MetadynHillsInfo(**kwargs)

	def testExpectedFromMergingTwo(self):
		expObj = self._loadExpObjCombinedAB()
		actObj = self._runTestFunct()
		self.assertEqual(expObj, actObj)


	def testExpectedFromMergingTwo_withCopy(self):
		#
		self.copyData = True
		expObj = copy.deepcopy( self._loadExpObjCombinedAB() )
		actObj = self._runTestFunct()
		self.assertEqual(expObj, actObj)

		#Change one of the input lists; if datas been copied it shouldnt affect the output
		self.testObjA.times[0] += 1
		self.assertEqual(expObj, actObj)
		

class TestMetadynHillsInfoClass(unittest.TestCase):

	def setUp(self):
		self.times = [2,1,3]
		self.positions = [ [3,4], [5,6], [7,8] ]
		self.scales =  [ [4,5], [6,7], [8,9] ]
		self.heights = [ [2,3], [4,5], [6,7] ]
		self.sortTimes = False
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"times":self.times, "positions":self.positions, "scales":self.scales,
		              "heights":self.heights, "sortTimes":self.sortTimes}
		self.testObjA = tCode.MetadynHillsInfo(**currKwargs)

	def testGetTimesWithinRange_noRangeSpecified(self):
		expTimes = self.times
		actTimes = self.testObjA.getTimesWithinRange()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expTimes,actTimes)]

	def testGetTimesWithinRange_rangeSpecified(self):
		timeRange = [0.8,2.9]
		expTimes = [2,1]
		actTimes = self.testObjA.getTimesWithinRange(timeRange=timeRange)
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expTimes,actTimes)]

	@mock.patch("gen_basis_helpers.analyse_md.analyse_metadyn_hills._getIterOfMultiDimHills")
	def testCreateMultiDimHillsCallsExpected_noTimeRange(self, mockGetMultiDimHills):
		expOutput = mock.Mock()
		mockGetMultiDimHills.side_effect = lambda *args,**kwargs: expOutput
		actOutput = self.testObjA.createMultiDimHills()
		mockGetMultiDimHills.assert_called_with(self.heights, self.scales, self.positions)
		self.assertEqual(expOutput, actOutput)

	@mock.patch("gen_basis_helpers.analyse_md.analyse_metadyn_hills._getIterOfMultiDimHills")
	def testCreateMultiDimHillsCallsExpected_timeLimitsToSingleIdx(self, mockGetMultiDimHills):
		timeRange = [1.5,2.5]
		expOutput = mock.Mock()
		mockGetMultiDimHills.side_effect = lambda *args, **kwargs: expOutput
		actOutput = self.testObjA.createMultiDimHills(timeRange=timeRange)

		currArgs = [ [self.heights[0]], [self.scales[0]], [self.positions[0]] ]
		mockGetMultiDimHills.assert_called_with(*currArgs)
		self.assertEqual(expOutput, actOutput)

	@mock.patch("gen_basis_helpers.analyse_md.analyse_metadyn_hills.GroupedMultiDimGaussHills")
	@mock.patch("gen_basis_helpers.analyse_md.analyse_metadyn_hills.MetadynHillsInfo.createMultiDimHills")
	def testCreateGroupedHills(self, mockCreateMultiDimHills,  mockGroupedHillsCls):
		#Setup mocks
		expHillsObj, expGroupedHillsObj = mock.Mock(), mock.Mock()
		mockCreateMultiDimHills.side_effect = lambda *args,**kwargs: expHillsObj
		mockGroupedHillsCls.side_effect = lambda *args,**kwargs: expGroupedHillsObj

		#Test
		actOutput = self.testObjA.createGroupedHills(timeRange=None, timeTol=1)
		mockCreateMultiDimHills.assert_called_with(timeRange=None, timeTol=1)
		mockGroupedHillsCls.assert_called_with(expHillsObj)
		self.assertEqual(expGroupedHillsObj, actOutput)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA,objB)

	def testUnequalObjsCompareUnequal_diffTimes(self):
		objA = copy.deepcopy(self.testObjA)
		self.times[1] += 2
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffScales(self):
		objA = copy.deepcopy(self.testObjA)
		self.scales[0][1] += 2
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testSortTimesFunction(self):
		#Figure out expected obj
		expIdxOrder, kwargDict = [ 1,0,2 ], dict()
		for attr in ["times","positions","scales","heights"]:
			startVals = getattr(self,attr)
			kwargDict[attr] = [copy.deepcopy(startVals[idx]) for idx in expIdxOrder]
		expObj = tCode.MetadynHillsInfo(**kwargDict)

		self.assertNotEqual(expObj, self.testObjA)
		self.testObjA.sortTimes()
		self.assertEqual(expObj, self.testObjA)

	def testMultiplyHeightsByConstant(self):
		objA = copy.deepcopy(self.testObjA)
		multFactor = 2
		self.heights = [ [multFactor*x for x in [2,3]], [multFactor*x for x in [4,5]], [multFactor*x for x in [6,7]] ]
		self.createTestObjs()
		objB = self.testObjA
		objA.multiplyHeightsByFactor(multFactor)
		self.assertEqual(objA,objB)

class TestMetadynHillsInfoClass_indicesWithinTimeRanges(unittest.TestCase):

	def setUp(self):
		self.times = [2,1,4]
		self.timeRange = [0.5,5]
		self.timeTol = 1e-3
		self.sortTimes = False
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.MetadynHillsInfo(times=self.times, sortTimes=self.sortTimes)

	def _runTestFunct(self):
		args = [self.timeRange, self.timeTol]
		return self.testObjA._getIndicesWithinTimeRange(*args)

	def testExpectedCase_allWithinRange(self):
		expIndices = [0,1,2]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testExpectedCase_oneLessThanMinTime(self):
		self.timeRange = [1.5,5]
		expIndices = [0,2]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testExpectedCase_oneGreaterThanMaxTime(self):
		self.timeRange = [0,3]
		expIndices = [0,1]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testExpectedCase_oneJustOutOfTolerance(self):
		self.timeRange = [0, 4-(1.1*self.timeTol)]
		expIndices = [0,1]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testExpectedCase_oneJustWithinTolerance(self):
		self.timeRange = [0, 4-(0.9*self.timeTol)]
		expIndices = [0,1,2]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testExpected_reverseTimeRange(self):
		self.timeRange = sorted(self.timeRange, reverse=True)
		with self.assertRaises(ValueError):
			self._runTestFunct()

class TestGroupedMultiDimGaussHills(unittest.TestCase):

	def setUp(self):
		self.gauA = tCode.MultiDimGaussHill.fromIters(heights=[4], scales=[3], positions=[2])
		self.gauB = tCode.MultiDimGaussHill.fromIters(heights=[5], scales=[4], positions=[3])
		self.positions = [[-1],[0],[1]]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.GroupedMultiDimGaussHills([self.gauA, self.gauB])

	def testExpectedIndividualContribsArray(self):
		expVals = [ [2.42612263885053, 3.03265329856317],
		            [3.20294961166723, 3.77419800994504],
		            [3.78383787562706, 4.41248451292298] ]

		actVals = self.testObjA.getContribsAtPositions(self.positions)
		self.assertTrue( np.allclose( np.array(expVals),np.array(actVals) ) )

	def testExpectedTotalVals(self):
		expVals = [5.4587759374137, 6.97714762161227, 8.19632238855004]
		actVals = self.testObjA.evalFunctAtVals(self.positions)
		self.assertTrue( np.allclose( np.array(expVals),np.array(actVals) ) )

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffLengthMultiDimGaus(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		self.testObjA.multiDimGaus.append(self.gauA)
		objB = self.testObjA
		self.assertNotEqual(objA, objB)		

	def testUnequalObjsCompareUnequal_diffGauB(self):
		objA = copy.deepcopy(self.testObjA)
		self.gauA = tCode.MultiDimGaussHill.fromIters(heights=[4+2], scales=[3], positions=[2])
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)


class TestMultiDimGauHillFunction(unittest.TestCase):

	def setUp(self):
		self.scales = [3,4]
		self.gauPos = [2,3]
		self.heights = [4,5]
		self.positions = [[-1,0], [0,1] , [1,2] ]
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"height":self.heights[0], "scale":self.scales[0], "pos":self.gauPos[0]}
		self.gauA = tCode.OneDimGaussianHill(**currKwargs)
		currKwargs = {"height":self.heights[1], "scale":self.scales[1], "pos":self.gauPos[1]}
		self.gauB = tCode.OneDimGaussianHill(**currKwargs)
		self.multiGau = tCode.MultiDimGaussHill([self.gauA,self.gauB])

	def _runTestFunct(self):
		return self.multiGau.evalFunctAtVals(self.positions)

	def testExpValsA(self):
		expVals = [9.15666723543229, 14.1329655571543, 18.3371071146406]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expVals,actVals)]

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.multiGau)
		self.createTestObjs()
		objB = self.multiGau
		self.assertEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffNumbGauHills(self):
		objA = copy.deepcopy(self.multiGau)
		self.createTestObjs()
		self.multiGau.oneDimHills.append(self.gauA)
		objB = self.multiGau
		self.assertNotEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffPositions(self):
		objA = copy.deepcopy(self.multiGau)
		self.gauPos[0] += 2
		self.createTestObjs()
		objB = self.multiGau
		self.assertNotEqual(objA, objB)

	def testFromIters(self):
		testObj = tCode.MultiDimGaussHill.fromIters(heights=self.heights, positions=self.gauPos, scales=self.scales)
		self.assertEqual(self.multiGau, testObj)


#Want this to have properties related to the scale/pos etc. vars
class TestOneDimGauHillFunction(unittest.TestCase):

	def setUp(self):
		self.scale = 3 
		self.gauPos = 1
		self.height = 2
		self.positions = [-1,0,1]
		self.createTestObjs()

	def createTestObjs(self):
		kwargs = {"height":self.height, "scale":self.scale, "pos":self.gauPos}
		self.gauA = tCode.OneDimGaussianHill(**kwargs)

	def _runTestFunct(self):
		return self.gauA.evalFunctAtVals(self.positions)

	def testExpectedValA(self):
		expVals = [1.60147480583362, 1.89191893781353, 2] #Just calculated in excel
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expVals, actVals)]

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.gauA)
		self.createTestObjs()
		objB = self.gauA
		self.assertEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffPosition(self):
		objA = copy.deepcopy(self.gauA)
		self.gauPos += 1
		self.createTestObjs()
		objB = self.gauA
		self.assertNotEqual(objA, objB)



