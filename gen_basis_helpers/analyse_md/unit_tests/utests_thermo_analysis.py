
import copy
import itertools as it
import math
import unittest
import unittest.mock as mock


import gen_basis_helpers.analyse_md.thermo_data as thermoDataHelp
import gen_basis_helpers.analyse_md.analyse_thermo as tCode



class TestGetSimpleCentralWindowAverageFromThermoData(unittest.TestCase):

	def setUp(self):
		self.temp = [10,15,20]
		self.step = [10,20,30]
		self.time = [20,40,60]
		self.prop = "temp"
		self.startIdx = 0
		self.widthEachSide = 3
		self.fillVal = None
		self.createTestObjs()

	def createTestObjs(self):
		thermoDataDict = {"step":self.step, "time":self.time, "temp":self.temp}
		self.thermoDataObjA = thermoDataHelp.ThermoDataStandard(thermoDataDict)

	def _runTestFunct(self):
		args = [self.thermoDataObjA, self.prop, self.widthEachSide]
		kwargs = {"fillVal":self.fillVal, "startIdx":self.startIdx} 
		return tCode.getSimpleCentralWindowAverageFromThermoData(*args,**kwargs)

	@mock.patch("gen_basis_helpers.analyse_md.analyse_thermo._getSimpleCentralWindowAverageFromIter")
	def testExpectedCallsA(self, mockGetMovingAvg):
		expVal = mock.Mock()
		mockGetMovingAvg.side_effect = lambda *args,**kwargs: expVal
		actVal = self._runTestFunct()
		mockGetMovingAvg.assert_called_with(self.temp, self.widthEachSide, fillValWhenNotCalc=self.fillVal)
		self.assertEqual(expVal, actVal)

	@mock.patch("gen_basis_helpers.analyse_md.analyse_thermo._getSimpleCentralWindowAverageFromIter")
	def testExpectedCalls_startIdxOne(self, mockGetMovingAvg):
		expVal = mock.Mock()
		mockGetMovingAvg.side_effect = lambda *args,**kwargs: expVal
		self.startIdx = 1
		actVal = self._runTestFunct()
		mockGetMovingAvg.assert_called_with(self.temp[self.startIdx:], self.widthEachSide, fillValWhenNotCalc=self.fillVal)
		self.assertEqual(expVal, actVal)


class TestGetSimpleCentralWindowAverageFromIter(unittest.TestCase):

	def setUp(self):
		self.width = 2
		self.fillValWhenNotCalc = None
		self.iterA = [1,2,3,4,5,6,7,8,10]

	def _runTestFunct(self):
		args = [self.iterA, self.width]
		kwargs = {"fillValWhenNotCalc":self.fillValWhenNotCalc}
		return tCode._getSimpleCentralWindowAverageFromIter(*args, **kwargs)

	def testSimpleCaseA(self):
		expVals = [None, None, 15/5, 20/5, 25/5, 30/5, 36/5, None, None]
		actVals = self._runTestFunct()
		self.assertAlmostEqual(expVals, actVals)

	def testExpectedForSuperShortList(self):
		self.iterA = [1,2]
		expVals = [None,None]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

class TestGetMovingAverageFromThermoData(unittest.TestCase):


	def setUp(self):
		self.temp = [10,15,20]
		self.step = [10,20,30]
		self.time = [20,40,60]
		self.prop = "temp"
		self.startIdx = 0
		self.createTestObjs()

	def createTestObjs(self):
		thermoDataDict = {"step":self.step, "time":self.time, "temp":self.temp}
		self.thermoDataObjA = thermoDataHelp.ThermoDataStandard(thermoDataDict)

	def _runTestFunct(self):
		return tCode.getMovingAverageFromThermoData(self.thermoDataObjA, self.prop, startIdx=self.startIdx)

	def testExpectedValsSimpleCaseA(self):
		expVals = [10, 25/2, 45/3]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expVals,actVals)]

	def testExpectedForStartIdxOne(self):
		self.startIdx=1
		self.createTestObjs()
		expVals = [15, 35/2]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expVals,actVals)]

class TestGetStandardStatsDictForThermoProps(unittest.TestCase):

	def setUp(self):
		self.temp = [10,15,20]
		self.step = [10,20,30]
		self.time = [20,40,60]
		self.startTime = None
		self.endTime = None
		self.startStep = None #zero based numbering i guess
		self.endStep = None
		self.timeTol = 1e-3
		self.props = ["temp"]
		self.createTestObjs()

	def createTestObjs(self):
		thermoDataDict = {"step":self.step, "time":self.time, "temp":self.temp}
		self.thermoDataObjA = thermoDataHelp.ThermoDataStandard(thermoDataDict)
		kwargDict = {"startTime":self.startTime, "endTime":self.endTime,
		             "startStep":self.startStep, "endStep":self.endStep,
		             "props":self.props, "timeTol":self.timeTol}
		self.testObjA = tCode.GetStatsForThermoProps(**kwargDict)

	def testExpectedForSimpleDataA(self):
		actDict = self.testObjA.create(data=self.thermoDataObjA)
		expDict = self._loadExpectedDictSimpleDataA()
		self.assertEqual(expDict, actDict)

	def testExpectedForSimpleData_startTime40(self):
		self.startTime = 40 - 1e-5 #Really 40; just checking we handle this for float errors
		self.createTestObjs()
		expStdDev = math.sqrt( (2*(2.5**2)) /2)
		kwargDict = {"mean":17.5, "standardDev":expStdDev, "minVal":15,
		             "maxVal":20, "nVals":2}
		expDict = {"temp":tCode.StandardStatsObject(**kwargDict)}
		actDict = self.testObjA.create(data=self.thermoDataObjA)
		self.assertEqual(expDict, actDict)

	def testExpectedForSimpleData_endTime40(self):
		self.endTime = 40 - 1e-5 ##Really 40; just checking we handle this for float errors
		self.createTestObjs()
		expStdDev = math.sqrt( (2*(2.5**2)) /2)
		kwargDict = {"mean":12.5, "standardDev":expStdDev, "minVal":10,
		             "maxVal":15, "nVals":2}
		expDict = {"temp":tCode.StandardStatsObject(**kwargDict)}
		actDict = self.testObjA.create(data=self.thermoDataObjA)
		self.assertEqual(expDict,actDict)

	def testStartStepOveridesStartTime(self):
		self.startTime = 40
		self.startStep = 0
		self.createTestObjs()
		expDict = self._loadExpectedDictSimpleDataA()
		actDict = self.testObjA.create(data=self.thermoDataObjA)
		self.assertEqual(expDict,actDict)

	def testEndStepOveridesEndTime(self):
		self.endTime = 40
		self.endStep = 100
		self.createTestObjs()
		expDict = self._loadExpectedDictSimpleDataA()
		actDict = self.testObjA.create(data=self.thermoDataObjA)
		self.assertEqual(expDict, actDict)

	def _loadExpectedDictSimpleDataA(self):
		expStdDev = math.sqrt(50/3)
		kwargDict = {"mean":15, "standardDev":expStdDev,"minVal":10,
		             "maxVal":20, "nVals":3}
		expDict = {"temp":tCode.StandardStatsObject(**kwargDict)}
		return expDict

class TestStandardStatsObject(unittest.TestCase):

	def setUp(self):
		self.mean = 4
		self.stdDev = 20
		self.nVals = 9
		self.minVal = 2
		self.maxVal = 10
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"mean":self.mean, "standardDev":self.stdDev,
		             "nVals":self.nVals, "minVal":self.minVal, "maxVal":self.maxVal}
		self.testObjA = tCode.StandardStatsObject(**kwargDict)

	def testRangeGivesExpectedVal(self):
		expRange = abs(self.maxVal-self.minVal)
		actRange = self.testObjA.range
		self.assertEqual(expRange,actRange)

	def testRangeReturnsNoneIfOneOfValsNotSet(self):
		self.maxVal = None
		self.createTestObjs()
		expRange = None
		actRange = self.testObjA.range
		self.assertEqual(expRange, actRange)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA,objB)

	def testEqualObjsCompareEqual_nValsBothNone(self):
		self.nVals = None
		self.createTestObjs()
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA,objB)

	def testUnequalObjsCompareUnequal_diffMean(self):
		objA = copy.deepcopy(self.testObjA)
		self.mean += 1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)
		self.assertNotEqual(objB,objA)

	def testUnequalObjsCompareUnequal_nStepsNoneOnOne(self):
		objA = copy.deepcopy(self.testObjA)
		self.nVals = None
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)

	

