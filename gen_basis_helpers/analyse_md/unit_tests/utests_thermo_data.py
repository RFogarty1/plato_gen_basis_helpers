
import copy
import unittest
import unittest.mock as mock

import numpy as np

import gen_basis_helpers.analyse_md.thermo_data as tCode

class TestThermoDataStandard(unittest.TestCase):

	def setUp(self):
		self.step = [10,20,30]
		self.temp = [20.5,21.6,22.7]
		self.pressure = [10,11,12]
		self.createTestObjs()

	def createTestObjs(self):
		self.inpKwargDictA = {"step":self.step, "temp":self.temp, "pressure":self.pressure}
		self.testObjA = tCode.ThermoDataStandard(self.inpKwargDictA)

	def testExpectedPropsReturned(self):
		expProps = self.inpKwargDictA.keys()
		actProps = self.testObjA.props
		for key in expProps:
			self.assertTrue(key in actProps)

	def testExpectedArrayReturned_twoProps(self):
		props = ["step", "temp"]
		expVals = np.array( (self.step,self.temp) ).transpose()
		actVals = self.testObjA.getPropsArray(props)
		self.assertTrue(np.allclose(expVals,actVals))

	def testExpectedArrayReturned_threeProps(self):
		props = ["step","temp","pressure"]
		expVals = np.array( (self.step,self.temp,self.pressure) ).transpose()
		actVals = self.testObjA.getPropsArray(props)
		self.assertTrue( np.allclose(expVals,actVals) )

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = copy.deepcopy(self.testObjA)
		self.assertEqual(objA,objB)
		
	def testUnequalObjsCompareUnequal_diffPropList(self):
		objA = tCode.ThermoDataStandard({"keyA":[1,2],"keyB":[3,4]})
		objB = tCode.ThermoDataStandard({"keyA":[1,2],"keyC":[3,4]})
		self.assertNotEqual(objA,objB)

	def testUnequalObjsCompareUnequal_extraPropOnSecond(self):
		objA = tCode.ThermoDataStandard({"keyA":[1,2],"keyB":[3,4]})
		objB = tCode.ThermoDataStandard({"keyA":[1,2],"keyB":[3,4], "keyC":[5,6]})
		self.assertNotEqual(objA,objB)

	def testUnequalObjsCompareUnequal_diffNumbs(self):
		objA = copy.deepcopy(self.testObjA)
		self.step[-1] += 1
		self.createTestObjs()
		objB = copy.deepcopy(self.testObjA)
		self.assertNotEqual(objA,objB)

	def testFromSelectedStdKwargs(self):
		objA = tCode.ThermoDataStandard({"step":[1,2], "temp":[3,4]})
		objB = tCode.ThermoDataStandard.fromStdKwargs(temp=[3,4], step=[1,2])
		self.assertEqual(objA,objB)

	def testToAndFromDictConsistent(self):
		objA = self.testObjA
		dictA = objA.toDict()
		objB = tCode.ThermoDataStandard.fromDict(dictA)
		self.assertEqual(objA,objB)


class TestThermoDataStandard_listLengthsChecker(unittest.TestCase):

	def setUp(self):
		self.steps = [1,2,3]
		self.time = [4,5,6]
		self.eTotal = [7,8,9]
		self.createTestObjs()

	def createTestObjs(self):
		self.propDict = {"step":self.steps, "time":self.time, "eTotal":self.eTotal}
		self.testObjA = tCode.ThermoDataStandard(self.propDict)

	def testReturnsTrueWhenAllListsEqual(self):
		self.assertTrue( self.testObjA.dataListLengthsAllEqual )

	def testReturnsFalseWhenOneListLonger(self):
		self.steps.append(4)
		self.createTestObjs()
		self.assertFalse( self.testObjA.dataListLengthsAllEqual )



class TestMergingThermoDataStandard(unittest.TestCase):

	def setUp(self):
		self.stepsA = [0,5]
		self.stepsB = [10,15]
		self.tempA = [200,210]
		self.tempB = [220,230]
		self.overlapStrat = None
		self.createTestObjs()

	def createTestObjs(self):
		kwargDictA = {"step":self.stepsA, "temp":self.tempA}
		kwargDictB = {"step":self.stepsB, "temp":self.tempB}
		self.thermoObjA = tCode.ThermoDataStandard(kwargDictA)
		self.thermoObjB = tCode.ThermoDataStandard(kwargDictB)
		self.dataList = [self.thermoObjA, self.thermoObjB]

	def testExpectedFromMergingOrderedStructs(self):
		expKwargDict = {"step": self.stepsA + self.stepsB,
		                "temp": self.tempA + self.tempB}
		expObj = tCode.ThermoDataStandard( expKwargDict )
		actObj = tCode.getMergedStandardThermoData( self.dataList, overlapStrat=self.overlapStrat )
		self.assertEqual(expObj,actObj)

	def testRaisesIfPropsNotTheSame(self):
		self.dataList[1].dataDict["fake_prop"] = list()
		with self.assertRaises(ValueError):
			tCode.getMergedStandardThermoData( self.dataList, overlapStrat=self.overlapStrat )

	def testRaisesIfDataListsNotAllTheSame(self):
		self.stepsB.append(20)
		with self.assertRaises(AssertionError):
			tCode.getMergedStandardThermoData( self.dataList, overlapStrat=self.overlapStrat )

	def testExpectedStepsWithOverhang_trimStratSimple(self):
		self.stepsA = [1,2,3,4,5]
		self.stepsB = [4,6] #Deleting step 5 makes it more likely we really are taking the step 4 from this set of trajs
		self.tempA = [2*x for x in self.stepsA]
		self.tempB = [3*x for x in self.stepsB]
		self.createTestObjs()

		expKwargDict = {"step": self.stepsA[:3] + self.stepsB,
		                "temp": self.tempA[:3]  + self.tempB }
		expObj = tCode.ThermoDataStandard( expKwargDict )
		actObj = tCode.getMergedStandardThermoData( self.dataList, overlapStrat=self.overlapStrat )
		self.assertEqual(expObj, actObj)





