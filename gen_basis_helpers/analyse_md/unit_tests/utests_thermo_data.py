
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

