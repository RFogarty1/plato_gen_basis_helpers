#!/usr/bin/python3
import types
import unittest

import numpy as np

import gen_basis_helpers.elemental_eos.multi_cryst_eos as tCode

class TestSingleCrystEos(unittest.TestCase):

	def setUp(self):
		self.testObjA = createSingleCrystEosStubA()

	def testPlotDataDictProperty(self):
		expDict = {"tstructA": np.array( ([1,3],[2,6]) )}
		actDict = self.testObjA.getPlotData()
		for key in expDict:
			self.assertTrue( np.allclose(expDict[key],actDict[key]) )

	def testE0(self):
		expVal = {"tstructA":4}
		actVal = self.testObjA.e0
		self.assertEqual(expVal,actVal)

	def testGetTableData(self):
		expDict = {"tstructA": ["Exact_Dft","1.000","2.000","4.000"]}
		actDict = self.testObjA.getTableData()
		for key in expDict:
			self.assertEqual(expDict[key],actDict[key])

def createSingleCrystEosStubA():
	v0, b0, e0 = 1,2,4
	fitData = np.array( ([1,2],[2,4]) )
	actData = np.array( ([1,3],[2,6]) )
	structLabel = "tstructA"
	elementLabel = "Mg"
	methodLabel = "Exact_Dft"
	return tCode.SingleCrystEosResult(v0=v0,b0=b0,e0=e0,fitData=fitData,actData=actData,structLabel=structLabel,elementLabel=elementLabel,methodLabel=methodLabel)

class TestMultiCrystEos(unittest.TestCase):

	def setUp(self):
		self.testLabelA = "testLabel"
		self.singCrystA = createSingleCrystForCompositeTestStubA()
		self.singCrystB = createSingleCrystForCompositeTestStubB()
		self.createTestObj()

	def createTestObj(self):
		self.testObjA = tCode.MultiCrystEosResult([self.singCrystA,self.singCrystB])

	def testDeltaE0(self):
		expDict = {"testStructA": 8,
		           "testStructB": 0}
		actDict = self.testObjA.deltaE0
		self.assertEqual(expDict,actDict)

	def testInitiationWithDuplicateStructLabels(self):
		self.singCrystB.label[0].structKey = self.singCrystA.label[0].structKey
		with self.assertRaises(ValueError):
			self.createTestObj()

	def testInitiationWithVaryingMethodLabels(self):
		self.singCrystB.label[0].methodKey = self.singCrystB.label[0].methodKey + "_now_different"
		with self.assertRaises(ValueError):
			self.createTestObj()

	def testInitialisationWithVaryingElementLabels(self):
		self.singCrystB.label[0].eleKey = self.singCrystA.label[0].eleKey + "now_different"
		with self.assertRaises(ValueError):
			self.createTestObj()

	def testGetTableData(self):
		expTableData = {"testStructA": ["fake_method","1.000","2.000","3.000"],
		                "testStructB": ["fake_method","4.000","5.000","6.000"]}
		actTableData = self.testObjA.getTableData()
		self.assertEqual(expTableData, actTableData)

	def testGetPlotData(self):
		expPlotData = {"testStructA": np.array( ([1,2],[2,4]) ),
		               "testStructB": np.array( ([1,4],[2,6]) )}
		actPlotData = self.testObjA.getPlotData()
		for key in expPlotData.keys():
			self.assertTrue( np.allclose(expPlotData[key],actPlotData[key]) )


def createSingleCrystForCompositeTestStubA():
	eleLabel, structLabel, methodLabel = "Zr", "testStructA", "fake_method"
	fakeLabel = tCode.StandardLabel(eleKey=eleLabel, structKey=structLabel, methodKey=methodLabel)
	fakeTableData = {structLabel: [1,2,3]}
	fakePlotData = {structLabel: np.array( ([1,2],[2,4]) )}

	def tableDataFunct(numbDp=None):
		return {"testStructA": ["fake_method"] + ["{:.3f}".format(x) for x in [1,2,3]]}

	outDict = { "e0": {structLabel:22},
	            "structLabels": [structLabel],
	            "methodLabel": methodLabel,
	            "elementLabel": eleLabel,
	            "getTableData": tableDataFunct,
	            "getPlotData": lambda: fakePlotData,
	            "label": [fakeLabel] }
	outObj = types.SimpleNamespace(**outDict)
	return outObj

def createSingleCrystForCompositeTestStubB():
	eleLabel, structLabel, methodLabel = "Zr", "testStructB", "fake_method"
	fakeLabel = tCode.StandardLabel(eleKey=eleLabel, structKey=structLabel, methodKey=methodLabel)
	fakeTableData ={structLabel: [4,5,6]}
	fakePlotData = {structLabel: np.array( ([1,4],[2,6]) )}

	def tableDataFunct(numbDp=None):
		return {"testStructB": ["fake_method"] + ["{:.3f}".format(x) for x in [4,5,6]]}

	outDict = {"e0": {structLabel:14},
	           "structLabels": [structLabel],
	           "methodLabel": methodLabel,
	           "elementLabel": eleLabel,
	           "getTableData": tableDataFunct,
	           "getPlotData": lambda: fakePlotData,
	           "label": [fakeLabel] }
	outObj = types.SimpleNamespace(**outDict)
	return outObj

if __name__ == '__main__':
	unittest.main()

