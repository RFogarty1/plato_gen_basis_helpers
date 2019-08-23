#!/usr/bin/python3

import itertools as it
import sys
import unittest

sys.path.append('..')
import ref_elemental_objs as tCode

class TestShellAngMomMapper(unittest.TestCase):

	def setUp(self):
		self.MgShellToAngMom = [2,3,0,1,2]
		self.ZrShellToAngMom = [2,0,1]
		self.eleToAngMomDict = {"MG":self.MgShellToAngMom, "zR": self.ZrShellToAngMom}
		self.createTestObj()

	def createTestObj(self):
		self.testObj = tCode.ShellAngMomMapper(self.eleToAngMomDict)

	def testCorrectAngMomFromIdx(self):
		testArgs = [("mg",2), ("zr",0), ("mg",1)]
		expVals =  [ 0      , 2,        3       ]
		actVals = [self.testObj.getShellIdxToAngMom(*currArgs) for currArgs in testArgs]
		[self.assertEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	def testCorrectSingleIndexFromGivenAngMom(self):
		testEle, testAngMom = "mg",3
		expIndex = 1
		actIndex = self.testObj.getShellIndicesForAngMom(testEle,testAngMom)[0]
		self.assertEqual(expIndex,actIndex)

	def testCorrectMultipleElementsForGivenAngMom(self):
		testEle, testAngMom = "mg", 2
		expIndices = [0,4]
		actIndices = self.testObj.getShellIndicesForAngMom(testEle,testAngMom)
		self.assertEqual(expIndices,actIndices)

	def testGetIndexFromAngMomSpecificZetaSecondVal(self):
		testEle, testAngMom = "mg", 2
		testNthVal = 1
		expShellIndex = 4
		actShellIndex = self.testObj.getSingleShellIndexFromAngMom(testEle, testAngMom, nthVal=1)
		self.assertEqual(expShellIndex,actShellIndex)

	def testGetIndexFromAngMomSpecificZetaRaisesWhenNTooHigh(self):
		testEle, testAngMom = "Zr", 2
		with self.assertRaises(tCode.AngMomMissing):
			self.testObj.getSingleShellIndexFromAngMom(testEle, testAngMom, nthVal=1)

	def testGetIndexFromAngMomSpecificZetaRaisesWhenAngMomMissing(self):
		testEle, testAngMom = "Zr", 4
		with self.assertRaises(tCode.AngMomMissing):
			self.testObj.getSingleShellIndexFromAngMom(testEle,testAngMom, nthVal=0)


	def testShellToAngMomDict(self):
		testEle = "Mg"
		expDict = {0:2, 1:3, 2:0, 3:1, 4:2}
		actDict = self.testObj.getShellToAngMomDict(testEle)
		self.assertEqual(expDict,actDict)

if __name__ == '__main__':
	unittest.main()

