#!/usr/bin/python3

import itertools as it
import types
import unittest

import gen_basis_helpers.job_utils.convergers as tCode

class TestConvRunnerComposite(unittest.TestCase):

	def setUp(self):
		self.convRunnerStubA = createConvRunnerStubA()
		self.convRunnerStubB = createConvRunnerStubB()
		self.testCompositeA = tCode.PropConvJobRunComposite([self.convRunnerStubA, self.convRunnerStubB])

	def createSimpleComposite(self):
		self.testComposite = tCode.PropConvJobRunComposite([self.convRunnerStubA, self.convRunnerStubB])

	def testGetRunCommsCaseA(self):
		self.createSimpleComposite()
		expectedOutput = [1,2,3,4]
		actualOutput = self.testComposite.runComms
		self.assertEqual(expectedOutput, actualOutput)	

	def testGetRunCommsWithOneCalcJobsFalse(self):
		self.convRunnerStubA.shouldWeRunCalcs = False
		self.createSimpleComposite()
		expectedOutput = [3,4]
		actualOutput = self.testComposite.runComms
		self.assertEqual(expectedOutput, actualOutput)


def createConvRunnerStubA():
	outDict = { "shouldWeRunCalcs": True,
	            "runComms": [1,2] }

	return types.SimpleNamespace( **outDict ) 


def createConvRunnerStubB():
	outDict = {"shouldWeRunCalcs": True,
	           "runComms": [3,4] }

	return types.SimpleNamespace( **outDict )



class TestConvAnalyserComposite(unittest.TestCase):

	def setUp(self):
		self.convAnalyserStubA = createConvAnalyserStubA()
		self.convAnalyserStubB = createConvAnalyserStubB()

	def createAnalyserComposite(self):
		self.testAnalyser = tCode.PropConvAnalyserComposite( [self.convAnalyserStubA, self.convAnalyserStubB] )

	def testGetData(self):
		self.createAnalyserComposite()
		expData = [ [5,6], [6,7], [8,9] ]
		actData = self.testAnalyser.data
		for exp,act in it.zip_longest(expData,actData):
			self.assertTrue( exp, act )


def createConvAnalyserStubA():
	outDict = { "data": [[5,6], [6,7]],
	            "plotData": [7] }
	return types.SimpleNamespace( **outDict )

def createConvAnalyserStubB():
	outDict = {"data": [ [8,9] ],
	           "plotData": [10] }
	return types.SimpleNamespace( **outDict )


if __name__ == '__main__':
	unittest.main()

