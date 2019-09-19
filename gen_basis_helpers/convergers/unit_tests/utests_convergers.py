#!/usr/bin/python3

import itertools as it
import types
import unittest
import unittest.mock as mock

import numpy as np

import gen_basis_helpers.convergers.convergers as tCode

class TestConvRunnerComposite(unittest.TestCase):

	def setUp(self):
		self.convRunnerStubA = createConvRunnerStubA()
		self.convRunnerStubB = createConvRunnerStubB()
		self.createSimpleComposite()

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

	def testWorkFolderAsExpected(self):
		self.convRunnerStubA.workFolders.append(None)
		self.createSimpleComposite()
		expWorkFolders = ["testA","testB"]
		actWorkFolders = self.testComposite.workFolders
		self.assertEqual(expWorkFolders, actWorkFolders)


	def testCantCreateWithDuplicateWorkFolders(self):
		self.convRunnerStubB.workFolders.append("testA")
		with self.assertRaises(tCode.DuplicateWorkFoldersError):
			self.createSimpleComposite()


def createConvRunnerStubA():
	outDict = { "shouldWeRunCalcs": True,
	            "runComms": [1,2],
	            "workFolders": ["testA"] }

	return types.SimpleNamespace( **outDict ) 


def createConvRunnerStubB():
	outDict = {"shouldWeRunCalcs": True,
	           "runComms": [3,4],
	           "workFolders":["testB"] }

	return types.SimpleNamespace( **outDict )



class TestConvAnalyserComposite(unittest.TestCase):

	def setUp(self):
		self.convAnalyserStubA = createConvAnalyserStubA()
		self.convAnalyserStubB = createConvAnalyserStubB()
		self.createAnalyserComposite()

	def createAnalyserComposite(self):
		self.testAnalyser = tCode.PropConvAnalyserComposite( [self.convAnalyserStubA, self.convAnalyserStubB] )

	def testGetData(self):
		expData = [ [5,6], [6,7], [8,9] ]
		actData = self.testAnalyser.data
		for exp,act in it.zip_longest(expData,actData):
			self.assertTrue( exp, act )

	def testPlotData(self):
		self.convAnalyserStubA.plotData.return_value = list()
		self.convAnalyserStubB.plotData.return_value = list()
		self.testAnalyser.plotData()
		self.convAnalyserStubA.plotData.assert_called_once_with()
		self.convAnalyserStubB.plotData.assert_called_once_with()


def createConvAnalyserStubA():
	outDict = { "data": [[5,6], [6,7]],
	            "plotData": mock.Mock() }
	return types.SimpleNamespace( **outDict )

def createConvAnalyserStubB():
	outDict = {"data": [ [8,9] ],
	           "plotData": mock.Mock() }
	return types.SimpleNamespace( **outDict )


class TestBasicConvRunnerImplementation(unittest.TestCase):

	def setUp(self):
		self.fakeWorkFlowA = createStubWorkFlowA()
		self.fakeWorkFlowB = createStubWorkFlowB()
		self.shouldWeRunCalcs = True
		self.testLabelA = "labelA"
		self.varyParams = [1,2]
		self.createTestObject()

	def createTestObject(self):
		wFlowList = [self.fakeWorkFlowA, self.fakeWorkFlowB]
		self.testObj = tCode.PropConvJobRunnerStandard(wFlowList, self.varyParams, self.testLabelA, shouldWeRunCalcs=self.shouldWeRunCalcs)

	def testCantInitWithDuplWorkFolders(self):
		self.fakeWorkFlowA.workFolder = self.fakeWorkFlowB.workFolder
		with self.assertRaises(tCode.DuplicateWorkFoldersError):
			self.createTestObject()

	def testRunComms(self):
		expRunComms = ["commA","commB","commC"]
		actRunComms = self.testObj.runComms
		self.assertEqual(expRunComms, actRunComms)

	@mock.patch("gen_basis_helpers.convergers.convergers.jobRun.executeRunCommsParralel")
	def testRunJobsWhenShouldWeRunCalcsTrue(self, jobRunnerMock):
		expRunComms = ["commA","commB","commC"]
		self.testObj.runJobs()
		jobRunnerMock.assert_called_once_with(expRunComms)

	@mock.patch("gen_basis_helpers.convergers.convergers.jobRun.executeRunCommsParralel")
	def testRunJobsWhenShouldWeRunCalcsFalse(self, jobRunnerMock):
		self.testObj.shouldWeRunCalcs = False
		self.testObj.runJobs()
		jobRunnerMock.assert_not_called()

	@mock.patch("gen_basis_helpers.convergers.convergers.PropConvAnalyserStandard")
	def testCreateAnalyser(self, analyserMock):
		#Create a stub representing the state after run has been applied
		self.fakeWorkFlowA.output = types.SimpleNamespace(fake_output=22)
		self.fakeWorkFlowB.output = types.SimpleNamespace(fake_output=13)
		self.fakeWorkFlowA.run = mock.Mock()
		self.fakeWorkFlowB.run = mock.Mock()
		self.createTestObject()

		outAnalyser = self.testObj.createAnalyser()
		expConvData = [(1,22),(2,13)]
		analyserMock.assert_called_once_with(expConvData, self.testLabelA)
		self.fakeWorkFlowA.run.assert_any_call()
		self.fakeWorkFlowB.run.assert_any_call()

def createStubWorkFlowA():
	outDict = {"workFolder":"fakeA",
	           "namespaceAttrs":["energy"],
	           "preRunShellComms":["commA","commB"]}
	outFlow = types.SimpleNamespace( **outDict )
	return outFlow

def createStubWorkFlowB():
	outDict = {"workFolder":"fakeB",
	           "namespaceAttrs":["volume","b0"],
	           "preRunShellComms":["commC"]}
	outFlow = types.SimpleNamespace( **outDict )
	return outFlow



class TestBasicJobAnalyserImplementation(unittest.TestCase):

	def setUp(self):
		self.fakeDataA = np.array([(1,2), (3,4)])
		self.fakeLabelA = "labelA"
		self.dataPlotterA = mock.Mock()
		self.createTestObj()

	def createTestObj(self):
		self.testObj = tCode.PropConvAnalyserStandard(self.fakeDataA, self.fakeLabelA, dataPlotter=self.dataPlotterA)

	@mock.patch("gen_basis_helpers.convergers.convergers.PropConvAnalyserStandard.data")
	def testPlotDataFunct(self,fakeData):
		self.testObj.plotData()
		self.dataPlotterA.createPlot.assert_called_once_with( fakeData, titleStr=self.fakeLabelA )

	def testGetDataRelativeToMaxConvVal(self):
		self.testObj.useDataRelativeToLargestConvValue = True
		expData = [np.array( ([1,-2], [3,0]) )]
		actData = self.testObj.data

		print("expData:")
		print(expData)
		print("actData:")
		print(actData)
		for exp,act in it.zip_longest(expData,actData):
			self.assertTrue( np.allclose(exp,act) )



if __name__ == '__main__':
	unittest.main()

