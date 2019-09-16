#!/usr/bin/python3

import os
import unittest
import unittest.mock as mock

import sys
sys.path.append('..')

import gen_basis_helpers.job_utils.interstit_helpers as tCode #WARNING: tCode alias doesnt work for the patch decorator


def simpleMockedFunct(*args,**kwargs):
	response = mock.Mock()
	return response


class TestInterstitialTypeClass(unittest.TestCase):

	def setUp(self):
		self.relaxed = "relaxed_constant_p"
		self.structType = "hcp"
		self.interstitType = "octahedral"
		self.cellDims = [3,2,1]
		self.runCalcs = False
		self.createObj()

	def createObj(self):
		self.testObj = tCode.InterstitialType(self.relaxed, self.structType,self.interstitType, self.cellDims, self.runCalcs)

	def testCorrectQueryDefData(self):
		expQueryArgs = [self.structType, self.interstitType, self.relaxed, "3_2_1"]
		mockDb = mock.Mock()
		self.testObj.getInterStructFromRefDataStruct(mockDb)
		mockDb.getSelfInterstitialPlaneWaveStruct.assert_called_with(*expQueryArgs)

@mock.patch('gen_basis_helpers.job_utils.interstit_helpers.supCell.superCellFromUCell')
@mock.patch('gen_basis_helpers.job_utils.interstit_helpers.wFlowCoordinator.WorkFlowCoordinator')
@mock.patch('gen_basis_helpers.job_utils.interstit_helpers.createWFlows.CreateInterstitialWorkFlow')
class TestInterstitialCalcFactory(unittest.TestCase):

	def setUp(self):
		self.testInterstitObjA = createInterstitTypeObjA()
		self.testInterstitObjB = createInterstitTypeObjB()
		self.allInterstitInfoObjs = [self.testInterstitObjA, self.testInterstitObjB]
		self.optDict = {"dataset":"fake_dataset"}
		self.platoComm = "tb1"
		self.label = "test_label"
		self.nCores = 4
		self.workFolder = os.path.abspath(os.getcwd())
		self.refDataObj = mock.Mock()
		self.createFactory()

	def createFactory(self):
		self.testFactory = tCode.InterstitialCalcSetFactory(self.workFolder, self.platoComm, self.optDict, self.allInterstitInfoObjs, self.refDataObj, self.nCores, self.label)

	def testExpectedNumberOfCallsForIndividualWorkflows(self, workFlowMocked:"Makes decorator work", wFlowCoordMocked, supCellMocked):
		expCallCount = 2
		self.testFactory()
		self.assertTrue(workFlowMocked.call_count==expCallCount)



#This is a test for the leaf structure that runs calcs
class TestInterstitialCalculationsSet(unittest.TestCase):

	def setUp(self):
		self.wFlowCoord = mock.Mock()
		self.label = "test_label"
		self.runJobs = True
		self.createObj()

	def createObj(self):
		self.testObj = tCode.InterstitialCalculationSet(self.wFlowCoord, self.label, self.runJobs)

	def testCorrectRunCmdWhenRunCalcsTrue(self):
		self.testObj.runCalcs()
		self.wFlowCoord.run.assert_called_with(inclPreRun=True)

	def testCorrectRunCmdWhenCalcsFalse(self):
		self.runJobs = False
		self.createObj()
		self.testObj.runCalcs()
		self.wFlowCoord.run.assert_not_called()



def createInterstitTypeObjA():
	relaxed = "relaxed_constant_p"
	structType = "hcp"
	interstitType = "octahedral"
	cellDims = [3,2,1]
	runCalcs = False
	testObj = tCode.InterstitialType(relaxed, structType,interstitType, cellDims, runCalcs)
	return testObj


def createInterstitTypeObjB():
	relaxed = "unrelaxed"
	structType = "bcc"
	interstitType = "octahedral"
	cellDims = [3,2,2]
	runCalcs = False
	testObj = tCode.InterstitialType(relaxed, structType,interstitType, cellDims, runCalcs)
	return testObj



if __name__ == '__main__':
	unittest.main()

