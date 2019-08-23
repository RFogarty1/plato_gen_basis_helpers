#!/usr/bin/python3

import copy
import itertools as it
import os
import types
import unittest
import unittest.mock as mock

import sys
sys.path.append('..')
import inv_sk_workflows as tCode


class TestFetchForComposite(unittest.TestCase):
	
	def setUp(self):
		self.mockCompositeObj = createMockCompositeObjA()
		self.createTestObj()

	def createTestObj(self):
		self.testObj = tCode.CompositeInvSkWorkFlow([self.mockCompositeObj],label="top-level")

	def testCorrectInputLeaf(self):
		fetchArgs = ["Zr","perfectCrystals","hcp"]
		expLabel = "hcp"
		outObj = self.testObj.fetch(*fetchArgs)
		self.assertEqual(expLabel,outObj.label)

	def testCorrectInputNode(self):
		fetchArgs = ["Zr", "perfectCrystals"]
		expLabel = "perfectCrystals"
		outObj = self.testObj.fetch(*fetchArgs)
		self.assertEqual(expLabel,outObj.label)


class TestCompositeNamespaceAttrs(unittest.TestCase):

	def setUp(self):
		self.mockCompositeObj = createMockCompositeObjA()
		self.createTestObj()

	def createTestObj(self):
		self.testObj = tCode.CompositeInvSkWorkFlow([self.mockCompositeObj],label="top-level")

	def testExpectedNamespaceAttrs(self):
		expNamespaceAttrs = ["fake_generic","fake_ver_two"]
		actNamespaceAttrs = self.testObj.namespaceAttrs
		self.assertEqual(expNamespaceAttrs, actNamespaceAttrs)

def createMockCompositeObjA():
	standardNamespaceAtt = "fake_generic"
	altNamespaceAtt = "fake_ver_two"
	leafHcpPerfect = types.SimpleNamespace(label="hcp",namespaceAttrs=[standardNamespaceAtt]) #Only property leafs need for this
	leafBccPerfect = types.SimpleNamespace(label="bcc",namespaceAttrs=[standardNamespaceAtt])
	interstitialLeaf = types.SimpleNamespace(label="interstit", namespaceAttrs = [altNamespaceAtt])

	perfectCrystalComposite = tCode.CompositeInvSkWorkFlow([leafHcpPerfect,leafBccPerfect],"perfectCrystals")
	topLevel = tCode.CompositeInvSkWorkFlow([perfectCrystalComposite,interstitialLeaf],"Zr")

	return topLevel



class TestCompositeSkWorkFlow(unittest.TestCase):

	def setUp(self):
		self.mockWorkFlowA = mock.Mock()
		self.mockWorkFlowB = mock.Mock()
		self.allMockWorkFlows = [self.mockWorkFlowA, self.mockWorkFlowB]
		self.label = "test_label"
		self.createTestObj()

	def createTestObj(self):
		self.testObj = tCode.CompositeInvSkWorkFlow(self.allMockWorkFlows, self.label)

	def testPreRunCommsTwoVals(self):
		self.mockWorkFlowA.preRunShellComms = ["commA"]
		self.mockWorkFlowB.preRunShellComms = ["commB"]
		expComms = ["commA","commB"]
		actComms = self.testObj.preRunShellComms 
		self.assertEqual(expComms,actComms)

	def testPreRunCommsOneNone(self):
		self.mockWorkFlowA.preRunShellComms = ["commA"]
		self.mockWorkFlowB.preRunShellComms = None
		expComms = ["commA"]
		actComms = self.testObj.preRunShellComms
		self.assertEqual(expComms,actComms)



@mock.patch('inv_sk_workflows.InvSkWorkFlow._writeFiles')
class TestInvSkWorkFlow(unittest.TestCase):

	def setUp(self):
		self.structs = [None,None,None]
		self.optDict = {"dataset":"fake_data"}
		self.workFolder = "fake_relative_dir"
		self.elements = ["Mg","Zr"]
		self.label  = "useless_label"

	def initWorkFlow(self):
		self.testWorkFlow = tCode.InvSkWorkFlow(self.workFolder, self.optDict, self.elements, self.structs, self.label)

	def testInvSkWorkFlowInitSimpleAttrs(self,writeFilesMock):
		self.initWorkFlow() #Patch only mocks write_files in this for test functions		
		self.workFolder = os.path.abspath(self.workFolder) #we want abs-paths stored in the workFlow
		self.optDict["inversesk"] = 1
		for attr in ["structs","optDict","workFolder","elements"]:
			self.assertEqual( getattr(self,attr), getattr(self.testWorkFlow,attr) )

	def testWriteFilesCalledUponInitialisation(self, writeFilesMock):
		self.initWorkFlow()
		writeFilesMock.assert_called_once_with()

	def testNamespaceAttrs(self,writeFilesMock):
		self.initWorkFlow() 
		expAttrs = ["invSkResults_"+x for x in ["Mg_Mg","Mg_Zr","Zr_Mg","Zr_Zr"]]
		actAttrs = self.testWorkFlow.namespaceAttrs
		[self.assertEqual(exp,act) for exp,act in it.zip_longest(expAttrs,actAttrs)]

	@mock.patch('inv_sk_workflows.jobRun.invSkInputPathsToBashComms')
	def testPreRunShellComms(self,cmdConvMock,writeFilesMock):
		self.initWorkFlow() 
		expComms = "test_cmd_convStr"
		cmdConvMock.return_value = expComms
		actComms = self.testWorkFlow.preRunShellComms
		self.assertEqual(expComms, actComms)
		cmdConvMock.assert_called_once_with(self.testWorkFlow._inpFilePaths)

	def testGetOutFilePathsOneEleCombo(self,writeFilesMock):
		self.initWorkFlow() 		
		eleCombo = ["Mg","Zr"]
		basePath = os.path.abspath(self.workFolder)
		expFileNames = 	["inv_sk_inp_{}_Mg_Zr_SK.csv".format(x) for x in range(len(self.structs))]
		expVals = [os.path.join(basePath,x) for x in expFileNames]
		actVals = self.testWorkFlow.getOutFilePathsOneEleCombo(*eleCombo)
		for exp,act in it.zip_longest(expVals,actVals):
			self.assertEqual(exp,act)

if __name__ == '__main__':
	unittest.main()

