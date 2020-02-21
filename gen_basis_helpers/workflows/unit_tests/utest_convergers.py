
import itertools as it
import types

import unittest
import unittest.mock as mock

import gen_basis_helpers.workflows.convergers as tCode


class TestConvergerTemplate(unittest.TestCase):

	def setUp(self):
		self.runComms = ["fake_comm_a", "fake_comm_b"]
		self.convVals = [2,4]
		self.attrAVals = ["something","else"]
		self.mapFunction = _stubMapFunctionGridValsVsAttrA
		self.namespaceAttrs = ["fake_attr_a"]
		self.parsedFileA = types.SimpleNamespace( attrA=self.attrAVals[0] )
		self.parsedFileB = types.SimpleNamespace( attrA=self.attrAVals[1] )
		self.createTestObjs()

	def createTestObjs(self):
		self.mockCalcMethObjA = _createMockCalcMethObj(self.runComms[0], self.parsedFileA)
		self.mockCalcMethObjB = _createMockCalcMethObj(self.runComms[1], self.parsedFileB)
		self.calcObjs = [self.mockCalcMethObjA, self.mockCalcMethObjB]
		self.testWorkflowA = tCode.ConvergerWorkflowTemplate( self.calcObjs, self.convVals, self.mapFunction, self.namespaceAttrs )

	def testTemplateConvergerPreShellRunComms(self):
		expComms = self.runComms
		actComms = self.testWorkflowA.preRunShellComms
		self.assertEqual(expComms,actComms)

	def testTemplateConvergerNamespaceAttrs(self):
		self.assertEqual( self.namespaceAttrs, self.testWorkflowA.namespaceAttrs )

	def testRunCommand(self):
		self.testWorkflowA.run()
		expVals = [(x,y) for x,y in it.zip_longest(self.convVals, self.attrAVals)]
		actVals = self.testWorkflowA.output.fake_attr_a
		self.assertEqual(expVals, actVals) 


	def testWriteFileCalledUponInitialisation(self):
		self.mockCalcMethObjA.writeFile.assert_called_once_with()
		self.mockCalcMethObjB.writeFile.assert_called_once_with()

def _createMockCalcMethObj( runComm, parsedFile ):
	outObj = mock.Mock()
	outObj.runComm = runComm
	outObj.parsedFile = parsedFile
	return outObj

def _stubMapFunctionGridValsVsAttrA(workflow):
	parsedFiles = list()
	for x in workflow.calcObjs:
		parsedFiles.append(x.parsedFile)

	convVals = [x for x in workflow.convVals]
	attrAVals = [x.attrA for x in parsedFiles]

	attrAOutput = [(x,y) for x,y in it.zip_longest(convVals,attrAVals)]
	workflow.output = types.SimpleNamespace(fake_attr_a=attrAOutput)


class TestGridConvergenceWorkflow(unittest.TestCase):

	def setUp(self):
		self.convVals = [2,4,6]
		self.energies = [1,2,3]
		self.createTestObjects()

	def createTestObjects(self):
		energyObjs = [types.SimpleNamespace(electronicTotalE=x) for x in self.energies]
		parsedFileVals = [types.SimpleNamespace(energies=x) for x in energyObjs]
		calcObjs = [mock.Mock() for x in parsedFileVals]
		for cObj,pVal in it.zip_longest(calcObjs,parsedFileVals):
			cObj.parsedFile = pVal
		self.testObjA = tCode.GridConvergenceEnergyWorkflow(calcObjs, self.convVals)

	def testRunMethod(self):
		self.testObjA.run()
		expOutput = [(x,y) for x,y in it.zip_longest(self.convVals,self.energies)]
		actOutput = self.testObjA.output.convResults
		self.assertEqual(expOutput,actOutput)

