
import os

import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.elastic_constants as tCode
import gen_basis_helpers.shared.label_objs as labelHelp

class TestStandardInputFactory(unittest.TestCase):

	def setUp(self):
		self.baseWorkFolder = "fake_folder"
		self.extToWorkFolder = "elastic_hcp" #The original default value
		self.baseGeom = mock.Mock()
		self.creator = mock.Mock()
		self.eleKey, self.structKey, self.methodKey = mock.Mock(), mock.Mock(), mock.Mock()
		self.strainValues = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.HcpElasticStandardInputCreator(baseGeom=self.baseGeom, creator=self.creator,
		                                                     strainValues=self.strainValues, eleKey=self.eleKey,
		                                                     structKey=self.structKey, methodKey=self.methodKey,
		                                                     baseWorkFolder=self.baseWorkFolder)

	@mock.patch("gen_basis_helpers.job_helpers.elastic_constants.elasticFlow.HcpElasticWorkflowCreator")
	def testExpectedArgsPassedToWorkflowFactory(self, mockedWorkFlowFactory):
		expArgDict = {"baseGeom":self.baseGeom, "strainValues":self.strainValues,
		              "creator":self.creator}
		expArgDict["workFolder"] = os.path.join(self.baseWorkFolder, self.extToWorkFolder)

		self.testObjA.create()
		mockedWorkFlowFactory.assert_called_once_with(**expArgDict)

	def testExpectedLabelCreated(self):
		expLabel = labelHelp.StandardLabel(eleKey=self.eleKey, methodKey=self.methodKey, structKey=self.structKey)
		actLabel = self.testObjA._createLabel()
		self.assertEqual(expLabel,actLabel)

	@mock.patch("gen_basis_helpers.job_helpers.elastic_constants.calcRunners.StandardInputObj")
	@mock.patch("gen_basis_helpers.job_helpers.elastic_constants.HcpElasticStandardInputCreator._createWorkFlow")
	def testExpectedArgsPassedToStandardInput(self, mockedWorkflowCreator, mockedStandardInput):
		expWorkflow = "fake_workflow"
		expFinalObject = "fake_standard_input"
		expLabel = labelHelp.StandardLabel(eleKey=self.eleKey, methodKey=self.methodKey, structKey=self.structKey)

		mockedWorkflowCreator.side_effect = lambda : expWorkflow
		mockedStandardInput.side_effect = lambda *args, **kwargs : expFinalObject

		actFinalObject = self.testObjA.create()
		mockedStandardInput.assert_called_once_with(expWorkflow, expLabel)
	
		self.assertEqual(expFinalObject, actFinalObject)


