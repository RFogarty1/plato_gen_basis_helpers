
import os
import itertools as it
import types

from collections import OrderedDict

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
		self.eType = None
		self.strainValues = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.HcpElasticStandardInputCreator(baseGeom=self.baseGeom, creator=self.creator,
		                                                     strainValues=self.strainValues, eleKey=self.eleKey,
		                                                     structKey=self.structKey, methodKey=self.methodKey,
		                                                     baseWorkFolder=self.baseWorkFolder, eType=self.eType)

	@mock.patch("gen_basis_helpers.job_helpers.elastic_constants.elasticFlow.HcpElasticWorkflowCreator")
	def testExpectedArgsPassedToWorkflowFactory(self, mockedWorkFlowFactory):
		expArgDict = {"baseGeom":self.baseGeom, "strainValues":self.strainValues,
		              "creator":self.creator, "eType":self.eType}
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



class TestMapElasticToUsefulFormatHcpCase(unittest.TestCase):

	def setUp(self):
		self.structStr = "hcp"

		self.fitFunctXVals = [-2,1,0,1,2]
		self.elasticConstants = [x for x in range(1,6)]
		self.eleKey = "fake_ele_key"
		self.methodKey = "fake_method_key"
		self.structKey = "fake_struct_key"
		self.expHeadings = ["method","c11","c12","c13","c33","c44"]
		self.actPlotData = [mock.Mock() for x in self.elasticConstants]
		self.fitVals = [mock.Mock() for x in self.elasticConstants]
		self.strainObjs = [mock.Mock() for x in self.elasticConstants] #One strain per elastic obviously
		self.runMethod = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		labelA = labelHelp.StandardLabel(eleKey=self.eleKey,structKey=self.structKey,methodKey=self.methodKey)
		elasticStrs = ["c11","c12","c13","c33","c44"]
		elasticsInCorrectFormat = OrderedDict( [(strLabel,y) for strLabel,y in it.zip_longest(elasticStrs,self.elasticConstants)] )
		stressStrainData = [mock.Mock() for x in self.elasticConstants]
		self.testObjA = tCode.MapElasticflowOutputToUsefulFormatStandard(self.structStr, fitStrainVals=self.fitFunctXVals)
		self.testOutputA = types.SimpleNamespace(**{"elasticConsts":elasticsInCorrectFormat,"stressStrainData":stressStrainData,
		                                            "strains":self.strainObjs})
		testWorkflowA = types.SimpleNamespace(output=[self.testOutputA], run=self.runMethod)
		self.standardInpObjA = types.SimpleNamespace(workflow=testWorkflowA, label=[labelA])

		#Need to mock the fitFunction to retunr an iterable
		for idx,x in enumerate(self.standardInpObjA.workflow.output[0].stressStrainData):
			x.fitFunct.side_effect = lambda *args,**kwargs: [self.fitVals[idx]] #Just needs to be iterable

		#Also need to mock the actData to return self.actPlotData
		for idx,x in enumerate(self.standardInpObjA.workflow.output[0].stressStrainData):
			x.actVals = self.actPlotData[idx] #Mahybe sohuld be an iter or something?

	def runMapFunctOnDataA(self):
		return self.testObjA( self.standardInpObjA )

	def testCanOnlyCreateForHcpForNow(self):
		self.structStr = "some_fake_struct"
		with self.assertRaises(ValueError):
			self.createTestObjs()
	
	def testRunMethodIsCalled(self):
		self.runMapFunctOnDataA()
		self.runMethod.assert_called_once_with()

	def testExpTableDataPresent(self):
		outputObj = self.runMapFunctOnDataA()
		expTableData = [self.methodKey] + ["{:.2f}".format(x) for x in self.elasticConstants]
		actTableData = outputObj.tableData
		self.assertEqual(expTableData,actTableData)

	def testStrainStrsAsExpected(self):
		outputObj = self.runMapFunctOnDataA()
		expStrainStrs = [x.toStr() for x in self.strainObjs]
		actStrainStrs = outputObj.strainStrs
		self.assertEqual(expStrainStrs, actStrainStrs)

	def testCorrectValsPassedToFitVals(self):
		for x in self.standardInpObjA.workflow.output[0].stressStrainData:
			x.fitFunct.side_effect = lambda *args,**kwargs: [mock.Mock()] #Just needs to be iterable

		#Modify fit function mock to support iteration
		self.runMapFunctOnDataA()
		for x in self.standardInpObjA.workflow.output[0].stressStrainData:
			x.fitFunct.assert_called_once_with(self.fitFunctXVals)



	@unittest.skip("")
	def testExpectedFitStrainValsGeneratedAsDefault(self):
		self.fitFunctXVals = None
		self.standardInpObjA.workflow.strainValues = self.fitFunctXVals
		self.createTestObjs()
		self.runMapFunctOnDataA()


	@unittest.skip("Couldnt get to work properly; so sorta abandoning for now")
	def testExpectedPlotData(self):
		expFitData = [[[x,y]] for x,y in it.zip_longest(self.fitFunctXVals,self.fitVals)]
		expPlotData = [ [actData, fitData] for actData, fitData in it.zip_longest(self.actPlotData,expFitData) ]
		actOutput = self.runMapFunctOnDataA()
		actPlotData = actOutput.plotData
		self.assertEqual(expPlotData, actPlotData)




