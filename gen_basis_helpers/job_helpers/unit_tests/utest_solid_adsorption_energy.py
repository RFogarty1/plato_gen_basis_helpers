
import itertools as it

import os
import unittest
import unittest.mock as mock


import gen_basis_helpers.shared.label_objs as labelHelp
import gen_basis_helpers.job_helpers.solid_adsorption_energy as tCode

class TestSolidAdsorptionStandardInputCreatorTemplate(unittest.TestCase):

	def setUp(self):
		self.baseWorkFolder = "base_work_folder"
		self.eleKey = "test_ele"
		self.structKey = "test_struct"
		self.methodKey = "test_method"
		self.kPtsBulk = mock.Mock()
		self.geomBulkWithAdsorbed = mock.Mock()
		self.bulkWithoutAdsorbedGeom = mock.Mock()
		self.gasPhaseReactantGeoms = [mock.Mock for x in range(2)]

		self.gasPhaseReactantStoichiometries = [1,2]
		self.gasPhaseProductStoichiometries = [2]

		self.expOutFolder = os.path.join(self.baseWorkFolder, self.eleKey, self.structKey, self.methodKey)
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"eleKey":self.eleKey, "structKey":self.structKey, "methodKey":self.methodKey,
		             "baseWorkFolder":self.baseWorkFolder, "bulkWithAdsorbedGeom":self.geomBulkWithAdsorbed,
		             "kPtsBulk":self.kPtsBulk, "bulkWithoutAdsorbedGeom":self.bulkWithoutAdsorbedGeom,
		             "gasPhaseReactantGeoms":self.gasPhaseReactantGeoms,
		             "gasPhaseReactantStoichiometries":self.gasPhaseReactantStoichiometries}
		self.testObjA = tCode.CodeSpecificStandardInputCreatorTemplate( **kwargDict )

	def testOutFolderAsExpected(self):
		expOutFolder = self.expOutFolder
		actOutFolder = self.testObjA.outFolder
		self.assertEqual(expOutFolder,actOutFolder)

	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.CodeSpecificStandardInputCreatorTemplate._getBaseCreator")
	def testExpectedModsToAdsorbedCalcObj(self, mockedGetBaseCreator):
		baseCreatorObj = mock.Mock()
		mockedGetBaseCreator.side_effect = lambda *args,**kwargs: baseCreatorObj

		self.testObjA._getBulkWithAdsorbedGeomCreator()

		mockedGetBaseCreator.assert_called_once_with()

		expKwargDict = {"workFolder": self.expOutFolder,
		                "fileName": "bulk_with_adsorbed",
		                "kPts":self.kPtsBulk,
		                "geom":self.geomBulkWithAdsorbed} 

		for key in expKwargDict.keys():
			self.assertEqual( expKwargDict[key], getattr(baseCreatorObj,key) )

	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.CodeSpecificStandardInputCreatorTemplate._getBaseCreator")
	def testExpectedModsToBulkCalcObj(self, mockedGetBaseCreator):
		baseCreatorObj = mock.Mock()
		mockedGetBaseCreator.side_effect = lambda *args,**kwargs: baseCreatorObj
		self.testObjA._getBulkWithoutAdsorbedGeomCreator()

		mockedGetBaseCreator.assert_called_once_with()

		expKwargDict = {"workFolder": self.expOutFolder,
		                "fileName": "bulk_calculation",
		                "kPts": self.kPtsBulk,
		                "geom": self.bulkWithoutAdsorbedGeom}

		for key in expKwargDict.keys():
			self.assertEqual( expKwargDict[key], getattr(baseCreatorObj,key) )

	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.CodeSpecificStandardInputCreatorTemplate._getBaseCreator")
	def testExpectedModsForGasPhaseReactants(self, mockedGetBaseCreator):
		mockedGetBaseCreator.side_effect = lambda *args,**kwargs: mock.Mock()

		outCreatorObjs = self.testObjA._getGasPhaseReacantObjs()

		mockedGetBaseCreator.assert_any_call() #Should really test it has the correct number of calls but...
		
		expKwargDictShared = {"workFolder": self.expOutFolder,
		                      "kPts": [1,1,1]}

		fileNames = ["gas_phase_reactant_{}".format(x) for x in range(len(self.gasPhaseReactantGeoms))]

		#Check shared option as expected
		for creatorObj in outCreatorObjs:
			for key in expKwargDictShared.keys():
				self.assertEqual( expKwargDictShared[key], getattr(creatorObj,key) )

		#Check specific options as expected
		for expFileName, creatorObj in it.zip_longest(fileNames,outCreatorObjs):
			self.assertEqual( expFileName,creatorObj.fileName ) 

		for expGeom, creatorObj in it.zip_longest(self.gasPhaseReactantGeoms, outCreatorObjs):
			self.assertEqual( expGeom, creatorObj.geom )

	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.CodeSpecificStandardInputCreatorTemplate._getGasPhaseReacantObjs")
	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.CodeSpecificStandardInputCreatorTemplate._getBulkWithoutAdsorbedGeomCreator")
	def testGetReactantCalcObjs(self, mockedBulkCreator, mockedGasPhaseReactantsCreator):
		bulkCreatorObj = mock.Mock()
		bulkCalcObj = mock.Mock()
		gasPhaseCreatorObj = mock.Mock()
		gasPhaseCalcObj = mock.Mock()

		mockedBulkCreator.side_effect = lambda *args,**kwargs: bulkCreatorObj
		bulkCreatorObj.create.side_effect = lambda *args,**kwargs: bulkCalcObj
	
		mockedGasPhaseReactantsCreator.side_effect = lambda *args,**kwargs: [gasPhaseCreatorObj]
		gasPhaseCreatorObj.create.side_effect = lambda *args, **kwargs: gasPhaseCalcObj

		expOutput = [bulkCalcObj, gasPhaseCalcObj]
		actOutput = self.testObjA._getReactantCalcObjs()

		bulkCreatorObj.create.assert_called_once_with()
		gasPhaseCreatorObj.create.assert_called_once_with()
		self.assertEqual(expOutput,actOutput)


	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.CodeSpecificStandardInputCreatorTemplate._getGasPhaseProductObjs")
	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.CodeSpecificStandardInputCreatorTemplate._getBulkWithAdsorbedGeomCreator")
	def testGetProductCalcObjs(self, mockedBulkCreator, mockedGasPhaseProductsCreator):
		bulkCreatorObj = mock.Mock()
		bulkCalcObj = mock.Mock()
		gasPhaseCreatorObj = mock.Mock()
		gasPhaseCalcObj = mock.Mock()

		mockedBulkCreator.side_effect = lambda *args,**kwargs: bulkCreatorObj
		bulkCreatorObj.create.side_effect = lambda *args,**kwargs: bulkCalcObj
	
		mockedGasPhaseProductsCreator.side_effect = lambda *args,**kwargs: [gasPhaseCreatorObj]
		gasPhaseCreatorObj.create.side_effect = lambda *args, **kwargs: gasPhaseCalcObj

		expOutput = [bulkCalcObj, gasPhaseCalcObj]
		actOutput = self.testObjA._getProductCalcObjs()

		bulkCreatorObj.create.assert_called_once_with()
		gasPhaseCreatorObj.create.assert_called_once_with()
		self.assertEqual(expOutput,actOutput)


	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.totEnergyFlow.TotalEnergyWorkflow")
	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.CodeSpecificStandardInputCreatorTemplate._getReactantCalcObjs")
	def testCreateReactantTotalEnergyWorkflows(self, mockedCalcObjGetter, mockedTotEnergyFlow):
		reactantObjA = mock.Mock()
		reactantObjB = mock.Mock()
		reactantObjs = [reactantObjA, reactantObjB]
		outFlows = [mock.Mock(),mock.Mock()]

		mockedCalcObjGetter.side_effect = lambda *args,**kwargs: reactantObjs
		mockedTotEnergyFlow.side_effect = outFlows

		outObjs = self.testObjA._getReactantTotalEnergyWorkflows()

		mockedTotEnergyFlow.assert_any_call(reactantObjA)
		mockedTotEnergyFlow.assert_any_call(reactantObjB)

		self.assertEqual(outFlows, outObjs)

	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.totEnergyFlow.TotalEnergyGroupWorkflow")
	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.CodeSpecificStandardInputCreatorTemplate._getReactantTotalEnergyWorkflows")
	def testGetTotalReactantsWorkflow(self, reactantFlowGetters, mockedTotEnergyGroupFlow):
		bulkTotalEnergyFlow = mock.Mock()
		reactantAFlow = mock.Mock()
		reactantBFlow = mock.Mock()
		allReactantFlows = [bulkTotalEnergyFlow, reactantAFlow, reactantBFlow]
		allReactantStoics = [1] + self.gasPhaseReactantStoichiometries #Bulk is always 1st and walsy stoic=1
		expOutflow = mock.Mock()

		mockedTotEnergyGroupFlow.side_effect = expOutflow
		reactantFlowGetters.side_effect = lambda *args,**kwargs : [bulkTotalEnergyFlow, reactantAFlow, reactantBFlow]
		
		actOutflow = self.testObjA._getReactantsWorkflow()
		mockedTotEnergyGroupFlow.assert_called_once_with(allReactantFlows, allReactantStoics)
		self.assertTrue(expOutflow,actOutflow)

	def testGasPhaseProductStoicsReturnEmptyListWhenNotSet(self):
		self.gasPhaseProductStoichiometries = None
		self.createTestObjs()
		expProductStoics = list()
		actProductStoics = self.testObjA.gasPhaseProductStoichiometries
		self.assertEqual(expProductStoics, actProductStoics)

	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.reactFlow.ReactionEnergyWorkflow")
	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.CodeSpecificStandardInputCreatorTemplate._getProductsWorkflow")	
	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.CodeSpecificStandardInputCreatorTemplate._getReactantsWorkflow")
	def testCreateReactionWorkflow(self, mockedReactantWorkflowGetter, mockedProductWorkflowGetter, mockedReactionEnergyFlow):
		reactantsWorkflow = mock.Mock()
		productsWorkflow = mock.Mock()
		expOutflow = mock.Mock()

		mockedReactantWorkflowGetter.side_effect = [reactantsWorkflow]
		mockedProductWorkflowGetter.side_effect = [productsWorkflow]
		mockedReactionEnergyFlow.side_effect = [expOutflow]

		actOutflow = self.testObjA._getReactionWorkflow()
		mockedReactionEnergyFlow.assert_called_once_with(reactantsWorkflow, productsWorkflow)
		self.assertEqual(expOutflow, actOutflow)

	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.calcRunners.StandardInputObj")
	@mock.patch("gen_basis_helpers.job_helpers.solid_adsorption_energy.CodeSpecificStandardInputCreatorTemplate._getReactionWorkflow")
	def testCreateFromSelf(self, mockedWorkFlowGetter, mockedStdInp):
		expReactionFlow = mock.Mock()
		expOutObj = mock.Mock()
		expLabel = labelHelp.StandardLabel(eleKey=self.eleKey,structKey=self.structKey,methodKey=self.methodKey)

		mockedWorkFlowGetter.side_effect = [expReactionFlow]
		mockedStdInp.side_effect = [expOutObj]
		actOutObj = self.testObjA.create()

		mockedStdInp.assert_called_once_with(expReactionFlow, expLabel)
		self.assertEqual(expOutObj, actOutObj)




