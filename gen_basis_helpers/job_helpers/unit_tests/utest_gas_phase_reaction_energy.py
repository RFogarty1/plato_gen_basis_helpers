
import itertools as it
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.gas_phase_reaction_energy as tCode

import gen_basis_helpers.shared.calc_runners as calcRunners
import gen_basis_helpers.shared.label_objs as labelHelp

class TestCreateStandardInputForGasPhaseReactionEnergy(unittest.TestCase):


	def setUp(self):
		self.baseWorkFolder = "fake_folder"
		self.reactantGeoms = [mock.Mock(), mock.Mock()]
		self.reactantStoichiometries = [x for x in range(len(self.reactantGeoms))]
		self.reactantCharges = [2 for x in range(len(self.reactantGeoms))]
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"baseWorkFolder":self.baseWorkFolder,
		             "reactantGeoms":self.reactantGeoms,
		              "reactantStoichiometries":self.reactantStoichiometries,
		             "reactantCharges":self.reactantCharges}
		self.testObjA = tCode.CodeSpecificStandardInputCreatorTemplate(**kwargDict)

	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.CodeSpecificStandardInputCreatorTemplate.outFolder", new_callable=mock.PropertyMock)
	def testExpectedSharedModificationsApplied(self, mockedOutFolder):
		mockCreator = mock.Mock()
		expFolder = "fake_folder_path"
		mockedOutFolder.return_value = expFolder
		expKPts = [1,1,1]
		self.testObjA._applySharedOptionsToCreator(mockCreator)
		self.assertEqual(expKPts,mockCreator.kPts)
		self.assertEqual(expFolder, mockCreator.workFolder)

	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.CodeSpecificStandardInputCreatorTemplate._applySharedOptionsToCreator")
	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.CodeSpecificStandardInputCreatorTemplate._getBaseCreator")
	def testGetCreatorsForSetOfGeoms(self, mockedGetBaseCreator, mockedApplySharedOptions):
		mockedBaseCreators = [mock.Mock() for x in self.reactantGeoms]
		mockedGetBaseCreator.side_effect = mockedBaseCreators #Each call gets a different mock

		outCreators = self.testObjA._getCreatorsForSetOfGeoms(self.reactantGeoms)

		self.assertEqual( mockedBaseCreators, outCreators)

		for creator,unused in it.zip_longest(outCreators,self.reactantGeoms):
			mockedApplySharedOptions.assert_any_call(creator)

	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.CodeSpecificStandardInputCreatorTemplate._getCreatorsForSetOfGeoms")
	def testExpectedReactantCreators(self, mockedCreatorsForGeomSet):
		expCreators = [mock.Mock() for x in range(len(self.reactantGeoms))]
		expFileNames = ["reactant_{}".format(x) for x in range(len(self.reactantGeoms))]
		expCharges = self.reactantCharges
		mockedCreatorsForGeomSet.side_effect = lambda *args,**kwargs: expCreators

		outCalcObjs = self.testObjA._getReactantCalcObjs()
		actFileNames = [x.fileName for x in expCreators]
		actCharges = [x.charge for x in expCreators]

		mockedCreatorsForGeomSet.assert_called_once_with(self.reactantGeoms)
		for x in expCreators:
			x.create.assert_called_once_with()

		self.assertEqual(expFileNames,actFileNames)
		self.assertEqual(expCharges,actCharges)

	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.CodeSpecificStandardInputCreatorTemplate._getReactantCalcObjs")
	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.totEnergyFlow.TotalEnergyGroupWorkflow")
	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.totEnergyFlow.TotalEnergyWorkflow")
	def testExpectedReactantWorkflow(self, mockTotEnergyFlow, mockGroupedEnergyFlow,mockedGetCalcObjs):
		expCalcObjs = [mock.Mock() for x in self.reactantGeoms]
		expTotEnergyFlows = [mock.Mock() for x in self.reactantGeoms]
		expOutObj = mock.Mock()

		mockedGetCalcObjs.side_effect = lambda *args,**kwargs:expCalcObjs
		mockTotEnergyFlow.side_effect = expTotEnergyFlows
		mockGroupedEnergyFlow.side_effect = lambda *args,**kwargs:expOutObj

		outFlow = self.testObjA._getReactantsTotalEnergyWorkflow()

		for totEWorkflow in expCalcObjs:
			mockTotEnergyFlow.assert_any_call(totEWorkflow)
		mockGroupedEnergyFlow.assert_called_once_with(expTotEnergyFlows, self.reactantStoichiometries)
		self.assertEqual(expOutObj,outFlow)


	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.reactFlow.ReactionEnergyWorkflow")
	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.CodeSpecificStandardInputCreatorTemplate._getProductsTotalEnergyWorkflow")
	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.CodeSpecificStandardInputCreatorTemplate._getReactantsTotalEnergyWorkflow")
	def testGetReactionEnergyWorkflow(self, mockedReactantsGetter, mockedProductsGetter, mockedReactionEFlow):
		expReactants = mock.Mock()
		expProducts = mock.Mock()
		expOutput = mock.Mock()

		mockedReactantsGetter.side_effect = [expReactants]
		mockedProductsGetter.side_effect = [expProducts]
		mockedReactionEFlow.side_effect = [expOutput]

		actOutput = self.testObjA._getReactionEnergyWorkflow()
		mockedReactionEFlow.assert_called_once_with(expReactants, expProducts)
		self.assertEqual(expOutput,actOutput)

	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.calcRunners.StandardInputObj")
	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.CodeSpecificStandardInputCreatorTemplate.label", new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_helpers.gas_phase_reaction_energy.CodeSpecificStandardInputCreatorTemplate._getReactionEnergyWorkflow")
	def testCreateStandardInput(self, mockedReactionFlow, mockedLabel, mockedStdInp):
		expReactionFlow = mock.Mock()
		expLabel = mock.Mock()
		expOutput = mock.Mock()

		mockedReactionFlow.side_effect = [expReactionFlow]
		mockedLabel.return_value = expLabel
		mockedStdInp.side_effect = [expOutput]

		actOutput = self.testObjA.create()
		mockedStdInp.assert_called_once_with( expReactionFlow, expLabel )
		self.assertEqual( expOutput, actOutput )	





class TestMapFunction(unittest.TestCase):

	def setUp(self):
		self.eleKey = "fake_ele"
		self.structKey = "fake_struct"
		self.methodKey = "fake_method"
		self.reactEnergyA = 40
		self.reactEnergyFmt = "{:4f}"
		self.createTestObjs()

	def createTestObjs(self):
		labelA = labelHelp.StandardLabel( eleKey=self.eleKey, methodKey=self.methodKey, structKey=self.structKey )
		outputObj = types.SimpleNamespace( energy=self.reactEnergyA )

		self.workflowA = types.SimpleNamespace( output=[outputObj], run=mock.Mock() )
		self.testObjA = tCode.MapGasPhaseReactionEnergyWorkflowToUsefulFormatStandard(reactEnergyFmt=self.reactEnergyFmt)
		self.standardInpObjA = calcRunners.StandardInputObj( self.workflowA, labelA )

	def _runTestFunct(self):
		return self.testObjA( self.standardInpObjA )

	def testRunMethodCalled(self):
		self._runTestFunct()
		self.workflowA.run.assert_called_once_with()

	def testExpectedSimpleTableData(self):
		self._runTestFunct()
		expData = [self.methodKey] + [self.reactEnergyFmt.format(self.reactEnergyA)]
		actData = self._runTestFunct().tableData
		self.assertEqual(expData,actData)







