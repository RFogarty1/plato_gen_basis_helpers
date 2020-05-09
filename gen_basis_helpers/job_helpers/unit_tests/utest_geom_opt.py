

import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.geom_opt as tCode

class TestCreateStandardInputForGeomOpt(unittest.TestCase):

	def setUp(self):
		self.baseCreator = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.CodeSpecificStandardInputCreatorTemplate( baseCreator=self.baseCreator )

	@mock.patch("gen_basis_helpers.job_helpers.geom_opt.CodeSpecificStandardInputCreatorTemplate.outFolder", new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_helpers.geom_opt.CodeSpecificStandardInputCreatorTemplate._getBaseCreator")
	def testGeomOptSetForCreator(self, mockedGetBaseCreator, mockedOutFolder):
		mockedGetBaseCreator.side_effect = [ self.baseCreator ]
		expRunType = "geomOpt"
		actCreator = self.testObjA._getCreatorObj()
		actRunType = actCreator.runType
		self.assertEqual(expRunType.lower(), actRunType.lower())

	@mock.patch("gen_basis_helpers.job_helpers.geom_opt.CodeSpecificStandardInputCreatorTemplate.outFolder", new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_helpers.geom_opt.CodeSpecificStandardInputCreatorTemplate._getBaseCreator")
	def testOutFolderSetForCreator(self, mockedBaseCreator, mockedOutFolderProp):	
		expOutFolder = "fake_folder_path"
		mockedOutFolderProp.return_value = expOutFolder
		actCreator = self.testObjA._getCreatorObj()
		actOutFolder = actCreator.workFolder
		self.assertEqual(expOutFolder, actOutFolder)


	@mock.patch("gen_basis_helpers.job_helpers.geom_opt.goptFlow.GeomOptWorkflow")
	@mock.patch("gen_basis_helpers.job_helpers.geom_opt.CodeSpecificStandardInputCreatorTemplate._getCreatorObj")
	def testExpWorkflowCreated(self, mockedGetCreator, mockedWorkflow):
		expWorkflow, expCreator, expCalcObj = mock.Mock(), mock.Mock(), mock.Mock()
		mockedGetCreator.side_effect = [expCreator]
		expCreator.create.side_effect = [expCalcObj]
		mockedWorkflow.side_effect = [expWorkflow]
		actWorkflow = self.testObjA._getWorkflow()
		mockedWorkflow.assert_called_once_with(expCalcObj)
		self.assertEqual(expWorkflow,actWorkflow)

	@mock.patch("gen_basis_helpers.job_helpers.geom_opt.CodeSpecificStandardInputCreatorTemplate.label", new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_helpers.geom_opt.calcRunners.StandardInputObj")
	@mock.patch("gen_basis_helpers.job_helpers.geom_opt.CodeSpecificStandardInputCreatorTemplate._getWorkflow")
	def testExpStandardInpCallMade(self, mockedWorkflowCreatorFunct, mockedStdInp, mockedLabel):
		expWorkflow, expOutput, expLabel = mock.Mock(), mock.Mock(), mock.Mock()
		mockedWorkflowCreatorFunct.side_effect = [expWorkflow]
		mockedStdInp.side_effect = [expOutput]
		mockedLabel.return_value = expLabel
		actOutput = self.testObjA.create()
		mockedStdInp.assert_called_once_with(expWorkflow, expLabel)
		self.assertEqual(expOutput, actOutput)



