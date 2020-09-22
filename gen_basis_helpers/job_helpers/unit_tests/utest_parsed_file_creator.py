
import copy
import itertools as it
import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.parsed_out_files as tCode

class TestParsedFileStdInpCreator(unittest.TestCase):

	def setUp(self):
		self.geomA = mock.Mock()
		self.geomB = mock.Mock()
		self.geoms = [self.geomA, self.geomB]
		self.baseCreatorA = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		kwargDictA = {"geoms":self.geoms, "baseCreator":self.baseCreatorA}
		self.testObjA = tCode.ParsedFileObjsForMultiGeomsStandardInputCreator( **kwargDictA )

	def testGetFileNames(self):
		expFileNames = ["inp_file_1", "inp_file_2"]
		actFileNames = self.testObjA._getFileNamesForAll()
		self.assertEqual(expFileNames, actFileNames)

	@mock.patch("gen_basis_helpers.job_helpers.parsed_out_files.ParsedFileObjsForMultiGeomsStandardInputCreator.outFolder", new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_helpers.parsed_out_files.ParsedFileObjsForMultiGeomsStandardInputCreator._getFileNamesForAll")
	def testGetCreators(self, mockFileNames, mockedOutFolderProp):
		expFileNames = ["nameA", "nameB"]
		expFolder = "fake_folder"
		mockFileNames.side_effect = lambda : expFileNames
		mockedOutFolderProp.return_value = expFolder

		expCreators = [copy.deepcopy(self.baseCreatorA) for x in self.geoms]
		expCreators[0].fileName, expCreators[1].fileName = expFileNames
		expCreators[0].geom, expCreators[1].geom = self.geoms

		actCreators = self.testObjA._getAllCreators()

		for exp,act in it.zip_longest(expCreators, actCreators):
			self.assertEqual(expFolder, act.workFolder)
			self.assertEqual(exp.geom,act.geom)
			self.assertEqual(exp.fileName, act.fileName)

	@mock.patch("gen_basis_helpers.job_helpers.parsed_out_files.ParsedFileObjsForMultiGeomsStandardInputCreator._getAllCreators")
	def testGetCalcObjs(self, mockedGetCreators):
		expCreators = [mock.Mock(), mock.Mock()]
		expCalcObjs = [mock.Mock() for x in expCreators]
		expCreators[0].create.side_effect = lambda :expCalcObjs[0]
		expCreators[1].create.side_effect = lambda :expCalcObjs[1]

		mockedGetCreators.side_effect  = lambda : expCreators
		actCalcObjs = self.testObjA._getAllCalcObjs()
		self.assertEqual(expCalcObjs,actCalcObjs)

	@mock.patch("gen_basis_helpers.job_helpers.parsed_out_files.baseFlow.StandardLabelledWorkflowComposite")
	@mock.patch("gen_basis_helpers.job_helpers.parsed_out_files.parsedFileFlow.ParsedFileWorkflow")
	@mock.patch("gen_basis_helpers.job_helpers.parsed_out_files.ParsedFileObjsForMultiGeomsStandardInputCreator._getAllCalcObjs")
	def testGetWorkflow(self, mockGetCalcObjs, mockSingleWorkflow, mockedCompositeWorkflow):
		expCalcObjs = [mock.Mock(), mock.Mock()]
		expOutput = mock.Mock()
		mockGetCalcObjs.side_effect = lambda : expCalcObjs
		mockSingleWorkflow.side_effect = lambda x, **kwargs: x
		mockedCompositeWorkflow.side_effect = lambda *args: expOutput

		actOutput = self.testObjA._getWorkflow()

		for x in expCalcObjs:
			mockSingleWorkflow.assert_any_call(x,catchParserErrors=False)
		mockedCompositeWorkflow.assert_called_with(expCalcObjs) #Since single mocked workflow just returns the mock calcObjs
		self.assertEqual(expOutput, actOutput)

	@mock.patch("gen_basis_helpers.job_helpers.parsed_out_files.calcRunners.StandardInputObj")
	@mock.patch("gen_basis_helpers.job_helpers.parsed_out_files.ParsedFileObjsForMultiGeomsStandardInputCreator._getWorkflow")
	@mock.patch("gen_basis_helpers.job_helpers.parsed_out_files.ParsedFileObjsForMultiGeomsStandardInputCreator.label", new_callable=mock.PropertyMock)
	def testCreateFromSelf(self, mockLabel, mockGetWorkflow, mockStdInpObj):
		expLabel, expWorkflow = mock.Mock(), mock.Mock()
		expOutput = mock.Mock()

		mockLabel.return_value = expLabel
		mockGetWorkflow.side_effect = lambda : expWorkflow
		mockStdInpObj.side_effect = lambda *args: expOutput

		actOutput = self.testObjA.create()

		mockGetWorkflow.assert_called_with()
		mockStdInpObj.assert_called_with(expWorkflow, expLabel)
		self.assertEqual(expOutput, actOutput)


