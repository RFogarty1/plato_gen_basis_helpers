
import os
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.calc_runners as calcRunnersHelp
import gen_basis_helpers.shared.label_objs as labelHelp

import gen_basis_helpers.db_help.nudged_band_files_help as tCode


class TestReadInPathwayFromStandardRecord(unittest.TestCase):

	def setUp(self):
		self.startDir = "fake_start_dir"
		self.extPath = "fake_ext_path"
		self.fileName = "fake_file_name"
		self.createTestObjs()

	def createTestObjs(self):
		self.recordA = {"neb_path_ext":self.extPath, "neb_pathway_filename":self.fileName}

	@mock.patch("gen_basis_helpers.db_help.nudged_band_files_help.pathwayHelp.readNudgedBandPathwayFromJsonFileStandard")
	def testExpectedCallsMade(self, mockedReadInpFile):
		expOutObj = mock.Mock()
		expPath = os.path.join(self.startDir, self.extPath, self.fileName)
		mockedReadInpFile.side_effect = lambda *args,**kwargs: expOutObj
		actOutObj = tCode.getNebPathwayFromRecord(self.startDir, self.recordA)
		mockedReadInpFile.assert_called_with(expPath)
		self.assertEqual(expOutObj,actOutObj)


class TestDumpNebFiles(unittest.TestCase):

	def setUp(self):
		self.eleKey = "ele_key"
		self.structKey = "struct_key"
		self.methodKey = "method_key"
		self.startDir = "fake_dir"
		self.createTestObjs()

	def createTestObjs(self):
		self.labelA = labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)

		self.dataDictA = {"fake_key":"fake_val"}

		self.parsedFileA = types.SimpleNamespace(**self.dataDictA)
		dataObj = types.SimpleNamespace( parsedFile=self.parsedFileA )
		self.stdOutObjA = calcRunnersHelp.StandardOutputObj([dataObj], self.labelA)

	@mock.patch("gen_basis_helpers.db_help.nudged_band_files_help.parseNebHelp.getNebPathFromParsedFileObj")
	@mock.patch("gen_basis_helpers.db_help.nudged_band_files_help.sharedIoHelp.dumpObjWithToDictToJson")
	def testExpectedCallsToDump(self, mockDump, mockGetPathway):
		expNebPathway = mock.Mock()
		expPath = os.path.join(self.startDir, self.eleKey, self.structKey, self.methodKey, "neb_path.json")
		mockGetPathway.side_effect = lambda *args,**kwargs: expNebPathway

		tCode.dumpNebFilesFromStdOut_simple(self.stdOutObjA, self.startDir)
		mockGetPathway.assert_called_with(self.parsedFileA)
		mockDump.assert_called_with(expNebPathway, expPath)


