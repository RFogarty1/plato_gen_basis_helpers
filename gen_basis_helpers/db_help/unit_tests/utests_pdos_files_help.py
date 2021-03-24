
import os
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.calc_runners as calcRunnersHelp
import gen_basis_helpers.shared.label_objs as labelHelp

import gen_basis_helpers.cp2k.parse_pdos_files as parsePdosHelp
import gen_basis_helpers.db_help.pdos_files_help as tCode


class TestGetPdosFromRecord(unittest.TestCase):

	def setUp(self):
		self.pdosFileName = "fake.json"
		self.pdosPath = "fake_dir"
		self.startDir = "start_dir"
		self.dbExtPath = "ext_a"
		self.atomKindDictA = {"eigenValues":[2,3]}
		self.atomListDictA = {"eigenValues":[4,5,6]}
		self.createTestObjs()

	def createTestObjs(self):
		self.pdosDictFromFile = {"atomKinds": [self.atomKindDictA] , "atomLists": [self.atomListDictA]}
		self.expPdosDict = {"atomKinds": [parsePdosHelp.PdosFragmentStandard.fromDict(self.atomKindDictA)],
		                    "atomLists": [parsePdosHelp.PdosFragmentStandard.fromDict(self.atomListDictA)]}
		self.recordA = {"pdos_path_ext":self.pdosPath, "pdos_filename":self.pdosFileName, "db_ext_path":self.dbExtPath}


	@mock.patch("gen_basis_helpers.db_help.pdos_files_help._readSinglePdosIntoListFromJsonFile")
	def testExpectedFilePath(self, mockedReadFromJson):
		mockedReadFromJson.side_effect = lambda *args,**kwargs: self.pdosDictFromFile
		expPath = os.path.join(self.startDir, self.dbExtPath, self.pdosPath, self.pdosFileName)
		expDict = self.expPdosDict
		actDict = tCode.getPdosDictFromRecord(self.startDir, self.recordA)
		mockedReadFromJson.assert_called_with(expPath)
		self.assertEqual(expDict, actDict)


class TestDumpPdosFiles(unittest.TestCase):

	def setUp(self):
		self.eleKey = "ele_key"
		self.structKey = "struct_key"
		self.methodKey = "method_key"
		self.startDir = "fake_dir"

		self.atomPdos = [mock.Mock(), mock.Mock()]
		self.listPdos = [mock.Mock()]

		self.atomDicts = [{"key_a":"val_a"}, {"key_b":"val_b"}]
		self.listDicts = [{"key_c":"val_c"}]
		self.createTestObjs()

	def createTestObjs(self):
		self.labelA = labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)

		self.dataDictA = {"pdos":{"atomKinds":self.atomPdos, "atomLists":self.listPdos}}

		self.parsedFileA = types.SimpleNamespace(**self.dataDictA)
		dataObj = types.SimpleNamespace( parsedFile=self.parsedFileA )
		self.stdOutObjA = calcRunnersHelp.StandardOutputObj([dataObj], self.labelA)


	@mock.patch("gen_basis_helpers.db_help.pdos_files_help.dictsDbHelp.dumpDictsToFilePath")
	@mock.patch("gen_basis_helpers.db_help.pdos_files_help._turnListOfPdosObjsToListOfDicts")
	def testExpectedCallsMade(self, mockGetDictList, mockDumper):
		#Setup mock function
		def _getListOfDictsStub(inpList):
			if inpList==self.atomPdos:
				return self.atomDicts
			elif inpList==self.listPdos:
				return self.listDicts
			else:
				raise ValueError("")

		mockGetDictList.side_effect = _getListOfDictsStub
		tCode.dumpPdosFilesFromStdOut_simple(self.stdOutObjA, self.startDir)
		expDict = {"atomKinds":self.atomDicts, "atomLists":self.listDicts}
		expFilePath = os.path.join(self.startDir,self.eleKey, self.structKey, self.methodKey, "pdos.json")
		mockDumper.assert_called_with( [expDict], expFilePath) 





