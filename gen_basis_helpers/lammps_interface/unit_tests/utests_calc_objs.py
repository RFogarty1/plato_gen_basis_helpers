
import copy
import collections
import os
import unittest
import unittest.mock as mock

import gen_basis_helpers.lammps_interface.lammps_calc_objs as tCode


class TestCalcObjStandard(unittest.TestCase):

	def setUp(self):
		self.baseFolder = "fake_path"
		self.baseFileName = "test_file"
		self.dataFileOrderedDict = collections.OrderedDict([["fake_header","fake_val"]])
		self.scriptFileOrderedDict = collections.OrderedDict([["read_data","datafile"],["commA","valA"]])
		self.createTestObjs()

	def createTestObjs(self):
		argList = [self.baseFolder, self.baseFileName, self.dataFileOrderedDict, self.scriptFileOrderedDict]
		self.testObjA = tCode.LammpsCalcObjStandard(*argList)

	def testScriptFilePathAsExpected(self):
		expPath = os.path.join(self.baseFolder, self.baseFileName) + ".in"
		actPath = self.testObjA.scriptFilePath
		self.assertEqual(expPath,actPath) 

	def testDataFilePathAsExpected(self):
		expPath = os.path.join(self.baseFolder, self.baseFileName) + ".data"
		actPath = self.testObjA.dataFilePath
		self.assertEqual(expPath,actPath)


	@mock.patch("gen_basis_helpers.lammps_interface.lammps_calc_objs.fileIoHelp")
	@mock.patch("gen_basis_helpers.lammps_interface.lammps_calc_objs.pathlib")
	@mock.patch("gen_basis_helpers.lammps_interface.lammps_calc_objs.LammpsCalcObjStandard.dataFilePath", new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.lammps_interface.lammps_calc_objs.LammpsCalcObjStandard.scriptFilePath", new_callable=mock.PropertyMock)
	def testWriteFileCallsExpected(self, mockScriptPathProp, mockDataPathProp, mockedPathLib, mockedFileIo):
		expScriptPath, expDataPath = "script_path", "data_path"
		mockScriptPathProp.return_value = expScriptPath
		mockDataPathProp.return_value = expDataPath
		expScriptDict  = copy.deepcopy(self.scriptFileOrderedDict)
		expScriptDict["read_data"] = expDataPath #Really only needs the fileName rather than path
		self.testObjA.writeFile()
		mockedPathLib.Path.assert_called_with(self.baseFolder)
		mockedFileIo.writeScriptFileFromTokens.assert_called_with(expScriptPath, self.scriptFileOrderedDict)
		mockedFileIo.writeDataFileFromTokens.assert_called_with(expDataPath, self.dataFileOrderedDict)

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_calc_objs.LammpsCalcObjStandard.scriptFilePath", new_callable=mock.PropertyMock)
	def testExpectedRunComm(self, mockedScriptPath):
		expScriptPath = "some_dir/test_script_path"
		mockedScriptPath.return_value = expScriptPath
		expRunComm = "lmp -in test_script_path"
		actRunComm = self.testObjA.runComm
		self.assertEqual(expRunComm,actRunComm)


class TestScriptFileOptsStandard(unittest.TestCase):

	def setUp(self):
		self.initOpts = collections.OrderedDict([["initOptKey","initOptVal"]])
		self.setupBox = collections.OrderedDict([["setupBoxKey","setupBoxVal"]])
		self.setupAtoms = collections.OrderedDict([["setupAtomsKey","setupAtomsVal"]])
		self.forceFieldOpts = collections.OrderedDict([["forceFieldKey","forceFieldVal"]])
		self.settingsOpts = collections.OrderedDict([["settingsKey","settingsVal"]])
		self.fixOpts = collections.OrderedDict([["fixKey","fixVals"]])
		self.outputSection = collections.OrderedDict([["outputKey","outputVal"]])
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"initOpts": self.initOpts, "setupBox": self.setupBox,
		             "setupAtoms": self.setupAtoms, "forceFieldOpts": self.forceFieldOpts,
		             "settingsOpts": self.settingsOpts, "fixSection": self.fixOpts,
		             "outputSection": self.outputSection}

		self.testObjA = tCode.ScriptFileOptionsStandard(**kwargDict)

	def testExpectedDict_allOptsSet(self):
		expDict = self._loadExpectedDictA()
		actDict = self.testObjA.getOutputDict()
		self.assertEqual(expDict,actDict)

	def testExpectedDict_fixOptsNotSet(self):
		self.fixOpts = None
		self.createTestObjs()
		expDict = self._loadExpectedDictA()
		expDict.pop("fixKey")
		actDict = self.testObjA.getOutputDict()
		self.assertEqual(expDict, actDict)

	def _loadExpectedDictA(self):
		expDict = collections.OrderedDict()
		expDict["initOptKey"] = "initOptVal"
		expDict["setupBoxKey"] = "setupBoxVal"
		expDict["setupAtomsKey"] = "setupAtomsVal"
		expDict["forceFieldKey"] = "forceFieldVal"
		expDict["settingsKey"] = "settingsVal"
		expDict["fixKey"] = "fixVals"
		expDict["outputKey"] = "outputVal"
		return expDict

