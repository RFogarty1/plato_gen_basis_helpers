
import os
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.analyse_md.analyse_metadyn_hills as aMetadynHillsHelp
import gen_basis_helpers.shared.calc_runners as calcRunnersHelp
import gen_basis_helpers.shared.label_objs as labelHelp

import gen_basis_helpers.db_help.md_files_help as tCode


class TestGetObjsFromBaseFolderAndRecord(unittest.TestCase):

	def setUp(self):
		self.baseDbFolder = "base_db_folder"
		self.mdPathExt = "fake_md_path_ext"
		self.trajFileName = "fake_traj_name"
		self.thermoFileName = "fake_thermo_name"
		self.createTestObjs()

	def createTestObjs(self):
		self.recordA = {"md_path_ext": self.mdPathExt, "md_thermo_filename": self.thermoFileName,
		                "md_traj_filename": self.trajFileName}

	@mock.patch("gen_basis_helpers.db_help.md_files_help.thermoDataHelp.readThermoDataFromFile")
	@mock.patch("gen_basis_helpers.db_help.md_files_help.getMdFilePathFromRecord")
	def testGetThermoObj(self, mockedGetPath, mockedReadThermoFile):
		expPath, expObj = "fake_out_path", mock.Mock()
		mockedGetPath.side_effect = lambda *args,**kwargs: expPath
		mockedReadThermoFile.side_effect = lambda *args,**kwargs: expObj

		actObj = tCode.getMdThermoObjFromBaseDbFolderAndRecord(self.baseDbFolder, self.recordA)
		mockedGetPath.assert_called_with(self.baseDbFolder, self.recordA, "md_thermo_filename")
		mockedReadThermoFile.assert_called_with(expPath)
		self.assertEqual(expObj, actObj)

	@mock.patch("gen_basis_helpers.db_help.md_files_help.trajHelp.readTrajObjFromFileToTrajectoryInMemory")
	@mock.patch("gen_basis_helpers.db_help.md_files_help.getMdFilePathFromRecord")
	def testGetTrajObj(self, mockedGetPath, mockedReadTrajFile):
		expPath, expObj = "fake_out_path", mock.Mock()
		mockedGetPath.side_effect = lambda *args,**kwargs: expPath
		mockedReadTrajFile.side_effect=  lambda *args,**kwargs: expObj

		actObj = tCode.getMdTrajInMemFromBaseDbFolderAndRecord(self.baseDbFolder, self.recordA)
		mockedGetPath.assert_called_with(self.baseDbFolder, self.recordA, "md_traj_filename")
		mockedReadTrajFile.assert_called_with(expPath)
		self.assertEqual(expObj,actObj)


class TestGetPathsFromMdRecordStandard(unittest.TestCase):

	def setUp(self):
		self.baseDbFolder = "fake_base_db_folder"
		self.dbExtPath = "fake_db_ext_path"
		self.mdPathExt = "fake_md_path_ext"
		self.mdThermoFileName = "fake_thermo_name"
		self.trajFileName = "fake_traj_name"
		self.restartFileName = "fake_restart_file"
		self.wfnFileName = "fake_wfn_file"
		self.createTestObjs()

	def createTestObjs(self):
		self.recordA = {"md_path_ext": self.mdPathExt, "md_thermo_filename": self.mdThermoFileName,
		                "md_traj_filename": self.trajFileName, "db_ext_path": self.dbExtPath,
		                "wfn_filename":self.wfnFileName, "restart_filename": self.restartFileName}

	def testExpPathForTrajFile(self):
		expPath = os.path.join(self.baseDbFolder, self.dbExtPath, self.mdPathExt, self.trajFileName)
		actPath = tCode.getMdFilePathFromRecord(self.baseDbFolder,self.recordA, "md_traj_filename")
		self.assertEqual(expPath,actPath)

	def testExpPathWfn(self):
		expPath = os.path.join(self.baseDbFolder, self.dbExtPath, self.mdPathExt, self.wfnFileName)
		actPath = tCode.getMdFilePathFromRecord(self.baseDbFolder, self.recordA, "wfn_filename")
		self.assertEqual(expPath,actPath)


class TestGetOutDictForMdSimple(unittest.TestCase):

	def setUp(self):
		self.stdOutObj = mock.Mock()
		self.writeFiles = True
		self.startDir = "fake_dir"
		self.maxWfnBackups = 6

	def _runTestFunct(self):
		return tCode.getOutDictForMDFromStdOutObj_simple(self.startDir, self.stdOutObj,
		                                                 writeFiles=self.writeFiles, maxNumbWfnBackups=self.maxWfnBackups)

	@mock.patch("gen_basis_helpers.db_help.md_files_help.MDFilesFromOutObjsFileDumperStandard")
	def testExpectedCallsA_writeFile(self, mockedDumpClass):
		expInstance = mock.Mock()
		mockedDumpClass.side_effect = lambda *args, **kwargs: expInstance
		self._runTestFunct()
		mockedDumpClass.assert_called_with(copyWfnFile=True, copyRestartFile=True, maxNumbWfnBackups=self.maxWfnBackups)
		expInstance.dumpFiles.assert_called_with(self.stdOutObj, self.startDir)


class TestGetStdMdFilesOutputDictFromStdInpCreatorAndOutput(unittest.TestCase):

	def setUp(self):
		self.runFolder = "fake_folder_path"
		self.thermoObj = mock.Mock()
		self.trajObj = mock.Mock()
		#Using these really means we're restricted to standard inp creators which use the StandardLabel class; though not really an issue
		self.eleKey = "ele_key"
		self.structKey = "struct_key"
		self.methodKey = "method_key"

		self.copyWfnFiles = True
		self.copyRestartFiles = True
		self.maxNumbWfnBackups = 1

		self.createTestObjs()

	def createTestObjs(self):
		self.workflowA = mock.Mock()
		self.labelA = labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)
		self.stdInpObjA = calcRunnersHelp.StandardInputObj(self.workflowA, self.labelA)

		self.dataDictA = {"thermo_data":self.thermoObj, "trajectory":self.trajObj, "finalRunFolder":self.runFolder}
		self.dataA = [types.SimpleNamespace(parsedFile=types.SimpleNamespace(**self.dataDictA))]
		self.stdOutObjA = calcRunnersHelp.StandardOutputObj(self.dataA, self.labelA)

		currKwargs = {"copyWfnFile":self.copyWfnFiles, "copyRestartFile":self.copyRestartFiles, "maxNumbWfnBackups":self.maxNumbWfnBackups}
		self.testObjA = tCode.MDFilesFromOutObjsFileDumperStandard(**currKwargs)
		self.expPathExtA = os.path.join(self.eleKey, self.structKey, self.methodKey)


	def _loadExpPathExtDictA(self):
		return {"md_path_ext": self.expPathExtA,
		           "md_traj_filename": "out_traj.traj",
		           "md_thermo_filename":"out_thermo.thermo",
		           "restart_filename":"inp_file_1-1.restart",
	               "wfn_filename":"inp_file_1-RESTART.wfn",
		           "final_run_dir": self.runFolder}


	@mock.patch("gen_basis_helpers.db_help.md_files_help.os.listdir")
	def testGetPathExtDict(self, mockedListDir):
		mockDirContents = ["inp_file_1-1.restart", "inp_file_1-2.restart", "inp_file_1-1.restart.bak-1",
		                 "inp_file_1-RESTART.wfn",  "inp_file_1-RESTART.wfn.bak-1"]
		mockedListDir.side_effect = lambda *args,**kwargs: mockDirContents
		expDict = self._loadExpPathExtDictA()
		actDict = self.testObjA._getFilePathExtensionDict(self.stdOutObjA)
		self.assertEqual(expDict,actDict)

	@mock.patch("gen_basis_helpers.db_help.md_files_help.pathlib")
	@mock.patch("gen_basis_helpers.db_help.md_files_help.shutil.copy2")
	def testExpCallsTCopyRestartFilesMade(self, mockCopy2, mockPathlib):
		startDir = "fake_dir"
		pathExtDict = self._loadExpPathExtDictA()
		self.testObjA._copyRestartFiles(startDir, pathExtDict)


		mockCopy2.assert_any_call( os.path.join(self.runFolder,  "inp_file_1-1.restart") , os.path.join(startDir, self.expPathExtA, "inp_file_1-1.restart") )
		mockCopy2.assert_any_call( os.path.join(self.runFolder, "inp_file_1-RESTART.wfn"), os.path.join(startDir, self.expPathExtA, "inp_file_1-RESTART.wfn") )

	@mock.patch("gen_basis_helpers.db_help.md_files_help.pathlib")
	@mock.patch("gen_basis_helpers.db_help.md_files_help.trajHelp.dumpTrajObjToFile")
	@mock.patch("gen_basis_helpers.db_help.md_files_help.thermoDataHelp.dumpStandardThermoDataToFile")
	def testDumpInfoFilesCallsMade(self, dumpThermoMock, dumpTrajMock, mockPathlib):
		startDir = "fake_dir"
		expThermoPath, expTrajPath = os.path.join(startDir, self.expPathExtA, "out_thermo.thermo"), os.path.join(startDir, self.expPathExtA, "out_traj.traj")
		self.testObjA._dumpMdFiles(self.stdOutObjA, startDir, self._loadExpPathExtDictA())
		dumpThermoMock.assert_called_with(self.thermoObj, expThermoPath)
		dumpTrajMock.assert_called_with(self.trajObj, expTrajPath)

	@mock.patch("gen_basis_helpers.db_help.md_files_help.MDFilesFromOutObjsFileDumperStandard._dumpMdFiles")
	@mock.patch("gen_basis_helpers.db_help.md_files_help.MDFilesFromOutObjsFileDumperStandard._copyRestartFiles")
	@mock.patch("gen_basis_helpers.db_help.md_files_help.MDFilesFromOutObjsFileDumperStandard._getFilePathExtensionDict")
	def testExpectedCallsForDumpFiles(self, mockGetExtDict, mockCopyRestart, mockDumpMdFiles):
		#Set mock function return values
		expPathExtDict, expStartDir, expStdOut = {"key":"val"}, "fake-dir", mock.Mock()
		mockGetExtDict.side_effect = lambda *args,**kwargs: expPathExtDict

		actPathExtDict = self.testObjA.dumpFiles(expStdOut, expStartDir)
		mockGetExtDict.assert_called_with(expStdOut)
		mockCopyRestart.assert_called_with(expStartDir, expPathExtDict)
		mockDumpMdFiles.assert_called_with(expStdOut, expStartDir, expPathExtDict)

		self.assertEqual(expPathExtDict, actPathExtDict)


	@mock.patch("gen_basis_helpers.db_help.md_files_help.MDFilesFromOutObjsFileDumperStandard._dumpMdFiles")
	@mock.patch("gen_basis_helpers.db_help.md_files_help.MDFilesFromOutObjsFileDumperStandard._copyRestartFiles")
	@mock.patch("gen_basis_helpers.db_help.md_files_help.MDFilesFromOutObjsFileDumperStandard._getFilePathExtensionDict")
	def testExpectedCallsForDumpFilesWhenWriteFilesOff(self, mockGetExtDict, mockCopyRestart, mockDumpMdFiles):
		self.copyWfnFiles = False
		self.copyRestartFiles = False
		self.createTestObjs()

		#Set mock function return values
		expPathExtDict, expStdOut = {"key":"val"}, mock.Mock()
		mockGetExtDict.side_effect = lambda *args,**kwargs: expPathExtDict

		actPathExtDict = self.testObjA.dumpFiles(expStdOut, None)
		mockGetExtDict.assert_called_with(expStdOut)
		mockCopyRestart.assert_not_called()
		mockDumpMdFiles.assert_not_called()

		self.assertEqual(expPathExtDict, actPathExtDict)


	@mock.patch("gen_basis_helpers.db_help.md_files_help.shutil.copy2")
	@mock.patch("gen_basis_helpers.db_help.md_files_help.os.listdir")
	def testCopyWfnRestart_multipleBackupFilesReqd(self, mockListDir, mockCopy):
		runDir = "fake_run_dir"
		outDir = "fake_out_dir"
		baseFileName = "fake_file.wfn"
		mockDirContents = [baseFileName, baseFileName+".bak-2", baseFileName+".bak-1", baseFileName+".bak-3"]
		mockListDir.side_effect = lambda *args, **kwargs: mockDirContents

		self.testObjA.maxNumbWfnBackups = 2
		self.testObjA._copyWfnRestartFiles(runDir, outDir, baseFileName)

		mockCopy.assert_any_call( os.path.join(runDir,baseFileName), os.path.join(outDir,baseFileName) )
		mockCopy.assert_any_call( os.path.join(runDir,baseFileName+".bak-1"), os.path.join(outDir,baseFileName+".bak-1") )
		with self.assertRaises(AssertionError):
			mockCopy.assert_any_call( os.path.join(runDir,baseFileName+".bak-2"), os.path.join(outDir, baseFileName+".bak-2") )


	@mock.patch("gen_basis_helpers.db_help.md_files_help.shutil.copy2")
	@mock.patch("gen_basis_helpers.db_help.md_files_help.os.listdir")
	def testCopyWfnRestart_reqdNumberFilesGreaterThanActual(self, mockListDir, mockCopy):
		runDir = "fake_run_dir"
		outDir = "fake_out_dir"
		baseFileName = "fake_file.wfn"
		mockDirContents = [baseFileName, baseFileName+".bak-2", baseFileName+".bak-1"]
		mockListDir.side_effect = lambda *args, **kwargs: mockDirContents

		self.testObjA.maxNumbWfnBackups = 10
		self.testObjA._copyWfnRestartFiles(runDir, outDir, baseFileName)

		mockCopy.assert_any_call( os.path.join(runDir, baseFileName), os.path.join(outDir, baseFileName) )
		mockCopy.assert_any_call( os.path.join(runDir, baseFileName+".bak-1"), os.path.join(outDir, baseFileName+".bak-1") )
		mockCopy.assert_any_call( os.path.join(runDir, baseFileName+".bak-2"), os.path.join(outDir, baseFileName+".bak-2") )


class TestDumpMetadynHillsFile(unittest.TestCase):

	def setUp(self):
		self.startDir = "fake_dir"
		self.hillsInfoObj = mock.Mock()
		self.pathExt = "path_extension"
		self.fileName = "fake_filename.file"
		self.createTestObjs()

	def createTestObjs(self):
		pass

	def _runTestFunct(self):
		args = [self.startDir, self.hillsInfoObj]
		kwargs = {"pathExt":self.pathExt, "fileName":self.fileName}
		return tCode.dumpMetadynHillsFileAndGetFileDictStd(*args, **kwargs)

	@mock.patch("gen_basis_helpers.db_help.md_files_help._dumpMetadynHillsToFile")
	def testExpectedCallsToDumpFunctA(self, mockedDumpToFile):
		expPath = os.path.join(self.startDir, self.pathExt, self.fileName)
		expDict = {"metadyn_hills_path_ext": self.pathExt, "metadyn_hills_filename": self.fileName}
		actDict = self._runTestFunct()
		mockedDumpToFile.assert_called_with(expPath, self.hillsInfoObj)
		self.assertEqual(expDict,actDict)

class TestDumpHillsToFileBackend(unittest.TestCase):

	def setUp(self):
		self.inpDir = "test_dir"
		self.fileName = "fake_filename"
		self.hillsObj = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.inpPath = os.path.join(self.inpDir, self.fileName)

	def _runTestFunct(self):
		tCode._dumpMetadynHillsToFile(self.inpPath, self.hillsObj)

	@mock.patch("gen_basis_helpers.db_help.md_files_help.ioHelp.dumpObjWithToDictToJson")
	def testExpectedCallsMade(self, mockDumpToJson):
		self._runTestFunct()
		mockDumpToJson.assert_called_with(self.hillsObj, self.inpPath)


class TestGetMetadynHillsObj(unittest.TestCase):

	def setUp(self):
		self.baseDbFolder = "test_folder_path"
		self.dbExt = "db_ext"
		self.hillPathExt = "hill_ext"
		self.hillFilename = "hill_file.json"
		self.createTestObjs()

	def createTestObjs(self):
		self.inpRecord = {"db_ext_path":self.dbExt, "metadyn_hills_path_ext": self.hillPathExt,
		                  "metadyn_hills_filename": self.hillFilename}

	def _runTestFunct(self):
		return tCode.getMetadynHillsInfoObjFromBaseDbFolderAndRecord(self.baseDbFolder, self.inpRecord)

	@mock.patch("gen_basis_helpers.db_help.md_files_help.ioHelp.readObjWithFromDictFromJsonFile")
	def testExpectedCallsMade(self, mockReadObjFromFile):
		expObj = mock.Mock()
		expFilePath = os.path.join(self.baseDbFolder, self.dbExt, self.hillPathExt, self.hillFilename)
		mockReadObjFromFile.side_effect = lambda *args, **kwargs: expObj

		actObj = self._runTestFunct()

		mockReadObjFromFile.assert_called_with( aMetadynHillsHelp.MetadynHillsInfo, expFilePath )
		self.assertEqual(expObj, actObj)



