
import copy
import itertools as it
import os

import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.cp2k_creator as tCode

import gen_basis_helpers.shared.geom_constraints as geomConstr

class TestStandardCreationObj(unittest.TestCase):

	def setUp(self):
		self.kPts = [20,20,20]
		self.addedMOs = 5
		self.geom = "fake_geom"
		self.basisObj = "fake_basis_obj"
		self.xcFunctional = None
		self.workFolder = "fake_folder"
		self.fileName = "fake_file_name"
		self.methodStr = "cp2k_test_object"
		self.grimmeDisp = None
		self.fragmentsBSSE = None
		self.runType = None
		self.printAOMullikenPop = False
		self.mdOpts = None
		self.walltime = None
		self.extrapolationMethod = None
		self.print_every_n_md_steps = None
		self.print_every_n_scf_steps = None
		self.restart_file_every_n_md_steps = None
		self.prefDiagLib = None
		self.epsDef = None
		self.nGrids = None
		self.colVars = None
		self.metaDynOpts = None
		self.thermostatOpts = None
		self.rsGridDistrib = None
		self.scfMixAlpha = None
		self.scfMixMethod = None
		self.scfOTMinimizer = None
		self.scfOTEnergies = None
		self.scfOTRotation = None
		self.scfGuess = None
		self.scfPrintRestartHistoryOn = None
		self.scfPrintRestartHistory_eachMD = None
		self.scfPrintRestartHistory_eachSCF = None

		self.createTestObjs()

	#Note we pass the None value for workFolder as a test essentially; if EITHER folderPath or workFolder are set to a real (not None) value then we take that one for both
	def createTestObjs(self):
		self.testCreatorObjA = tCode.CP2KCalcObjFactoryStandard(methodStr=self.methodStr, kPts=self.kPts,
		                                                        addedMOs=self.addedMOs, geom=self.geom,
		                                                        basisObjs=self.basisObj, folderPath=self.workFolder,
		                                                        fileName=self.fileName, workFolder=None, printAOMullikenPop=self.printAOMullikenPop,
		                                                        runType=self.runType, fragmentsBSSE=self.fragmentsBSSE, xcFunctional=self.xcFunctional,
		                                                        grimmeDisp=self.grimmeDisp, mdOpts=self.mdOpts, walltime=self.walltime,
		                                                        extrapolationMethod=self.extrapolationMethod, print_every_n_md_steps=self.print_every_n_md_steps,
		                                                        print_every_n_scf_steps=self.print_every_n_scf_steps,
		                                                        restart_file_every_n_md_steps=self.restart_file_every_n_md_steps,
		                                                        prefDiagLib=self.prefDiagLib, epsDef=self.epsDef, nGrids=self.nGrids,
		                                                        colVars=self.colVars, metaDynOpts=self.metaDynOpts, thermostatOpts=self.thermostatOpts,
		                                                        rsGridDistrib=self.rsGridDistrib, scfMixAlpha=self.scfMixAlpha, scfMixMethod=self.scfMixMethod,
		                                                        scfOTMinimizer=self.scfOTMinimizer, scfOTEnergies=self.scfOTEnergies,scfOTRotation=self.scfOTRotation,
		                                                        scfGuess=self.scfGuess, scfPrintRestartHistoryOn=self.scfPrintRestartHistoryOn,
		                                                        scfPrintRestartHistory_eachMD=self.scfPrintRestartHistory_eachMD,
		                                                        scfPrintRestartHistory_eachSCF=self.scfPrintRestartHistory_eachSCF)

	def testWrongKwargCaughtByInit(self):
		with self.assertRaises(KeyError):
			tCode.CP2KCalcObjFactoryStandard(fake_kwarg=None)

	def testWrongKwargCaughtByCreate(self):
		with self.assertRaises(KeyError):
			self.testCreatorObjA.create(fake_kwarg=None)

	def testReqArgsToBeSetAsExpected(self):
		baseReqArgs = set(self.testCreatorObjA._baseReqArgsToBeSet)
		additional = ["very_fake_arg"]
		expReqArgs = baseReqArgs.union(additional)
		self.assertNotEqual( baseReqArgs.union(additional), baseReqArgs )
		self.testCreatorObjA.requiredArgsToBeSet = additional
		actReqArgs = self.testCreatorObjA.requiredArgsToBeSet
		self.assertEqual( list(expReqArgs), list(actReqArgs) )

	def testCreateRaisesWhenReqArgIsMissing(self):
		self.assertTrue( "methodStr" in self.testCreatorObjA.requiredArgsToBeSet )
		self.methodStr = None
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testCreatorObjA.create()

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.methRegister")
	def testCorrectArgsPassedToMethodReg(self, methRegisterMock, mockFileHelpers):
		self.testCreatorObjA.create()
		methRegisterMock.createCP2KObjFromMethodStr.assert_called_once_with(self.methodStr)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.methRegister")
	def testSelectedArgsPassedToFileHelpers(self, mockMethReg, mockFileHelpers):
		self.xcFunctional = "BLYP"
		self.createTestObjs()
		expArgDict = {"kpts":self.kPts, "addedMOs".lower():self.addedMOs, "xcFunctional".lower():"BLYP"}
		self.testCreatorObjA.create()
		args,kwargs = mockFileHelpers.modCp2kObjBasedOnDict.call_args
		actArgDict = {k.lower():args[1][k] for k in expArgDict.keys()}
		self.assertEqual(expArgDict, actArgDict)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.methRegister")
	def testGeomAndBasisPassedToFileHelpers(self, mockMethReg, mockFileHelpers):
		self.testCreatorObjA.create()
		geoArgs  , kwargs = mockFileHelpers.addGeomInfoToSimpleCP2KObj.call_args
		basisArgs, kwargs = mockFileHelpers.addBasisInfoToSimpleCP2KObj.call_args
		self.assertEqual(self.geom,geoArgs[1])
		self.assertEqual(self.basisObj,basisArgs[1]) 

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.methRegister")
	def testPrintMullikenPopsPassedToFileHlpers(self, mockMethReg, mockFileHelpers):
		expRelevantArgDict = {"printAOMullikenPop".lower():True}
		self.testCreatorObjA.create(printAOMullikenPop=True)
		args,kwargs = mockFileHelpers.modCp2kObjBasedOnDict.call_args
		actRelevantArgDict = {k.lower():args[1][k] for k in expRelevantArgDict.keys()}
		self.assertEqual(expRelevantArgDict, actRelevantArgDict)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.methRegister")
	def testChagePassedToFileHelpers(self, mockMethReg, mockFileHelpers):
		expRelevantArgDict = {"charge".lower():2}
		self.testCreatorObjA.create(charge=2)
		args,kwargs = mockFileHelpers.modCp2kObjBasedOnDict.call_args
		actRelevantArgDict = {k.lower():args[1][k] for k in expRelevantArgDict.keys()}
		self.assertEqual(expRelevantArgDict, actRelevantArgDict)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.methRegister")
	def testGrimmeCoorsPassedToFileHelpers(self, mockMethReg, mockFileHelpers):
		expDict = {k.lower():v for k,v in {"keyA":"valA", "keyB":"valB"}.items()}
		stubCorrObj = types.SimpleNamespace( modPyCP2KDict=expDict )
		self.testCreatorObjA.create(grimmeDisp=stubCorrObj)
		args, kwargs = mockFileHelpers.modCp2kObjBasedOnDict.call_args
		actModDict = args[1]
		actRelevantArgDict = {k.lower():actModDict[k] for k in expDict.keys()}
		self.assertEqual(expDict, actRelevantArgDict)


	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.calcObjs.os.path.abspath")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	def testOutputFileSetProperlyBySettingFolderPath(self, mockedFileHelpers, absPathMock):
		absPathMock.side_effect = lambda x: x
		outObj = self.testCreatorObjA.create()
		expBasePath = os.path.join( self.workFolder, self.fileName ) #NOTE: This should be path WITHOUT the file extension
		actBasePath = outObj.basePath
		self.assertEqual(expBasePath, actBasePath)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.calcObjs.os.path.abspath")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	def testOutputFileSetProperlyBySettingWorkfolder(self, mockedFileHelpers, absPathMock):
		absPathMock.side_effect = lambda x: x
		testFolderPath = "fake_folder_b"
		self.assertNotEqual(testFolderPath, self.testCreatorObjA.folderPath)
		self.testCreatorObjA.workFolder = testFolderPath
		outObj = self.testCreatorObjA.create()
		expBasePath = os.path.join( testFolderPath, self.fileName)
		actBasePath = outObj.basePath
		self.assertEqual(expBasePath, actBasePath)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.calcObjs.os.path.abspath")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	def testOutputFileSetProperlyBySettingWorkfolderAtCreateTime(self, mockedFileHelpers, absPathMock):
		absPathMock.side_effect = lambda x: x
		testFolderPath = "fake_folder_b"
		self.assertNotEqual(testFolderPath, self.testCreatorObjA.folderPath)
		outObj = self.testCreatorObjA.create(workFolder=testFolderPath) 
		expBasePath = os.path.join( testFolderPath, self.fileName)
		actBasePath = outObj.basePath
		self.assertEqual(expBasePath, actBasePath)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.methRegister")
	def testCellOptPassedToFileHelpers(self, mockMethReg, mockFileHelpers):
		self.runType = "geomOpt"
		self.createTestObjs()
		expArgDict = {"runtype":"cell_opt"}
		self.testCreatorObjA.create()
		args,kwargs = mockFileHelpers.modCp2kObjBasedOnDict.call_args
		actArgDict = {k.lower():args[1][k] for k in expArgDict.keys()}
		self.assertEqual(expArgDict, actArgDict)

	def testExpectedModDictForRunTypeBSSE(self):
		self.runType = "bsse"
		self.fragmentsBSSE = mock.Mock()
		self.createTestObjs()
		expArgDict = {"runType".lower():"bsse", "fragmentsBSSE".lower():self.fragmentsBSSE}
		actArgDict = self.testCreatorObjA._getModDictBasedOnRunType()
		self.assertEqual(expArgDict, actArgDict)

	def testExpectedExtraMDOptsForMdRun(self):
		self.runType = "md"
		self.mdOpts = mock.Mock()
		mdDict = {"fake_key_a":"fake_val_a"}
		self.mdOpts.optDict = mdDict
		expDict = copy.deepcopy(mdDict)
		expDict["runType".lower()] = "md"
		self.createTestObjs()
		actDict = self.testCreatorObjA._getModDictBasedOnRunType()
		self.assertEqual(expDict, actDict)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.methRegister")
	def testExpectedMiscOptsAPassedToFileHelpers(self, mockMethReg, mockFileHelpers):
		self.walltime = 20
		self.extrapolationMethod = "fake_extrapolation_method"
		self.print_every_n_md_steps = 45
		self.print_every_n_scf_steps = 20
		self.restart_file_every_n_md_steps = 40
		self.epsDef = 1e-14
		self.prefDiagLib = "SL"
		self.nGrids = 8
		self.colVars = mock.Mock()
		self.thermostatOpts = mock.Mock()
		self.rsGridDistrib = [-1,-1,24]
		self.scfMixAlpha = 5
		self.scfMixMethod = "pulay"

		self.scfOTMinimizer = "DIIS"
		self.scfOTEnergies = True
		self.scfOTRotation = True
		self.scfGuess = "restart"
		self.scfPrintRestartHistoryOn = True
		self.scfPrintRestartHistory_eachMD = 5
		self.scfPrintRestartHistory_eachSCF = 7

		self.createTestObjs()
		expArgDict = {"qsExtrapolationMethod".lower(): self.extrapolationMethod, "walltime": self.walltime,
		               "trajPrintEachMd".lower(): self.print_every_n_md_steps,
		               "trajPrintEachScf".lower(): self.print_every_n_scf_steps,
		               "restartPrintEachMd".lower(): self.restart_file_every_n_md_steps,
		               "prefDiagLib".lower(): "SL",
		               "epsDef".lower():1e-14,
		               "nGrids".lower():self.nGrids,
		               "colVars".lower():self.colVars,
		               "mdThermoStatOpts".lower():self.thermostatOpts,
		               "rsGrid_distrib".lower():self.rsGridDistrib,
		               "scfMixAlpha".lower(): self.scfMixAlpha,
		               "scfMixMethod".lower(): self.scfMixMethod,
		               "scfOTMinimizer".lower():self.scfOTMinimizer,
		               "scfOTRotation".lower():self.scfOTRotation,
		               "scfGuess".lower():self.scfGuess,
		               "scfPrintRestartHistoryOn".lower():self.scfPrintRestartHistoryOn,
		               "scfPrintRestartHistory_eachMD".lower():self.scfPrintRestartHistory_eachMD,
		               "scfPrintRestartHistory_eachSCF".lower():self.scfPrintRestartHistory_eachSCF}

		self.testCreatorObjA.create()
		args,kwargs = mockFileHelpers.modCp2kObjBasedOnDict.call_args
		actArgDict = {k.lower():v for k,v in args[1].items()}

		for key in expArgDict:
			self.assertEqual( expArgDict[key], actArgDict[key] )

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.methRegister")
	def testExpectedDictFromMetaDynOptions(self, mockMethReg, mockFileHelpers):
		self.metaDynOpts = mock.Mock()
		self.createTestObjs()
		expArgDict = {"fake_arg_a":"fake_val_a", "fake_arg_b":"fake_val_b"}
		self.metaDynOpts.optDict = expArgDict
		self.testCreatorObjA.create()
		args,kwargs = mockFileHelpers.modCp2kObjBasedOnDict.call_args
		actArgDict = {k.lower():v for k,v in args[1].items()}

		for key in expArgDict:
			self.assertEqual( expArgDict[key], actArgDict[key] )


#TODO: Probably test atomic positions and cell constraints separately
class TestModDictBasedOnGeomConstraints(unittest.TestCase):

	def runTestFunct(self):
		return tCode.getGeoOptCP2KModDictBasedOnCellConstraints(self.testConstrObj.cellConstraints)

	def testNoConstraintsPresent(self):
		expDict = {"runType".lower():"cell_opt"} #Full optimisation, all variables free
		self.testConstrObj = geomConstr.GeomConstraints.initWithNoConstraints()
		actDict = self.runTestFunct()
		self.assertEqual(expDict, actDict)

	def testCellAnglesConstrained(self):
		self.testConstrObj = geomConstr.GeomConstraints.initWithNoConstraints()
		self.testConstrObj.cellConstraints.anglesToFix = [True,True,True]
		expDict = {"runtype".lower():"cell_opt",
		           "geo_constrain_cell_angles": [True,True,True]}
		actDict = self.runTestFunct()
		self.assertEqual(expDict,actDict)

	def testFullCellConstraintsPresent(self):
		allCellConstraints = [True,True,True]
		cellConstraints = geomConstr.CellConstraints(allCellConstraints, allCellConstraints)
		atomicPosConstraints=  geomConstr.AtomicPositionConstraints.initWithNoConstraints()
		self.testConstrObj = geomConstr.GeomConstraints(atomicPosConstraints, cellConstraints)
		expDict = {"runtype":"geo_opt"}
		actDict = self.runTestFunct()
		self.assertEqual(expDict,actDict)

	#TODO: Implement ability to actually fix lattice parameters independent of angles
	def testRaisesWhenOnlyLattParamsFixed(self):
		cellConstrs = geomConstr.CellConstraints.initWithNoConstraints()
		cellConstrs.lattParamsToFix = [True,True,True]
		atomicPosConstraints = geomConstr.AtomicPositionConstraints.initWithNoConstraints()
		self.testConstrObj = geomConstr.GeomConstraints(atomicPosConstraints,cellConstrs)
		with self.assertRaises(ValueError):
			self.runTestFunct()


class TestGetCP2KModDictBasedOnAtomicPosConstraints(unittest.TestCase):

	def setUp(self):
		self.atomsToFix = [0,1,2,3] #NOTE: We want to convert these into index-1 form
		self.componentsToFix = ["X","XY","XYZ","Z"]
		self.createTestObjs()

	def createTestObjs(self):
		self.geomConstrA = geomConstr.GeomConstraints.initWithNoConstraints()
		atomicCartConstr = list()
		for idx, comps in it.zip_longest(self.atomsToFix, self.componentsToFix):
			fixX = True if "X" in comps else False
			fixY = True if "Y" in comps else False
			fixZ = True if "Z" in comps else False
			currConstr = geomConstr.AtomicCartesianConstraint(idx, fixX=fixX, fixY=fixY, fixZ=fixZ)
			atomicCartConstr.append(currConstr)
		self.geomConstrA.atomicPositionConstraints.atomicCartConstraints = atomicCartConstr

	def _runTestFunct(self):
		return tCode._getModDictBasedOnAtomicPosConstraints(self.geomConstrA.atomicPositionConstraints)

	def testExpectedDictA(self):
		expDict = {"atPosConstraint_fixIdxPositions":[x+1 for x in self.atomsToFix],
		           "atPosConstraint_fixComponents":self.componentsToFix}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)


class TestGetHooksForCopyRestartFiles(unittest.TestCase):

	def setUp(self):
		self.restartName = "fake_file.restart"
		self.filePath = "fake_path.restart"
		self.workFolder = "fake_work_folder"
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"workFolder":self.workFolder, "inpRestartName":self.restartName,
		              "inpRestartPath":self.filePath}
		self.testObjA = tCode.CP2KCalcObjFactoryStandard(**currKwargs)

	def _runTestFunct(self):
		postWriteHooks = self.testObjA._getPostWriteFileHooks()
		assert len(postWriteHooks)==1
		postWriteHooks[0](mock.Mock())

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.shutil.copy2")
	def testExpectedCopyCommandPassed_workfolderSet(self, mockedCopyFunct):
		expOutPath = os.path.join(self.workFolder, self.restartName)
		expInpPath = self.filePath
		self._runTestFunct()
		mockedCopyFunct.assert_called_with(expInpPath, expOutPath)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.shutil.copy2")
	def testExpected_restartNameNone(self, mockedCopyFunct):
		self.restartName = None
		self.createTestObjs()
		expOutPath = os.path.join(self.workFolder, os.path.split(self.filePath)[-1])
		expInpPath = self.filePath
		self._runTestFunct()
		mockedCopyFunct.assert_called_with(expInpPath, expOutPath)

	def testExpected_workfolderNone(self):
		self.workFolder = None
		self.createTestObjs()
		with self.assertRaises(NotImplementedError):
			self._runTestFunct()

	def testExpected_inpPathNone(self):
		self.filePath = None
		self.createTestObjs()
		expVal = None
		actVal = self.testObjA._getPostWriteFileHooks()
		self.assertEqual(expVal, actVal)

