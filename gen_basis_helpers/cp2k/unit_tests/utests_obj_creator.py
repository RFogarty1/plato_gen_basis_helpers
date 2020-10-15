
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
		self.createTestObjs()

	#Note we pass the None value for workFolder as a test essentially; if EITHER folderPath or workFolder are set to a real (not None) value then we take that one for both
	def createTestObjs(self):
		self.testCreatorObjA = tCode.CP2KCalcObjFactoryStandard(methodStr=self.methodStr, kPts=self.kPts,
		                                                        addedMOs=self.addedMOs, geom=self.geom,
		                                                        basisObjs=self.basisObj, folderPath=self.workFolder,
		                                                        fileName=self.fileName, workFolder=None, printAOMullikenPop=self.printAOMullikenPop,
		                                                        runType=self.runType, fragmentsBSSE=self.fragmentsBSSE, xcFunctional=self.xcFunctional,
		                                                        grimmeDisp=self.grimmeDisp)

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
		args,kwargs = mockFileHelpers.addGeomAndBasisInfoToSimpleCP2KObj.call_args
		self.assertEqual(self.geom,args[1])
		self.assertEqual(self.basisObj,args[2]) 

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


#TODO: Probably test atomic positions and cell constraints separately
class TestModDictBasedOnGeomConstraints(unittest.TestCase):


	def runTestFunct(self):
		return tCode.getCP2KModDictBasedOnGeomConstraints(self.testConstrObj)

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

		

