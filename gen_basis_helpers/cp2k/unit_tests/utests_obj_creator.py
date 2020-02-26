
import os

import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.cp2k_creator as tCode

class TestStandardCreationObj(unittest.TestCase):

	def setUp(self):
		self.kPts = [20,20,20]
		self.addedMOs = 5
		self.geom = "fake_geom"
		self.basisObj = "fake_basis_obj"
		self.workFolder = "fake_folder"
		self.fileName = "fake_file_name"
		self.methodStr = "cp2k_test_object"
		self.createTestObjs()

	def createTestObjs(self):
		self.testCreatorObjA = tCode.CP2KCalcObjFactoryStandard(methodStr=self.methodStr, kPts=self.kPts,
		                                                        addedMOs=self.addedMOs, geom=self.geom,
		                                                        basisObjs=self.basisObj, folderPath=self.workFolder,
		                                                        fileName=self.fileName)

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
		expArgDict = {"kpts":self.kPts, "addedMOs".lower():self.addedMOs}
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


	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.calcObjs.os.path.abspath")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator.fileHelpers")
	def testOutputFileSetProperly(self, mockedFileHelpers, absPathMock):
		absPathMock.side_effect = lambda x: x
		outObj = self.testCreatorObjA.create()
		expBasePath = os.path.join( self.workFolder, self.fileName ) #NOTE: This should be path WITHOUT the file extension
		actBasePath = outObj.basePath
		self.assertEqual(expBasePath, actBasePath)


