
import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.stacking_faults_two_structs as tCode

class TestStackingFaultTwoStructs(unittest.TestCase):

	def setUp(self):
		self.baseCreatorPerfectStruct = mock.Mock()
		self.baseCreatorStackFaultStruct = mock.Mock()
		self.perfectStructGeom = mock.Mock()
		self.stackFaultGeom = mock.Mock()
		self.eleKey = "eleKey"
		self.structKey = "structKey"
		self.methodKey = "methodKey"
		self.baseWorkFolder = "fake_folder"
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"baseCreatorPerfectStruct": self.baseCreatorPerfectStruct,
		             "baseCreatorStackFaultStruct": self.baseCreatorStackFaultStruct,
		             "baseWorkFolder":self.baseWorkFolder, "eleKey":self.eleKey, "methodKey":self.methodKey,
		             "perfectStructGeom":self.perfectStructGeom, "stackFaultGeom":self.stackFaultGeom,
		             "structKey":self.structKey}
		self.testObjA = tCode.StandardInputCreatorStackFaultTwoStructs(**kwargDict)

	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults_two_structs.StandardInputCreatorStackFaultTwoStructs.outFolder", new_callable=mock.PropertyMock)
	def testExpectedCreatorMods_perfectStruct(self, mockOutFolderProp):
		expOutFolder = mock.Mock()
		mockOutFolderProp.return_value = expOutFolder
		expDict = dict()
		expDict["geom"] = self.perfectStructGeom
		expDict["workFolder"] = expOutFolder
		expDict["fileName"] = "no_fault_struct"
		actDict = self.testObjA._getKwargDictForModdingPerfectStructCreator()
		for key in expDict.keys():
			self.assertEqual( expDict[key], actDict[key] )

	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults_two_structs.StandardInputCreatorStackFaultTwoStructs.outFolder", new_callable=mock.PropertyMock)
	def testExpectedCreatorMods_stackStruct(self, mockOutFolderProp):
		expOutFolder = mock.Mock()
		mockOutFolderProp.return_value = expOutFolder
		expDict = dict()
		expDict["geom"] = self.stackFaultGeom
		expDict["workFolder"] = expOutFolder
		expDict["fileName"] = "stack_fault_geom"
		actDict = self.testObjA._getKwargDictForModdingStackFaultCreator()
		for key in expDict.keys():
			self.assertEqual( expDict[key], actDict[key] )


	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults_two_structs.StandardInputCreatorStackFaultTwoStructs._getCalcObjStackFaultStruct")
	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults_two_structs.StandardInputCreatorStackFaultTwoStructs._getCalcObjPerfectStruct")
	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults_two_structs.stackFaultFlow.StackingFaultWorkflowTwoStructs")
	def testExpectedArgsPassedToWorkflow(self, mockedWorkflow, mockedGetPerfStruct, mockedGetFaultStruct):
		expPerfStructObj, expFaultStructObj, expWorkflow = mock.Mock(), mock.Mock(), mock.Mock()
		mockedGetPerfStruct.side_effect = lambda: expPerfStructObj
		mockedGetFaultStruct.side_effect = lambda: expFaultStructObj
		mockedWorkflow.side_effect = lambda *args:expWorkflow
		actWorkflow = self.testObjA._createWorkflow()	
		mockedWorkflow.assert_called_with(expPerfStructObj,expFaultStructObj)
		self.assertEqual(expWorkflow, actWorkflow)


	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults_two_structs.calcRunners.StandardInputObj")
	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults_two_structs.StandardInputCreatorStackFaultTwoStructs.label", new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults_two_structs.StandardInputCreatorStackFaultTwoStructs._createWorkflow")
	def testExpectedCallToStdInpCreator(self, mockedGetWorkflow, mockedLabel, mockedStandardInp):
		expWorkflow, expLabel, expStdInp = mock.Mock(), mock.Mock(), mock.Mock()
		mockedGetWorkflow.side_effect = lambda : expWorkflow
		mockedLabel.return_value = expLabel
		mockedStandardInp.side_effect = lambda *args:expStdInp
		actStdInp = self.testObjA.create()
		mockedStandardInp.assert_called_with(expWorkflow, expLabel)
		self.assertEqual(expStdInp, actStdInp)

