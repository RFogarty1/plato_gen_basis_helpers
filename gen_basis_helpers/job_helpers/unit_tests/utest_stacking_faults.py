
import itertools as it

import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.stacking_faults as tCode



class TestStandardInputCreatorObj(unittest.TestCase):

	def setUp(self):
		self.baseCreator = mock.Mock()
		self.perfectCellBaseCreator = mock.Mock()
		self.perfectCellGeom = mock.Mock()
		self.dispZeroGeom = mock.Mock()
		self.dispVals = [1,2] #Lists are anoying to use mocks for
		self.fitterObj = mock.Mock()
		self.stackingFaultGeomGenerator = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"baseCreator":self.baseCreator, "perfectCellBaseCreator":self.perfectCellBaseCreator,
		             "perfectCellGeom":self.perfectCellGeom,
		             "dispZeroGeom":self.dispZeroGeom, "dispVals":self.dispVals,
		             "fitterObj":self.fitterObj, "stackingFaultGeomGenerator":self.stackingFaultGeomGenerator}
		self.testObjA = tCode.StandardInputCreatorTemplate(**kwargDict)

	def testExpectedCallsToCreateGeoms(self):
		mainFunct = self.stackingFaultGeomGenerator.getGeomForGivenDisplacement
		mainFunct.side_effect = lambda argA,dispVal: dispVal
		expGeoms  = [x for x in self.dispVals]
		actGeoms = self.testObjA._createAllDispGeoms()

		for x in self.dispVals:
			mainFunct.assert_any_call(self.dispZeroGeom,x)

		[self.assertEqual(exp,act) for exp,act in it.zip_longest(expGeoms,actGeoms)]

	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults.StandardInputCreatorTemplate.outFolder", new_callable=mock.PropertyMock)
	def testGetKwargDictToModPerfectCell(self, mockOutFolderProp):
		expWorkfolder = mock.Mock()
		mockOutFolderProp.return_value = expWorkfolder
		expDict = dict()
		expDict["geom"] = self.perfectCellGeom
		expDict["workFolder"] = expWorkfolder
		expDict["fileName"] = "perfect_cell"
		actDict = self.testObjA._getKwargDictForModdingPerfectCellCreatorObj()
		for key in expDict.keys():
			self.assertEqual(expDict[key],actDict[key])

	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults.StandardInputCreatorTemplate._getKwargDictForModdingPerfectCellCreatorObj")
	def testCreatePerfectCalcObj(self, mockedGetKwargDict):
		expOutObj = mock.Mock()
		expDict = {"fake_key":1}
		mockedGetKwargDict.side_effect = lambda *args:expDict
		self.perfectCellBaseCreator.create.side_effect = lambda *args,**kwargs:expOutObj
		actOutObj = self.testObjA._createPerfectCalcObj()
		self.perfectCellBaseCreator.create.assert_called_with(**expDict)
		self.assertEqual(expOutObj,actOutObj)

	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults.StandardInputCreatorTemplate.outFolder", new_callable=mock.PropertyMock)
	def testGetKwargDictToModDispCell(self, mockedOutFolderProp):
		dispVal = 0.5
		expWorkfolder = mock.Mock()
		expGeom = mock.Mock()
		mockedOutFolderProp.return_value = expWorkfolder
		expDict = dict()
		expDict["geom"] = expGeom
		expDict["workFolder"] = expWorkfolder
		expDict["fileName"] = "disp_val_0pt5"
		actDict = self.testObjA._getKwargDictForModdingDispCreatorObj(expGeom,dispVal)
		for key in expDict.keys():
			self.assertEqual( expDict[key], actDict[key] )


	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults.StandardInputCreatorTemplate._getKwargDictForModdingDispCreatorObj")
	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults.StandardInputCreatorTemplate._createAllDispGeoms")
	def testCreateDispCalcObjs(self, mockedGetGeoms, mockedGetKwargDict):
		expGeoms = self.dispVals
		expDicts = {"testKey":"testVal"}
		mockedGetGeoms.side_effect = lambda *args: expGeoms
		mockedGetKwargDict.side_effect = lambda geom,dispVal:expDicts
		self.baseCreator.create.side_effect = self.dispVals

		expCalcObjs = self.dispVals
		actCalcObjs = self.testObjA._createAllDispCalcObjs()

		for dVal in self.dispVals:
			mockedGetKwargDict.assert_any_call(dVal, dVal)
		self.baseCreator.create.assert_any_call(**expDicts)

		[self.assertEqual(exp,act) for exp,act in it.zip_longest(expCalcObjs,actCalcObjs)]

	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults.StandardInputCreatorTemplate._createAllDispCalcObjs")
	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults.StandardInputCreatorTemplate._createPerfectCalcObj")
	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows.StackingFaultWorkflow")
	def testCreateWorkflow(self, mockedflow, mockGetPerfObj, mockGetDispObjs):
		mockPerfectObj, mockDispObjs, mockedDispVals = mock.Mock, mock.Mock(), mock.Mock()
		expOutflow = mock.Mock()
		mockGetPerfObj.side_effect = lambda *args: mockPerfectObj
		mockGetDispObjs.side_effect = lambda *args: mockDispObjs
		mockedflow.side_effect = lambda *args,**kwargs: expOutflow
		actOutflow = self.testObjA._createWorkflow()	
		mockedflow.assert_called_with(mockPerfectObj, mockDispObjs, self.dispVals, self.fitterObj)
		self.assertEqual(expOutflow, actOutflow)	

	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults.calcRunners.StandardInputObj")
	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults.StandardInputCreatorTemplate.label", new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_helpers.stacking_faults.StandardInputCreatorTemplate._createWorkflow")
	def testCreate(self, mockedGetWorkflow, mockedGetLabel, mockedStdInp):
		expLabel, expWorkflow = mock.Mock(), mock.Mock()
		expOutObj = mock.Mock()
		
		mockedGetLabel.return_value = expLabel
		mockedGetWorkflow.side_effect = [expWorkflow]
		mockedStdInp.side_effect = lambda *args: expOutObj
		actOutObj = self.testObjA.create()

		mockedStdInp.assert_called_with(expWorkflow,expLabel)
		self.assertEqual(expOutObj,actOutObj)




