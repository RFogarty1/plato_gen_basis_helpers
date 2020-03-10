
import os
import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.self_defects as tCode

import gen_basis_helpers.shared.label_objs as labelHelp

class TestSelfDefectsFactory(unittest.TestCase):

	def setUp(self):
		self.baseFolder = os.path.join("fake","path")
		self.methodKey = "fake_method"
		self.eleKey = "fake_ele"
		self.structKey = "fake_struct"
		self.bulkGeom = mock.Mock()
		self.defectGeom = mock.Mock()
		self.kPtsBulk = mock.Mock()
		self.kPtsDefect = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"baseWorkFolder":self.baseFolder, "eleKey":self.eleKey, "methodKey":self.methodKey,
		             "structKey":self.structKey, "kPtsBulk":self.kPtsBulk, "kPtsDefect":self.kPtsDefect,
		             "bulkGeom":self.bulkGeom, "defectGeom":self.defectGeom} 
		self.testObjA = tCode.CodeSpecificStandardInputCreatorTemplate(**kwargDict)

	def testOutfolder(self):
		expFolder = os.path.join( self.baseFolder, self.eleKey, self.structKey, self.methodKey )
		actFolder = self.testObjA.outFolder
		self.assertEqual(expFolder,actFolder)


	@mock.patch("gen_basis_helpers.job_helpers.self_defects.CodeSpecificStandardInputCreatorTemplate._createCalcObjCreatorBulk")
	@mock.patch("gen_basis_helpers.job_helpers.self_defects.CodeSpecificStandardInputCreatorTemplate.outFolder",new_callable=mock.PropertyMock)
	def testExpectedArgsPassedToCreatorForBulkCalc(self, mockedOutFolder, mockedBaseCreatorGetter):
		""" Make sure the args that this template control are passed to the basic creator object """
		baseCreator, expFolder = mock.Mock(), mock.Mock()
		mockedBaseCreatorGetter.side_effect = lambda *args,**kwargs: baseCreator
		mockedOutFolder.return_value = expFolder
		actCreator = self.testObjA._getBulkCreator()

		self.assertEqual(actCreator.workFolder, expFolder)
		self.assertEqual(actCreator.fileName, "bulk_calc")
		self.assertEqual(actCreator.kPts, self.kPtsBulk)
		self.assertEqual(actCreator.geom, self.bulkGeom)

	@mock.patch("gen_basis_helpers.job_helpers.self_defects.CodeSpecificStandardInputCreatorTemplate._createCalcObjCreatorDefect")
	@mock.patch("gen_basis_helpers.job_helpers.self_defects.CodeSpecificStandardInputCreatorTemplate.outFolder",new_callable=mock.PropertyMock)
	def testExpectedArgsPassedToCreatorForDefectCalc(self, mockedOutFolder, mockedBaseCreatorGetter):
		baseCreator, expFolder = mock.Mock(), mock.Mock()
		mockedBaseCreatorGetter.side_effect = lambda *args,**kwargs: baseCreator
		mockedOutFolder.return_value = expFolder
		actCreator = self.testObjA._getDefectCreator()

		self.assertEqual(actCreator.workFolder, expFolder)
		self.assertEqual(actCreator.fileName, "defect_calc")
		self.assertEqual(actCreator.kPts, self.kPtsDefect)
		self.assertEqual(actCreator.geom, self.defectGeom)

	@mock.patch("gen_basis_helpers.job_helpers.self_defects.calcRunners.StandardInputObj")
	@mock.patch("gen_basis_helpers.job_helpers.self_defects.defectFlow.SelfPointDefectWorkflow")
	@mock.patch("gen_basis_helpers.job_helpers.self_defects.CodeSpecificStandardInputCreatorTemplate._getDefectCreator")
	@mock.patch("gen_basis_helpers.job_helpers.self_defects.CodeSpecificStandardInputCreatorTemplate._getBulkCreator")
	def testCreateFromSelf(self, mockedBulkCreator, mockedDefectCreator, mockedDefectFlow, mockedStandardInput):
		bulkCreator, defectCreator = mock.Mock(), mock.Mock()
		expBulkObj, expDefectObj = mock.Mock(), mock.Mock()
		expWorkflow = mock.Mock()
		mockedBulkCreator.side_effect    = lambda *args,**kwargs: bulkCreator
		mockedDefectCreator.side_effect  = lambda *args,**kwargs: defectCreator
		bulkCreator.create.side_effect   = lambda *args,**kwargs: expBulkObj
		defectCreator.create.side_effect = lambda *args,**kwargs: expDefectObj
		mockedDefectFlow.side_effect = lambda *args,**kwargs: expWorkflow

		self.testObjA._createFromSelf()	
		expLabel = labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)
		mockedDefectFlow.assert_called_once_with(expDefectObj, expBulkObj)
		mockedStandardInput.assert_called_once_with(expWorkflow, expLabel)



