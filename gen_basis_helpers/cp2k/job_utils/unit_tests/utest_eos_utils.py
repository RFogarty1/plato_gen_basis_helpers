

import itertools as it
import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.job_utils.eos_utils as tCode
import gen_basis_helpers.shared.label_objs as labelHelp


class TestEosWorkflowCreator(unittest.TestCase):

	def setUp(self):
		self.calcObjCreator = mock.Mock()
		self.structStrParamMapper = mock.Mock()
		self.eosFitFunction = mock.Mock()

#		self.geomList = [mock.Mock(),mock.Mock()] 
		self.volList = [20,30] #Needed in the file names i guess
		self.structStrs = ["hcp","bcc"]
		self.kPtVals = [ [20,20,12], [20,20,20] ]
		self.eleKey, self.methodKey = "Mg", "some_method"
		self.createTestObjs()

	def createTestObjs(self):
		#Modify the mocks i guess is needed
		self.geomsA = [mock.Mock() for x in self.volList]
		self.geomsA[0].volume = mock.Mock()
		self.geomsA[0].volume.return_value = self.volList[0]
		self.geomsA[0].volume.side_effect = lambda : 20
		for idx,x in enumerate(self.geomsA):
			x.volume = self.volList[idx]

		getKwargsFunctA = lambda key: {k:{"kpts":v} for k,v in it.zip_longest(self.structStrs,self.kPtVals)}[key]

		self.structStrParamMapper.getGeomsForStructStr.side_effect = lambda key: self.geomsA
		self.structStrParamMapper.getKwargDictForStructStr.side_effect = getKwargsFunctA


		self.testObjA = tCode.CP2KEosWorkflowCreator( calcObjCreator=self.calcObjCreator,
		                                              structStrParamMapper=self.structStrParamMapper,
		                                              eosFitFunction=self.eosFitFunction, eleKey=self.eleKey,
		                                              methodKey=self.methodKey, structStrs= self.structStrs )

	@mock.patch("gen_basis_helpers.cp2k.job_utils.eos_utils.eosFlow")
	def testSingleWorkflowCreation_callsToCalcObjs(self, mockedEosFlow):
		testStructStr = self.structStrs[0]
		outWorkflow = self.testObjA._createSingleWorkflow(testStructStr)

		expKwargDicts = list()
		for vol,geom in it.zip_longest(self.volList,self.geomsA):
			expDict = {"geom": geom, "kpts":self.kPtVals[0], "fileName": "vol_{:.3f}".format(vol).replace(".","pt")}
			expKwargDicts.append( expDict )

		for kDict in expKwargDicts:
			self.calcObjCreator.create.assert_any_call(**kDict)

	@mock.patch("gen_basis_helpers.cp2k.job_utils.eos_utils.eosFlow")
	def testSingleWorkflowCreation_CorrectArgsPassedToWorkflow(self, mockedEosFlow):
		testStructStr = self.structStrs[0]
		expLabel = labelHelp.StandardLabel( eleKey=self.eleKey, structKey=testStructStr, methodKey=self.methodKey )
		expArgs = [self.eosFitFunction, expLabel]
		self.testObjA._createSingleWorkflow(testStructStr)
		actArgs, unused_kwargs = mockedEosFlow.EosWorkflow.call_args
		self.assertEqual(expArgs, list(actArgs[1:])) #Dont need or want to test what gets passed as calcObjs


	@mock.patch("gen_basis_helpers.cp2k.job_utils.eos_utils.baseFlow.StandardLabelledWorkflowComposite")
	@mock.patch("gen_basis_helpers.cp2k.job_utils.eos_utils.CP2KEosWorkflowCreator._createSingleWorkflow")
	def testCreateMultipleWorkflows(self, mockedSingleWorkflowCreator, mockedCompositeWorkflow):
		mockedSingleWorkflowCreator.side_effect = lambda x: x #so we get back our structStr
		self.testObjA.create()
		#Check that we called for a single workflow with correct args
		for structStr in self.structStrs:
			mockedSingleWorkflowCreator.assert_any_call(structStr)

		#Check that we properly merged the workflows into a composite
		expCompArgs = [x for x in self.structStrs] #Since this is returned by our mocked single workflow
		mockedCompositeWorkflow.assert_called_with(expCompArgs)



class TestStructStrToParamMapper(unittest.TestCase):

	def setUp(self):
		self.structDatabase = mock.Mock()
		self.convDatabase = mock.Mock()
		self.fakeGeomsA = ["geom_a","geom_b"]
		self.modDictA = {"kpts":[20,20,12]}
		self.structDatabase.getStructsForEos.side_effect = lambda x: self.fakeGeomsA
		self.convDatabase.kptGridVals.getKptsPrimCell.side_effect = lambda x: self.modDictA["kpts"]
		self.testStructStrA = "hcp"
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.StructStrToParamMapper(self.structDatabase, self.convDatabase)

	def testCorrectGeomsObtained(self):
		expGeoms = self.fakeGeomsA
		actGeoms = self.testObjA.getGeomsForStructStr( self.testStructStrA )
		self.structDatabase.getStructsForEos.assert_called_once_with(self.testStructStrA)
		self.assertEqual(expGeoms, actGeoms)

	def testCorrectModDictObtained(self):
		expOutDict= self.modDictA
		actOutDict = self.testObjA.getKwargDictForStructStr(self.testStructStrA)
		self.convDatabase.kptGridVals.getKptsPrimCell.assert_called_once_with(self.testStructStrA)
		self.assertEqual(expOutDict, actOutDict)

