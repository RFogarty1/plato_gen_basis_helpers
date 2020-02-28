
import os
import itertools as it
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.job_utils.eos_utils as tCode
import gen_basis_helpers.shared.label_objs as labelHelp

import gen_basis_helpers.cp2k.cp2k_creator as creatorCode

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
		self.folderPath = "fake/path"
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

		self.calcObjCreator.folderPath = self.folderPath
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
			expDict = {"geom": geom, "kpts":self.kPtVals[0], "fileName": "vol_{:.3f}".format(vol).replace(".","pt"),
			           "folderPath":os.path.join(self.folderPath,testStructStr)}
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
		self.modDictA = {"kPts":[20,20,12]}
		self.structDatabase.getStructsForEos.side_effect = lambda x: self.fakeGeomsA
		self.convDatabase.kptGridVals.getKptsPrimCell.side_effect = lambda x: self.modDictA["kPts"]
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



class TestStandardInputCreator(unittest.TestCase):

	def setUp(self):
		self.eleStr = "Mg"
		self.basisStrDict = {"Mg":"fake_basis_str_mg", "H":"fake_basis_str_h"}
		self.basisAlias = "fake_alias"
		self.cp2kMethodStr = mock.Mock()
		self.structStrParamMapper = mock.Mock()
		self.structStrs = ["hcp","bcc"]
		self.absGridConv = 500
		self.relGridConv = 200
		self.addedMOs = 20
		self.workFolder = "fake/path/to/folder"
		self.eosStr = "birch" #Not actually a valid value (i dont think) for the standard fit funct
		self.maxFev = 200
		self._regKwargs = set(creatorCode.CP2KCalcObjFactoryStandard.registeredKwargs)
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.CP2KEosStandardObjCreator(workFolder=self.workFolder, absGridCutoff=self.absGridConv,
		                                                relGridCutoff=self.relGridConv, cp2kMethodStr=self.cp2kMethodStr,
		                                                structStrs=self.structStrs, basisStrDict=self.basisStrDict,
		                                                basisAlias=self.basisAlias, addedMOs=self.addedMOs,
		                                                eosStr=self.eosStr, maxFunctEvals=self.maxFev, eleStr=self.eleStr,
		                                                structStrParamMapper=self.structStrParamMapper)

	@mock.patch("gen_basis_helpers.cp2k.job_utils.eos_utils.basReg")
	def testCorrectCallsToBasisReg(self, mockedBasRegister):
		mockedBasRegister.createCP2KBasisObjFromEleAndBasisStr.side_effect = lambda eleKey,basStr: eleKey
		outBasisObjs = self.testObjA._getBasisObjects()

		#Check for calls
		for key,val in self.basisStrDict.items():
			mockedBasRegister.createCP2KBasisObjFromEleAndBasisStr.assert_any_call(key, val)	

		#Check for correct output
		expBasisObjs = sorted([k for k in self.basisStrDict.keys()])
		self.assertEqual( sorted(expBasisObjs), sorted(outBasisObjs) )

	#This is an implementation detail at current, but will probably be a hook or similar if i abstract it
	@mock.patch("gen_basis_helpers.cp2k.job_utils.eos_utils.cp2kCreator.CP2KCalcObjFactoryStandard")
	@mock.patch("gen_basis_helpers.cp2k.job_utils.eos_utils.CP2KEosStandardObjCreator._getBasisObjects")
	def testCorrectCreatorInitialised(self, mockedBasisGetter, mockedCreator):
		fakeBasisObjs = ["fake_basis"]
		fakeCreator = "fake_creator"
		mockedCreator.registeredKwargs = self._regKwargs
		mockedBasisGetter.side_effect = lambda : fakeBasisObjs
		mockedCreator.side_effect = lambda **kwargs: fakeCreator

		outCreator = self.testObjA._getCreatorObject()

		expKwargDict = {"methodStr":self.cp2kMethodStr, "addedMOs":self.addedMOs, "basisObjs":fakeBasisObjs,
		                "folderPath":self.workFolder, "absGridCutoff":self.absGridConv, "relGridCutoff":self.relGridConv}

		mockedCreator.assert_called_once_with(**expKwargDict)	

		self.assertEqual(fakeCreator,outCreator)

	@mock.patch("gen_basis_helpers.cp2k.job_utils.eos_utils.eosFlow.StandardEosFitFunction")
	def testFitFunctCorrectlyBuiltFromStandard(self, mockedFitClass):
		self.testObjA._getFitFunction()
		mockedFitClass.assert_called_with(eosStr=self.eosStr, maxFev=self.maxFev)

	#Probably a pretty pointless test actually
	@mock.patch("gen_basis_helpers.cp2k.job_utils.eos_utils.CP2KEosWorkflowCreator")
	def testCreateWorkflowFromCreator(self, mockedWorkflowCreator):
		expOutObj = "fake_workflow"
		mockCreator, mockFitFunct = mock.Mock(), mock.Mock()
		mockFactory = mock.Mock()
		mockFactory.create.side_effect = lambda **kwargs: expOutObj
		mockedWorkflowCreator.side_effect = lambda **kwargs : mockFactory
		outObj = self.testObjA._createWorkflowFromCreatorAndFitFunct(mockCreator, mockFitFunct)
		mockedWorkflowCreator.assert_called_with( calcObjCreator=mockCreator, structStrs=self.structStrs,structStrParamMapper=self.structStrParamMapper,
		                                          eosFitFunction=mockFitFunct , eleKey=self.eleStr, methodKey=self.basisAlias)
		self.assertEqual(expOutObj, outObj)

	def testCreateLabelAsExpected(self):
		expLabel = labelHelp.StandardLabel(eleKey=self.eleStr, structKey="eos", methodKey=self.basisAlias)
		actLabel = self.testObjA._createLabel()
		self.assertEqual(expLabel,actLabel)


class TestStandardMapFunction(unittest.TestCase):

	def setUp(self):
		self.plotDataA = [ [1,2], [2,4] ] #Each represents a different structure but same method
		self.plotDataB = [ [1,3], [2,6] ]
		self.methStr = "methA"
		self.eleKey, self.structKey = "Mg", "eos" #Dont actually matter here

		self.v0ValsA = [20,30]
		self.b0ValsA = [40,50]
		self.e0ValsA = [4,6]

		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.MapEosWorkflowOutputToUsefulFormatStandard()
		self.testOutputA = {"v0":self.v0ValsA[0], "b0":self.b0ValsA[0], "e0":self.e0ValsA[0],
		                    "data":self.plotDataA, "fitData":None, "fitAtDataPoints":None}
		self.testOutputB = {"v0":self.v0ValsA[1], "b0":self.b0ValsA[1], "e0":self.e0ValsA[1],
		                    "data":self.plotDataB, "fitData":None, "fitAtDataPoints":None}
		self.testTotalOutputAB = [types.SimpleNamespace(data=self.testOutputA), types.SimpleNamespace(data=self.testOutputB)] #workflows always return a list for outputs

		labelA = labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methStr)
		testWorkflowA = types.SimpleNamespace(output=self.testTotalOutputAB, run=mock.Mock()) #We only need the .output and for run to be called
		self.standardInpObjA = types.SimpleNamespace(workflow=testWorkflowA, label=[labelA])

	def runMapFunctOnDataA(self):
		return self.testObjA( self.standardInpObjA )

	def testPlotDataMapsProperly(self):
		expPlotOutput = [self.plotDataA, self.plotDataB]
		actTotalOutput = self.runMapFunctOnDataA()
		actPlotOutput = actTotalOutput.plotData
		self.assertEqual(expPlotOutput, actPlotOutput)

	def testTableDataMapsProperly(self):
		deltaE0Vals = [x-min(self.e0ValsA) for x in self.e0ValsA]
		expTableOutputA = [self.methStr, self.v0ValsA[0], self.b0ValsA[0], deltaE0Vals[0]]
		expTableOutputB = [self.methStr, self.v0ValsA[1], self.b0ValsA[1], deltaE0Vals[1]]
		expTableOutput = [expTableOutputA, expTableOutputB]
		actTotalOutput = self.runMapFunctOnDataA()
		actTableOutput = actTotalOutput.tableData
		self.assertEqual(expTableOutput, actTableOutput)

class TestCreateStandardInputForPlanewaveData(unittest.TestCase):

	def setUp(self):
		self.structStrs = ["hcp", "bcc"]
		self.database = mock.Mock()
		self.database.getEosFitDict.side_effect = lambda *args,**kwargs: args[0] #Means return structStr
		self.eleKey= "test_ele"
		self.eosKey="eos_key"

	def runTestCode(self):
		return tCode.createStubStandardInputObjForPlanewaveData(self.database, self.structStrs,
		                                                 eosStr=self.eosKey, eleKey=self.eleKey)
	def testExpectedDatabaseCallsMade(self):
		self.runTestCode()
		for sStr in self.structStrs:
			self.database.getEosFitDict.assert_any_call(sStr, eos=self.eosKey)

	def testExpectedLabelProduced(self):
		expLabel = labelHelp.StandardLabel(eleKey=self.eleKey, structKey="eos", methodKey="plane-wave")
		outObj = self.runTestCode()
		actLabel = outObj.label[0]
		self.assertEqual(expLabel, actLabel)

	def testExpectedOutputProduced(self):
		expOutput = [[types.SimpleNamespace(data=x) for x in self.structStrs]] #Since the database is mocked to always return structStr from calls
		inpObj = self.runTestCode()
		outObj = inpObj.createOutputObj()
		actOutput = outObj.data
		self.assertEqual(expOutput, actOutput)




