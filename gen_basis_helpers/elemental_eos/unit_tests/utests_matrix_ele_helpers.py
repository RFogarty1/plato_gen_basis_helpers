#!/usr/bin/python3

import copy
import itertools as it
import types
import unittest
import unittest.mock as mock

import numpy as np

import gen_basis_helpers.elemental_eos.matrix_ele_helpers as tCode

class TestGetMatrixElementsFromWorkFlow(unittest.TestCase):

	def setUp(self):
		self.testWorkFlow = createStubWorkFlowA()

	def runTestFunct(self):
		tCode.volVsDiagMatrixElementsFromWorkFlow(self.testWorkFlow)

	@mock.patch("gen_basis_helpers.elemental_eos.matrix_ele_helpers.matHelp.getVolVsOnSiteDiagTerms")
	def testExpectedOutputGenerated(self, mockGetDiagEles):
		fakeReturnVal = ["fake_val_a", "fake_val_b"]
		mockGetDiagEles.return_value = fakeReturnVal
		self.runTestFunct()
		mockGetDiagEles.assert_called_once_with(["fake_path_a", "fake_path_b"], volPerAtom=True)
		expOutput = types.SimpleNamespace(volVsDiagOnSite=fakeReturnVal)
		actOutput = self.testWorkFlow.output
		self.assertTrue(expOutput, actOutput)


def createStubWorkFlowA():
	outDict = {"_inpFilePaths":["fake_path_a", "fake_path_b"]}
	return types.SimpleNamespace(**outDict)


class TestReplaceRunMethodInWorkFlow(unittest.TestCase):

	def setUp(self):
		self.stubWorkFlow = types.SimpleNamespace(prop="fake_run_ret_val")
		self.stubFunct = lambda x:  x.prop

	def testForMockFunctA(self):
		tCode.replaceRunMethodInWorkFlow(self.stubWorkFlow, self.stubFunct)
		expOutVal = "fake_run_ret_val"
		actOutVal = self.stubWorkFlow.run()
		self.assertEqual(expOutVal, actOutVal)


class TestDiagResultsComposite(unittest.TestCase):

	def setUp(self):
		self.stubObjA = createMatrixDiagResultsStubA()
		self.stubObjB = createMatrixDiagResultsStubB()
		self.createObj()

	def createObj(self):
		self.testObj = tCode.EleEosMatrixDiagElementsComposite([self.stubObjA,self.stubObjB])

	def testInitFailsWhenDuplicatedLabels(self):
		self.stubObjA.label[0] = copy.deepcopy((self.stubObjB.label[0]))
		with self.assertRaises(ValueError):
			self.createObj()

	def testGetPlotData(self):
		expData = [ [1,2], [3,4], [5,6] ]
		actData = self.testObj.getPlotData()
		self.assertEqual(expData, actData)


def createMatrixDiagResultsStubA():
	label = [tCode.StandardLabel(eleKey="Zr", structKey="hcp", methodKey="dft")]
	data = [[1,2]]

	outDict = { "label":label,
	            "getPlotData": lambda: data }
	return types.SimpleNamespace(**outDict)


def createMatrixDiagResultsStubB():
	labelA = tCode.StandardLabel(eleKey="Mg", structKey="hcp", methodKey="tb1")
	labelB = tCode.StandardLabel(eleKey="Hf", structKey="hcp", methodKey="tb1")
	dataA = [3,4]
	dataB = [5,6]

	outDict = { "label": [labelA, labelB], 
	            "getPlotData": lambda: [dataA,dataB] }
	return types.SimpleNamespace(**outDict)


class TestCreateMatrixDiagResultsFromWorkFlow(unittest.TestCase):

	def setUp(self):
		diagDataVolA = [ [1,2], [3,4,5 ], [5,5] ] #1 entry per atom, lists are the diag matrix elements for that atom
		diagDataVolB = [ [6,7], [8,9,10], [6,6] ]
		volA, volB = 50, 60
		self.fakeData = [ [volA, diagDataVolA], [volB, diagDataVolB] ]
		self.fakeOutput = types.SimpleNamespace( volVsDiagOnSite= self.fakeData )
		self.methodLabel = "fake_method"
		self.structLabel = "fake_struct"
		self.eleLabel = "fake_ele"
		self.createTestWorkFlowStub()


	def createTestWorkFlowStub(self):
		self.testWorkFlowStub = types.SimpleNamespace( output=self.fakeOutput,
		                                               element=self.eleLabel,
		                                               structKey=self.structLabel,
		                                               methodKey=self.methodLabel ) 

	def runTestFunct(self):
		self.testObj =  tCode.createMatrixEleAnalyserFromWorkFlow(self.testWorkFlowStub)

	def testExpDataPassedOver(self):
		expOutData = [ np.array( [ [50,1,2,3,4,5 ,5,5],
		                           [60,6,7,8,9,10,6,6] ] ) ]
		expLabel = tCode.StandardLabel(eleKey=self.eleLabel, methodKey=self.methodLabel, structKey=self.structLabel)

		self.runTestFunct()
		actOutData = self.testObj.getPlotData()

		self.assertTrue( np.allclose(expOutData, actOutData) )
		self.assertEqual( expLabel, self.testObj.label[0] )


class TestShellMapper(unittest.TestCase):

	def setUp(self):
		self.sShellStart = [8,9]
		self.pShellStart = [5]
		self.dShellStart = [0]
		self.createTestObj()

	def createTestObj(self):
		self.testObj = tCode.ShellMapper( sStart=self.sShellStart, pStart=self.pShellStart, dStart=self.dShellStart )

	def testGetDiagIndicesNthShellOfType_validInput(self):
		testArgs = [("s", 0),
		            ("s", 1),
		            ("p", 0),
		            ("d", 0)]

		testAnswers = [  [8],
		                 [9],
		                 [5,6,7],
		                 [0,1,2,3,4] ]

		actAnswers = [self.testObj.getDiagIndicesNthShellOfType(*inpArgs) for inpArgs in testArgs]
		for exp,act in it.zip_longest(testAnswers, actAnswers):
			self.assertEqual(exp,act)

	def testEqMethod_copiedObj(self):
		objA = self.testObj
		objB = copy.deepcopy(self.testObj)
		self.assertEqual( objA, objB )

	def testEqMethod_slightlyDiffObjs(self):
		objA = copy.deepcopy(self.testObj)
		self.pShellStart[0] += 1
		self.createTestObj()
		objB = copy.deepcopy(self.testObj)
		self.assertNotEqual(objA,objB)

	def testInitFromShellOrderings(self):
		testShellOrderings = ["d","p","s","s"]
		self.sShellStart, self.pShellStart, self.dShellStart = [8,9], [5], [0]
		self.createTestObj()
		expObj = self.testObj
		actObj = tCode.ShellMapper.fromShellOrderings(testShellOrderings)
		self.assertEqual(expObj, actObj)

	@mock.patch("gen_basis_helpers.elemental_eos.matrix_ele_helpers.parseTbint.parseAdtFile")
	def testCreateFromAdtFile(self, mockedParser):
		fakeFilePath = "random_string"
		mockedParser.return_value =  createStubParsedAdtFile()
		expObj = self.testObj
		actObj = tCode.createShellMapperFromAdtFile(fakeFilePath)
		self.assertEqual(expObj,actObj)

	def testNumbOrbs(self):
		expNumbOrbs = 10
		actNumbOrbs = self.testObj.numbOrbs
		self.assertEqual(expNumbOrbs, actNumbOrbs)


def createStubParsedAdtFile():
	#keys are shell indices, vals are ang-momenta
	shellToAngMomMap = { 1:1, 0:2, 2:0, 3:0}
	return { "shellToAngMom": shellToAngMomMap }





class TestDataPresenter(unittest.TestCase):

	def setUp(self):
		self.analyserA = mock.Mock()
		self.returnedAnalyser = mock.Mock()
		self.analyserA.getObjectsWithComponents.return_value = [self.returnedAnalyser]
		self.shellMapperA = tCode.ShellMapper.fromShellOrderings(["s","p","d"])
		self.shellMapNamespace = types.SimpleNamespace(zr=self.shellMapperA)
		self.dataPlotter = mock.Mock()
		self.createTestObj()

	def createTestObj(self):
		self.testObj = tCode.MatEleEosDataPresenter(self.analyserA, self.shellMapNamespace, dataPlotter=self.dataPlotter)

	def testGetShellMappedArray(self):
		testOrdering = ["p","s"]
		#Actual ordering set by self.shellMapperA and is vol | sVal |pValA |pValB | pValC | dValA...
		volRowA = [50, 1, 2, 3, 4, 5] #s,p1,p2,p3,d
		volRowB = [60, 6, 7, 8, 9, 10]
		inpArray = np.array( [volRowA, volRowB] )
		expArray = np.array( [ [50, 2, 3, 4, 1],
		                       [60, 7, 8, 9, 6] ] )

		actArray = self.testObj._getShellMappedVersionOfArray( inpArray, testOrdering, self.shellMapperA ) 

		self.assertTrue( np.allclose(expArray, actArray) )


	def testGetPlotDataSingleMethod(self):
		testMethods = ["testMA"]
		testStructLabel, testEleKey = "hcp", "zr"
		self.returnedAnalyser.getPlotData.return_value = ["fake_np_array"]
		self.createTestObj()

		expDict = {"testMA":"fake_np_array"}
		actDict = self.testObj.getPlotDataWithoutMapping(testMethods, testStructLabel, testEleKey)

		self.analyserA.getObjectsWithComponents.assert_called_once_with( ["testMA", testStructLabel, testEleKey],caseSensitive=False )
		self.assertEqual(expDict, actDict)

	@mock.patch("gen_basis_helpers.elemental_eos.matrix_ele_helpers.MatEleEosDataPresenter._getShellMappedVersionOfArray")
	@mock.patch("gen_basis_helpers.elemental_eos.matrix_ele_helpers.MatEleEosDataPresenter.getPlotDataWithoutMapping")
	def testGetMappedPlotDataOneEleAndStruct(self, mockedGetData, mockedMapFunct):
		testEle, testStruct = "zr", "irrelevent"
		testShellOrdering = ["s","p"] #Doesnt actually matter due to mocking out the map function
		testSubMethod = True
		testMethods = ["methA","methB"]
		testRefMethod = "methB"
		dataDict = {"methA": np.array(([1,2],[2,4])),
		            "methB": np.array(([1,3],[2,3]))}
		mockedGetData.return_value = dataDict
		mockedMapFunct.side_effect = lambda x,*args: x
		expOutList = [ np.array(([1,-1],[2,1])),
		               np.array(([1, 0],[2,0])) ]

		actOutList = self.testObj.getMappedPlotDataOneEleAndStruct(testMethods, testStruct, testEle, testShellOrdering, subRefData=True, refMethod=testRefMethod)


		mockedGetData.assert_called_once_with( testMethods, testStruct, testEle )
		self.assertTrue( mockedMapFunct.call_count==2 ) #verifying args is tricky in this case, you can an error from comparing np arrays i think
#		mockedMapFunct.assert_has_calls( [mock.call(dataDict["methA"], testShellOrdering, self.shellMapperA)] )
#		mockedMapFunct.assert_has_calls( [mock.call(dataDict["methB"], testShellOrdering, self.shellMapperA)] )
		for exp,act in it.zip_longest(expOutList,actOutList):
			self.assertTrue( np.allclose( exp, act ) )









if __name__ == '__main__':
	unittest.main()
