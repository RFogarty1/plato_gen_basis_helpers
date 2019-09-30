#!/usr/bin/python3

import copy
import os
import types
import unittest
import unittest.mock as mock

import numpy as np

import gen_basis_helpers.job_utils.dos_objs as tCode

class TestDosRunnerPlato(unittest.TestCase):

	def setUp(self):
		self.mockPlatoCalcObj = createMockPlatoCalcObjA()
		self.smearWidth = 0.5
		self.stepSize = 0.1
		self.eleKey = "element"
		self.structKey = "structure"
		self.methodKey = "method"
		self.runDosGenerating = True
		self.runEnergy = True
		self.createTestObj()

	def createTestObj(self):
		self.testLabel = tCode.DosLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)
		self.testObj = tCode.DosRunnerPlato(self.mockPlatoCalcObj, self.smearWidth, self.stepSize, self.testLabel,
		                                    runDosGenerating=self.runDosGenerating, runEnergy=self.runEnergy)

	def testExpectedEnergyComms(self):
		expectedRunComms = ["fake_run_comm_a"]
		actualRunComm = self.testObj.singlePointEnergyComms
		self.assertEqual(expectedRunComms, actualRunComm)

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.dosHelp.getDosRunComm_plato")
	def testExpectedDoSComms(self, getDosRunCommMock):
		getDosRunCommMock.return_value = "fake_dos_comm_a"
		expectedOutput = ["fake_dos_comm_a"]
		actualOutput = self.testObj.dosComms
		getDosRunCommMock.assert_called_once_with(os.path.join("heres","a","path.occ"), self.smearWidth, self.stepSize)
		self.assertEqual(expectedOutput, actualOutput)

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.DosRunnerPlato.singlePointEnergyComms",new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_utils.dos_objs.jobRun.executeRunCommsParralel")
	def testRunEnergyComms(self, mockedRunParralel, mockedComms):
		fakeComms = ["fake_energy_comms"]
		mockedComms.return_value = fakeComms
		self.testObj.runSinglePointEnergyCalcs()
		self.mockPlatoCalcObj.writeFile.assert_called_once_with()
		mockedRunParralel.assert_called_once_with(fakeComms,1)

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.DosRunnerPlato.singlePointEnergyComms",new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_utils.dos_objs.jobRun.executeRunCommsParralel")
	def testRunEnergyCommsWhenSetToFalse(self, mockedRunParralel, mockedComms):
		self.createTestObj()
		self.assertTrue(self.testObj.runEnergy)
		self.testObj.runEnergy = False
		self.testObj.runSinglePointEnergyCalcs()
		mockedRunParralel.assert_not_called()


	@mock.patch("gen_basis_helpers.job_utils.dos_objs.DosRunnerPlato.dosComms",new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_utils.dos_objs.jobRun.executeRunCommsParralel")
	def testRunDosComms(self, mockedRunParralel, mockedComms):
		fakeComm = ["fake_dos_comm"]
		mockedComms.return_value = fakeComm
		self.testObj.runDosGeneratingCalcs()
		mockedRunParralel.assert_called_once_with(fakeComm,1) 

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.DosRunnerPlato.dosComms",new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_utils.dos_objs.jobRun.executeRunCommsParralel")
	def testRunDosCommsWhenSetToFalse(self, mockedExecuteRunComms, mockedComms):
		self.createTestObj()
		self.assertTrue(self.testObj.runDosGenerating)
		self.testObj.runDosGenerating = False
		self.testObj.runDosGeneratingCalcs()
		mockedExecuteRunComms.assert_not_called()

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.DosAnalyserStandard")
	@mock.patch("gen_basis_helpers.job_utils.dos_objs.dosHelp.getDosPlotData")
	def testCreateAnalyser(self, mockedGetDosData, mockedAnalyserClass):
		fakeDosData = "fake_dos_data"
		fakeAnalyserObject = "fake_analyser_object"
		fakeEFermi = 7

		mockedGetDosData.return_value = {"dosdata":fakeDosData, "efermi":fakeEFermi}
		mockedAnalyserClass.return_value = fakeAnalyserObject
		outObj = self.testObj.createAnalyser()

		mockedAnalyserClass.assert_called_with(fakeDosData, fakeEFermi, self.testLabel)
		self.assertEqual(outObj,fakeAnalyserObject)


def createMockPlatoCalcObjA():
	outMock = mock.Mock()
	outMock.getRunComm.return_value = "fake_run_comm_a"
	outMock.filePath = os.path.join("heres","a","path.in")
	return outMock


class TestDosLabel(unittest.TestCase):


	def setUp(self):
		self.structKey = "struct"
		self.eleKey = "ele"
		self.methodKey = "method"
		self.createTestObj()

	def createTestObj(self):
		self.testObj = tCode.DosLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)

	def testInitFailsWhenMissingMethodKey(self):
		with self.assertRaises(ValueError):
			tCode.DosLabel(self.eleKey,self.structKey)

	def testEqualityMethodTwoEqual(self):
		testObjA = copy.deepcopy(self.testObj)
		self.createTestObj()
		self.assertFalse(testObjA is self.testObj)
		self.assertEqual(testObjA,self.testObj)

	def testEqualityMethodTwoNotEqual(self):
		testObjA = copy.deepcopy(self.testObj)
		self.methodKey += "_hi"
		self.createTestObj()
		self.assertFalse( testObjA is self.testObj )
		self.assertFalse(testObjA == self.testObj)


class TestDosRunnerComposite(unittest.TestCase):

	def setUp(self):
		self.runnerStubA = createRunnerStubA()
		self.runnerStubB = createRunnerStubB()
		self.createTestObj()

	def createTestObj(self):
		self.testObj = tCode.DosRunnerComposite( [self.runnerStubA, self.runnerStubB] )

	def testLabelProperty(self):
		expLabel = ["fake_labelA","fake_labelB", "fake_labelC"]
		actLabel = self.testObj.label
		self.assertEqual(expLabel, actLabel)

	def testInitFailsForDuplicateLabels(self):
		self.runnerStubA.label = self.runnerStubB.label
		with self.assertRaises(ValueError):
			self.createTestObj()

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.DosAnalyserComposite")
	def testCreateAnalyser(self, analyserMock):
		expInpObjs = ["fake_analyser_a", "fake_analyser_b"]
		self.testObj.createAnalyser()
		analyserMock.assert_called_once_with(expInpObjs)

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.jobRun.executeRunCommsParralel")
	def testRunDosComms(self, mockRunner):
		nCores=4
		self.testObj.runDosGeneratingCalcs(nCores)
		expDosComms = ["fake_dos_a","fake_dos_b","fake_dos_c","fake_dos_d"]
		mockRunner.assert_called_once_with(expDosComms, nCores)

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.jobRun.executeRunCommsParralel")
	def testRunEnergyComms(self,mockRunner):
		nCores=3
		self.testObj.runSinglePointEnergyCalcs(nCores)
		expEnergyComms = ["fake_e_a", "fake_e_b", "fake_e_c", "fake_e_d"]
		mockRunner.assert_called_once_with(expEnergyComms, nCores)

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.jobRun.executeRunCommsParralel")
	def testRunDosWhenDosGeneratingFalse(self, mockRunner):
		self.assertTrue(self.testObj.runDosGenerating)
		self.testObj.runDosGenerating = False
		self.testObj.runDosGeneratingCalcs()
		mockRunner.assert_not_called()

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.jobRun.executeRunCommsParralel")
	def testRunEnergyCalcsWhenRunEnergyFalse(self, mockRunner):
		self.assertTrue(self.testObj.runEnergy)
		self.testObj.runEnergy = False
		self.testObj.runSinglePointEnergyCalcs()
		mockRunner.assert_not_called()


	
def createRunnerStubA():
	runnerDictA = {"label":["fake_labelA"],
	               "createAnalyser": lambda : "fake_analyser_a",
	               "dosComms":["fake_dos_a","fake_dos_b"],
	               "singlePointEnergyComms":["fake_e_a", "fake_e_b"],
	               "runEnergy":True,
	               "runDosGenerating":True,
	               "writeFiles": lambda: None}
	return types.SimpleNamespace(**runnerDictA)

#This one should mimic a composite
def createRunnerStubB():
	runnerDictB = {"label":["fake_labelB", "fake_labelC"],
	               "createAnalyser": lambda : "fake_analyser_b",
	               "dosComms":["fake_dos_c","fake_dos_d"],
	               "singlePointEnergyComms":["fake_e_c", "fake_e_d"],
	               "runEnergy": True,
	               "runDosGenerating":True,
	               "writeFiles": lambda: None}
	return types.SimpleNamespace(**runnerDictB)



class TestDosAnalyserStandard(unittest.TestCase):

	def setUp(self):
		self.dosData = [ [1,2], [2,4] ]
		self.eFermi = 5.4
		self.eleKey = "Zr"
		self.methodKey = "method"
		self.structKey = "struct"
		self.createTestObj()

	def createTestObj(self):
		label = tCode.DosLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)
		self.testObj = tCode.DosAnalyserStandard(self.dosData, self.eFermi, label)

	def testGetObjectsWithComponentsAllMatch(self):
		components = [self.eleKey, self.methodKey, self.structKey]
		outObjs = self.testObj.getObjectsWithComponents(components)
		self.assertTrue( len(outObjs)==1 )

	def testGetObjectsWithComponentsEmptyList(self):
		components = []
		outObjs = self.testObj.getObjectsWithComponents(components)
		self.assertTrue( len(outObjs)==1 )
	
	def testGetObjectWithComponentsCaseInsensitiveMatch(self):
		components = [x.upper() for x in [self.eleKey, self.methodKey, self.structKey]]
		components[0] = components[0].lower()
		outObjs = self.testObj.getObjectsWithComponents(components, caseSensitive=False)
		self.assertTrue( len(outObjs)==1 )

	def testGetObjectWithComponentsNotMatching(self):
		components = [self.eleKey, self.methodKey,"this_doesnt_match"]
		outObjs = self.testObj.getObjectsWithComponents(components)
		self.assertTrue( len(outObjs)==0 )

	def testAttachRefDataCorrectComponents(self):
		components = [self.eleKey]
		fakeRefData = np.array( [ [2,3],[3,4] ] )
		self.testObj.attachRefData(fakeRefData,components)
		self.assertTrue( np.allclose( fakeRefData, self.testObj.refData ) )


class TestDosAnalyserComposite(unittest.TestCase):
	
	def setUp(self):
		self.analyserA = createAnalyserStubA()
		self.analyserB = createAnalyserStubB()
		self.analyserC = createAnalyserStubC()
		self.allAnalysers = [self.analyserA, self.analyserB, self.analyserC]
		for x in self.allAnalysers:
			x.plotData = mock.Mock()
			x.plotData.return_value = ["handle"]
		self.createTestObj()

	def createTestObj(self):
		self.testObj = tCode.DosAnalyserComposite( self.allAnalysers )

	def testGetObjectsWithComponents(self):
		fakeComponents = ["irrelevant"]
		actObjs = self.testObj.getObjectsWithComponents(fakeComponents)
		expObjs = ["fake_obj_a", "fake_obj_b","fake_obj_c"]
		[x.getObjectsWithComponents.assert_called_once_with(fakeComponents, caseSensitive=True) for x in self.allAnalysers]
		self.assertEqual(expObjs, actObjs)

	def testAttachRefDataMultipleMatches(self):
		fakeComponents = ["irrelevant"]
		refDataToAttach = "not_important"
		self.testObj.attachRefData(refDataToAttach, fakeComponents, errorIfNoMatches=True)
		for x in self.allAnalysers:
			x.attachRefData.assert_called_once_with(refDataToAttach,fakeComponents, errorIfNoMatches=False, caseSensitiveComponents=True)

	def testAttachRefDataErrorsWhenNoMatches(self):
		fakeComponents = ["irrelevant"]
		refDataToAttach = "unimportant"
		for x in self.allAnalysers:
			x.getObjectsWithComponents.return_value = list()
		self.createTestObj()
		with self.assertRaises(ValueError):
			self.testObj.attachRefData(refDataToAttach, fakeComponents, errorIfNoMatches=True)

	def testPlotData(self):
		fakeKwargs = {"xlim":[4,10], "ylim":[3,5]}
		expVals = ["handle", "handle", "handle"]
		actVals = self.testObj.plotData(**fakeKwargs)
		self.assertEqual(expVals, actVals)
		for x in self.allAnalysers:
			x.plotData.assert_called_once_with( xlim=[4,10], ylim=[3,5] )


def createAnalyserStubA():
	getObjsWithCompsFunct = mock.Mock()
	attachRefData = mock.Mock()
	getObjsWithCompsFunct.return_value = ["fake_obj_a"]
	outDict = {"getObjectsWithComponents": getObjsWithCompsFunct,
	           "attachRefData": attachRefData}
	outObj = types.SimpleNamespace(**outDict)
	return outObj

def createAnalyserStubB():
	getObjsWithCompsFunct = mock.Mock()
	attachRefData = mock.Mock()
	getObjsWithCompsFunct.return_value = []
	outDict = {"getObjectsWithComponents": getObjsWithCompsFunct,
	           "attachRefData": attachRefData}
	outObj = types.SimpleNamespace(**outDict)
	return outObj

def createAnalyserStubC():
	getObjsWithCompsFunct = mock.Mock()
	attachRefData = mock.Mock()
	getObjsWithCompsFunct.return_value = ["fake_obj_b","fake_obj_c"]
	outDict = {"getObjectsWithComponents": getObjsWithCompsFunct,
	           "attachRefData": attachRefData}
	outObj = types.SimpleNamespace(**outDict)
	return outObj

if __name__ == '__main__':
	unittest.main()

