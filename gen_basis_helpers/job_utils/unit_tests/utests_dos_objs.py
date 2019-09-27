#!/usr/bin/python3

import copy
import os
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.job_utils.dos_objs as tCode

class TestDosRunnerPlato(unittest.TestCase):

	def setUp(self):
		self.mockPlatoCalcObj = createMockPlatoCalcObjA()
		self.smearWidth = 0.5
		self.stepSize = 0.1
		self.eleKey = "element"
		self.structKey = "structure"
		self.methodKey = "method"
		self.createTestObj()

	def createTestObj(self):
		self.testLabel = tCode.DosLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)
		self.testObj = tCode.DosRunnerPlato(self.mockPlatoCalcObj, self.smearWidth, self.stepSize, self.testLabel)

	def testExpectedEnergyComms(self):
		expectedRunComms = ["fake_run_comm_a"]
		actualRunComm = self.testObj.singlePointEnergyComms
		self.assertEqual(expectedRunComms, actualRunComm)

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.dosHelp.getDosPlotData")
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

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.DosRunnerPlato.dosComms",new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_utils.dos_objs.jobRun.executeRunCommsParralel")
	def testRunDosComms(self, mockedRunParralel, mockedComms):
		fakeComm = ["fake_dos_comm"]
		mockedComms.return_value = fakeComm
		self.testObj.runDosGeneratingCalcs()
		mockedRunParralel.assert_called_once_with(fakeComm,1) 

	@mock.patch("gen_basis_helpers.job_utils.dos_objs.DosAnalyserPlato")
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
		pass
	@unittest.skip("")
	def testInitFailsForDuplicateLabels(self):
		self.assertTrue(False)



def createRunnerCompositeFromStubsA():
	runnerDictA = {"label": ["fake_labelA"]}
	runnerDictB = {"label": ["fake_labelB"]}
	#TODO: Use these stubs as input to a composite constructor; maybe need to make this part of the class actually (so i can set and mess with runnerDictA properties)
	return None


if __name__ == '__main__':
	unittest.main()

