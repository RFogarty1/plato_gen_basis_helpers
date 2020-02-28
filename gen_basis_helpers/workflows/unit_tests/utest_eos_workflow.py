
import itertools as it
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.method_objs as methodObjs
import gen_basis_helpers.workflows.eos_workflow as tCode


class DudCalcObj(methodObjs.CalcMethod):

	def __init__(self, runComm, energyObj, numbAtoms, volume):
		self._runComm = runComm
		self._eObj = energyObj
		self._numbAtoms = numbAtoms
		self._volume = volume

	@property
	def runComm(self):
		return self._runComm

	@property
	def outFilePath(self):
		raise NotImplementedError("")

	@property
	def nCores(self):
		raise NotImplementedError("")

	@property
	def parsedFile(self):
		uCell = types.SimpleNamespace(volume=self._volume)
		outStub = types.SimpleNamespace( energies=self._eObj, numbAtoms=self._numbAtoms, unitCell=uCell )
		return outStub

	def writeFile(self):
		pass

class TestEosWorkflow(unittest.TestCase):

	def setUp(self):
		self.runCommsA = ["fake_comm_a", "fake_comm_b"] #1 per object
		self.energiesA = [40,20]
		self.volumesA = [20,40]
		self.nAtomsA = [2,4]
		self.calcObjs = list()
		self.fitFunction = _stubFitEosFunctionA
		self.labelA = None
		self.createTestObjs()

	def createTestObjs(self):
		fakeEnergyObjsA = [types.SimpleNamespace(electronicTotalE=x) for x in self.energiesA]
		calcObjsA = [DudCalcObj(rComm,eObj,nAtoms,vol) for rComm,eObj,nAtoms,vol in it.zip_longest(self.runCommsA,fakeEnergyObjsA,self.nAtomsA,self.volumesA)]
		self.testObjA = tCode.EosWorkflow(calcObjsA, self.fitFunction, self.labelA)

	def testRunCommsAsExpected(self):
		expRunComms = self.runCommsA
		actRunComms = self.testObjA.preRunShellComms
		self.assertEqual(expRunComms, actRunComms)

	#Mocks cant handle division, hence use of stubs instead here
	def testOutputContainsResultsInExpectedPlace(self):
		self.testObjA.run()
		expOutput = types.SimpleNamespace( data=types.SimpleNamespace(testAttr="test_strA") )
		actOutput = self.testObjA.output[0]
		self.assertTrue( len(self.testObjA.output)==1 )
		self.assertEqual(expOutput, actOutput)


def _stubFitEosFunctionA(volumes,energy):
	return types.SimpleNamespace(testAttr="test_strA")


class TestStandardFitEosFunction(unittest.TestCase):
	""" We only need to test that it correctly passes its values on to the plato_pylib interface """

	def setUp(self):
		self.testVolumesA = [10,30]
		self.testEnergiesA = [20,50]
		self.eosStr = "fake_str"
		self.maxFev = 2
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.StandardEosFitFunction( self.eosStr, self.maxFev )

	@mock.patch("gen_basis_helpers.workflows.eos_workflow.fitEos")
	def testExpectedValuesPassedToInterface(self, mockedFitEos):
		expOutput = "fake_output"
		mockedFitEos.getBulkModFromVolsAndEnergiesBohrAndEvUnits.side_effect = lambda *args,**kwargs:expOutput
		actOutput = self.testObjA(self.testVolumesA, self.testEnergiesA)
		expArgs = (self.testVolumesA, self.testEnergiesA)
		expKwargs = {"eosModel":self.eosStr, "maxFev":self.maxFev}
		mockedFitEos.getBulkModFromVolsAndEnergiesBohrAndEvUnits.assert_called_once_with(*expArgs, **expKwargs)
		self.assertEqual(expOutput,actOutput)


