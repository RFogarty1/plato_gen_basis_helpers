

import itertools as it

import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.workflows.total_energies as tCode


class TestTotalEnergyGroupWorkflow(unittest.TestCase):

	def setUp(self):
		self.runCommsA = ["fake_comm_a","fake_comm_b"]
		self.testFlowsA = [mock.Mock(),mock.Mock()]
		self.testWeightsA = [1.0,2.0]
		self.energiesA = [20,30]
		self.createTestObjs()

	def createTestObjs(self):
		outputA = types.SimpleNamespace( energy=self.energiesA[0] )
		outputB = types.SimpleNamespace( energy=self.energiesA[1] )
		
		self.testFlowsA[0].output = [outputA]
		self.testFlowsA[1].output = [outputB]

		for obj,runComm in it.zip_longest(self.testFlowsA, self.runCommsA):
			obj.preRunShellComms = [runComm]

		self.testObjA = tCode.TotalEnergyGroupWorkflow(self.testFlowsA, self.testWeightsA)

	def testRaisesIfWrongNumberOfWeights(self):
		self.testWeightsA.append(5)
		self.assertNotEqual( len(self.testWeightsA), len(self.testFlowsA) )
		with self.assertRaises(AssertionError):
			self.createTestObjs()

	def testExpPreRunShellComms(self):
		expRunComms = self.runCommsA
		actRunComms = self.testObjA.preRunShellComms
		self.assertEqual(expRunComms,actRunComms)

	def testExpectedTotalEnergy(self):
		self.testObjA.run()
		for x in self.testFlowsA:
			x.run.assert_called_once_with()

		expEnergy = sum( [e*weight for e,weight in it.zip_longest(self.energiesA,self.testWeightsA)] )
		actEnergy = self.testObjA.output[0].energy
		self.assertTrue( len(self.testObjA.output)==1 )
		self.assertAlmostEqual(expEnergy,actEnergy)

class TestTotalEnergyWorkflow(unittest.TestCase):

	def setUp(self):
		self.runCommA = "fake-run-comms"
		self.eType = "electronicTotalE"
		self.outEnergy = 40
		self.numbAtoms = 2
		self.ePerAtom = False
		self.createTestObjs()

	def createTestObjs(self):
		self.calcObjA = mock.Mock()
		self.calcObjA.runComm = self.runCommA
		energyObj = types.SimpleNamespace( **{self.eType:self.outEnergy} )
		self.calcObjA.parsedFile = types.SimpleNamespace(energies=energyObj, numbAtoms=self.numbAtoms)
		self.testObjA = tCode.TotalEnergyWorkflow(self.calcObjA, eType=self.eType, ePerAtom=self.ePerAtom)

	def testWritesInpFilesUponInitialisation(self):
		self.calcObjA.writeFile.assert_called_once_with()

	def testExpPreRunShellComms(self):
		expRunComms = [self.runCommA]
		actRunComms = self.testObjA.preRunShellComms
		self.assertEqual(expRunComms, actRunComms)

	def testRunForTotalEnergy(self):
		self.testObjA.run()
		expEnergy = self.outEnergy
		actEnergy = self.testObjA.output[0].energy
		self.assertEqual(expEnergy,actEnergy)

	def testRunForEPerAtom(self):
		self.ePerAtom = True
		self.createTestObjs()
		self.testObjA.run()
		expEnergy = self.outEnergy / self.numbAtoms
		actEnergy = self.testObjA.output[0].energy
		self.assertAlmostEqual(expEnergy, actEnergy)

