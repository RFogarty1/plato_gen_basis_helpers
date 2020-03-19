
import copy
import types
import unittest
import unittest.mock as mock


import gen_basis_helpers.workflows.reaction_energies as tCode


class TestReactionEnergyWorkflows(unittest.TestCase):

	def setUp(self):
		self.runCommsReactants = ["fake_comm_reactants"]
		self.runCommsProducts = ["fake_comm_products"]
		self.productEnergy = 20
		self.reactantEnergy = 14
		self.createTestObjs()

	def createTestObjs(self):
		self.reactants, self.products = mock.Mock(), mock.Mock()
		self.reactants.preRunShellComms = self.runCommsReactants
		self.products.preRunShellComms = self.runCommsProducts 

		outputReactant = types.SimpleNamespace(energy=self.reactantEnergy) 
		outputProduct = types.SimpleNamespace(energy=self.productEnergy)
		self.reactants.output = [outputReactant]
		self.products.output = [outputProduct]

		self.testObjA = tCode.ReactionEnergyWorkflow(self.reactants, self.products)


	def testRunComms(self):
		expComms = self.runCommsReactants + self.runCommsProducts
		actComms = self.testObjA.preRunShellComms
		self.assertEqual(expComms,actComms)

	def testRunCallsReactantAndProductRunMethods(self):
		self.testObjA.run()
		self.reactants.run.assert_called_once_with()
		self.products.run.assert_called_once_with()

	def testExpReactionEnergyReturned(self):
		self.testObjA.run()
		expEnergy = self.productEnergy - self.reactantEnergy
		actEnergy = self.testObjA.output[0].energy
		self.assertAlmostEqual(expEnergy,actEnergy)
		self.assertTrue( len(self.testObjA.output)==1 )
