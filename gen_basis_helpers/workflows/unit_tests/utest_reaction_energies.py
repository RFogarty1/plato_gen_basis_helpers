
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

		self.reactantComponents = [self.reactantEnergy+2, self.reactantEnergy-2]
		self.productComponents = [self.productEnergy]

		outputReactant = types.SimpleNamespace(energy=self.reactantEnergy, componentEnergies=self.reactantComponents) 
		outputProduct = types.SimpleNamespace(energy=self.productEnergy, componentEnergies=self.productComponents)
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

	def testExpReactantAndProductContribs(self):
		self.testObjA.run()
		self.assertEqual( self.reactantComponents, self.testObjA.output[0].reactantEnergies )
		self.assertEqual( self.productComponents, self.testObjA.output[0].productEnergies )


class TestCreateWorkflowFromEnergiesAndStoichiometries(unittest.TestCase):

	def setUp(self):
		self.reactantEnergyA = 5
		self.reactantEnergyB = -3
		self.productEnergyA = 4
		self.reactantStoicA = 2
		self.reactantStoicB = 3
		self.productStoicA = 3
		self.createTestObjs()

	def createTestObjs(self):
		self.reactantEnergiesA =[self.reactantEnergyA, self.reactantEnergyB] 
		self.productEnergiesA = [self.productEnergyA]
		self.reactantStoicsA = [self.reactantStoicA, self.reactantStoicB]
		self.productStoicsA = [self.productStoicA]
		self.testObjA = self._runFunction()

	def _runFunction(self):
		inpArgs = (self.reactantEnergiesA, self.reactantStoicsA, self.productEnergiesA, self.productStoicsA)
		return tCode.createReactionEnergyWorkflowFromEnergiesAndStoics( *inpArgs )

	def testExpectedReactionEnergyGiven(self):
		self.testObjA.run()
		expReactEnergy =( (self.productStoicA*self.productEnergyA) - 
		                 ( (self.reactantStoicA*self.reactantEnergyA) + (self.reactantStoicB*self.reactantEnergyB) ) )
		actReactEnergy = self.testObjA.output[0].energy
		self.assertAlmostEqual(expReactEnergy, actReactEnergy)

	def testAssertRaisedIfStoicsAndEnergiesDiffLengths(self):
		self.reactantStoicsA.append(5)
		self.assertTrue( len(self.reactantEnergiesA) != len(self.reactantStoicsA) )
		with self.assertRaises(AssertionError):
			self._runFunction()

