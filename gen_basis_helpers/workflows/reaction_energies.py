
import itertools as it
import types

from . import base_flow as baseFlow
from . import total_energies as totEnergyFlow

class ReactionEnergyWorkflow(baseFlow.BaseLabelledWorkflow):
	"""Workflow for calculating reaction energies

	"""

	def __init__(self, reactantWorkflow, productWorkflow):
		""" Initializer
		
		Args:
			reactantWorkflow: (TotalEnergyWorkflowBase object)
			productWorkflow: (TotalEnergyWorkflowBase object)
				 
		"""
		self.reactantWorkflow = reactantWorkflow
		self.productWorkflow = productWorkflow
		self._output = types.SimpleNamespace(**{k:None for k in self.namespaceAttrs})


	@property
	def preRunShellComms(self):
		return self.reactantWorkflow.preRunShellComms + self.productWorkflow.preRunShellComms

	def run(self):
		self.reactantWorkflow.run()
		self.productWorkflow.run()
		assert len(self.reactantWorkflow.output)==1
		assert len(self.productWorkflow.output)==1
		reactionEnergy = self.productWorkflow.output[0].energy - self.reactantWorkflow.output[0].energy
		self._output.energy = reactionEnergy

		#Record the (stoichiometry weighted) contributions from individual reactants and products
		self._output.reactantEnergies = self.reactantWorkflow.output[0].componentEnergies
		self._output.productEnergies = self.productWorkflow.output[0].componentEnergies

	@property
	def output(self):
		return [self._output]

	@property	
	def namespaceAttrs(self):
		return ["energy", "reactantEnergies", "productEnergies"]


def createReactionEnergyWorkflowFromEnergiesAndStoics(reactantEnergies, reactantStoics, productEnergies, productStoics):
	""" Creates a ReactionEnergyWorkflow when given reactant/product energies and stoichiometries. This is non-trivial since the workflow components are designed such that QM codes are used to calculate the individual energies
	
	Args:
		reactantEnergies: (float iter) Energies for each reactant
		reactantStoics: (float iter) Stoichiometries for each reactant, in same order as reactant energies list 
		productEnergies: (float iter) Energies for each product
		productStoics: (float iter) Stoichiometries for each product, in same order as produce energies list

	Returns
		outFlow: (ReactionEnergyWorkflow object) Once run method has been called, this contains an output property which holds the total reaction energy as well as contribution from each component

	Raises:
		AssertionError: If lengths of energy and stoichiometry iters do not match
 
	"""
	reactantObjs = [_getTotalEnergyStubFromEnergy(e) for e in reactantEnergies]
	productObjs  =  [_getTotalEnergyStubFromEnergy(e) for e in productEnergies]
	reactants = totEnergyFlow.TotalEnergyGroupWorkflow( reactantObjs, reactantStoics )
	products = totEnergyFlow.TotalEnergyGroupWorkflow( productObjs, productStoics )
	return ReactionEnergyWorkflow(reactants, products)



def _getTotalEnergyStubFromEnergy(energy):
	outputObj = types.SimpleNamespace( energy=energy, compnentEnergies=[energy] )
	return types.SimpleNamespace( output=[outputObj], run=lambda:None )

