

import types
from . import base_flow as baseFlow



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

	@property
	def output(self):
		return [self._output]

	@property	
	def namespaceAttrs(self):
		return ["energy"]

