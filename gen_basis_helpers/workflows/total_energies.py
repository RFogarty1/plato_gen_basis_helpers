
import itertools as it
import types

from . import base_flow as baseFlow


class TotalEnergyWorkflowBase(baseFlow.BaseLabelledWorkflow):
	"""Base object for total energy workflows; this can return total energy of a single structure or a group of structures (e.g. all reactants)

	"""
	@property	
	def namespaceAttrs(self):
		return ["energy", "componentEnergies"]

	@property
	def output(self):
		""" Needs to return a len-1 iter containing a types.SimpleNamespace object with namespaceAttrs defined
		"""
		raise NotImplementedError("")

	#Usually the base class behaviour is to return an empty list rather than throw an error
	@property
	def preRunShellComms(self):
		raise NotImplementedError("")


class TotalEnergyGroupWorkflow(TotalEnergyWorkflowBase):
	"""A group of TotalEnergyWorkflow objects; this is intended to represent all reactants or all products for a reaction. Calculates a weighted sum (to deal with stoichiometry issues) of the input total energy workflows

	"""

	def __init__(self, totEnergyFlows, weights):
		""" Initializer
		
		Args:
			totEnergyFlows: (iter of TotalEnergyWorkflow objects) Each should represent energy of one structure
			weights: (iter of floats) Values to multiply each total energy workflow by. These should be the stoichiometries if usign the workflow for calculating reaction energies
		
		"""
		assert len(weights)==len(totEnergyFlows), "{} weights given for {} total energy workflows; these two numbers need to be equal".format( len(weights), len(totEnergyFlows) )
		self.totEnergyFlows = totEnergyFlows
		self.weights = weights
		self._output = types.SimpleNamespace(**{k:None for k in self.namespaceAttrs})

	@property
	def preRunShellComms(self):
		outList = list()
		for x in self.totEnergyFlows:
			outList.extend(x.preRunShellComms)
		return outList


	def run(self):
		for x in self.totEnergyFlows:
			x.run()

		assert [len(x.output)==1 for x in self.totEnergyFlows]
		energies = [x.output[0].energy for x in self.totEnergyFlows]
		componentEnergies = [e*w for e,w in it.zip_longest(energies,self.weights)]
		weightedEnergy = sum(componentEnergies)
		self._output.energy = weightedEnergy
		self._output.componentEnergies = componentEnergies

	@property
	def output(self):
		return [self._output]



class TotalEnergyWorkflow(TotalEnergyWorkflowBase):
	"""Workflow for getting the total energy of a single structure. This can be electronic/free energy etc, or it can be the energy per atom of a bulk structure. 

	"""

	def __init__(self, calcObj, label=None, eType="electronicTotalE", ePerAtom=False):
		""" Initializer
		
		Args:
			calcObj: (CalcMethod object) Object represents most aspects of this calculation; contains methods to write and parse the file + to run the job
			ePerAtom: (Bool, optional) If true return the energy per atom, rather than total energy.
 	 
		"""

		self.calcObj = calcObj
		self.eType = eType
		self.ePerAtom = ePerAtom
		self._writeInpFiles()

		self._output = types.SimpleNamespace(**{k:None for k in self.namespaceAttrs})

	def _writeInpFiles(self):
		self.calcObj.writeFile()

	@property
	def preRunShellComms(self):
		return [self.calcObj.runComm]

	@property
	def output(self):
		return [self._output]

	def run(self):
		parsedFile = self.calcObj.parsedFile
		totalEnergy = getattr(parsedFile.energies, self.eType)
		if self.ePerAtom:
			outEnergy = totalEnergy / parsedFile.numbAtoms
		else:
			outEnergy = totalEnergy

		#Component energy is really mainly for composite workflows, but needed here to keep the interface the same
		self._output.energy = outEnergy
		self._output.componentEnergies = [outEnergy]
