
import itertools as it
import types

import plato_fit_integrals.core.workflow_coordinator as wFlow

class ConvergerWorkflowTemplate(wFlow.WorkFlowBase):
	""" Workflow for testing convergence of some property with respect to input convergence values.
	"""

	def __init__(self, calcObjs, convVals, mapFunction, namespaceAttrs):
		""" Initializer for template converger workflow. In this case the functions that require defining/overwriting are passed as explicit arguments. The alternative (more standard) way would be to use inheritence to overwrite these functions.
		
		Args:
			calcObjs: (iter of CalcMethod objects) Each represents a calculation we need to carry out
			convVals: (iter of numbers) List of values we're testing for convergence, for example a number representing the tightness of an integration grid. Note the order must correspond to the order of calcObjs
			mapFunction: This is the run() method of the class; hence the interface is mapFunction(workflow). This takes the workflow, after all calculations have been run, and figures out what to put in workflow.output [i.e. how to go from completed calculations to some desired output]
			namespaceAttrs: (iter) List of strings defining the attributes that mapFunction puts on the self.output object
		
		"""
		self.calcObjs = list(calcObjs)
		self.convVals = convVals
		self.mapFunction = mapFunction
		self._namespaceAttrs = namespaceAttrs

		for x in self.calcObjs:
			x.writeFile()

	@property
	def preRunShellComms(self):
		outComms = list()
		for x in self.calcObjs:
			outComms.append(x.runComm)
		return outComms

	@property
	def namespaceAttrs(self):
		return self._namespaceAttrs

	#Setting to None could lead to bugs if using with a co-ordinator, which wont catch the case when we reuse a workfolder. OTOH it shouldnt
	#be dificult for users to make sure file paths dont overlap
	@property
	def workFolder(self):
		return None

	def run(self):
		self.mapFunction(self)


class GridConvergenceEnergyWorkflow(ConvergerWorkflowTemplate):
	""" Workflow for testing convergence of energy with respect to grid spacing for a structure """


	def __init__(self, calcObjs, convVals, energyStr="electronicTotalE", ePerAtom=True):
		""" Initializer for a workflow to check convergence of energy

		Args:
			calcObjs: (iter of CalcMethod objects) Each represents a calculation we need to carry out
			convVals: (iter of numbers) List of values we're testing for convergence, for example a number representing the tightness of an integration grid. Note the order must correspond to the order of calcObjs
			energyStr: (str, optional) The type of energy to request from plato_pylib Energies object. (options correspond to properties on that object)
			ePerAtom: (Bool, optional) Whether to report energy per atom (instead of total). Default is True

		"""
		self.calcObjs = calcObjs
		self.convVals = convVals
		self.energyStr = energyStr
		self.output = types.SimpleNamespace()
		self.ePerAtom = bool(ePerAtom)

		for x in self.calcObjs:
			x.writeFile()


	@property
	def namespaceAttrs(self):
		return ["convResults"]

	def run(self):
		energies = [getattr(x.parsedFile.energies,self.energyStr) for x in self.calcObjs]
		if self.ePerAtom:
			nAtoms = [x.parsedFile.numbAtoms for x in self.calcObjs]
			energies = [e/nAtoms for e,nAtoms in it.zip_longest(energies,nAtoms)]

		outVals = [(conv,e) for conv,e in it.zip_longest(self.convVals,energies)]
		self.output.convResults = outVals



