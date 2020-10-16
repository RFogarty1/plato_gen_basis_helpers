
from . import base_flow as baseFlow

import types


class GeomOptWorkflow(baseFlow.BaseLabelledWorkflow):
	""" Workflow for carrying out a geometry optimisation (using whatever optimisers are available in the program used"""

	def __init__(self, calcObj):
		""" Initializer
		
		Args:
			calcObj: (CalcMethod object). Object represents most aspects of this calculation; contains methods to write and parse the file + to run the job. NOTE: the workflow has no way to check this is set to actually do a geometry optimisation
				 
		"""

		self.calcObj = calcObj
		self._output = types.SimpleNamespace(**{k:None for k in self.namespaceAttrs})
		self._writeInpFiles()

	def _writeInpFiles(self):
		self.calcObj.writeFile()

	@property
	def preRunShellComms(self):
		return [self.calcObj.runComm]

	@property
	def namespaceAttrs(self):
		return ["geom", "energy", "parsedFile"]

	@property
	def output(self):
		return [self._output]

	def run(self):
		parsedFile = self.calcObj.parsedFile
		self._output.geom = parsedFile.unitCell
		self._output.parsedFile = parsedFile
		self._output.energy = getattr(parsedFile.energies, "electronicTotalE")

