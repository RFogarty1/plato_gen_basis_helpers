
import types
from . import base_flow as baseFlow

class ParsedFileWorkflow(baseFlow.BaseLabelledWorkflow):
	""" Basic workflow to run a single job and get the parsed output file from it

	"""

	def __init__(self, calcObj):
		""" Initializer
		
		Args:
			calcObj (CalcMethod object): Object contains methods to parse the file and run the job
		"""
		self.calcObj = calcObj
		self._writeInpFiles()
		self._output =  types.SimpleNamespace( **{k:None for k in self.namespaceAttrs[0]} )

	def _writeInpFiles(self):
		self.calcObj.writeFile()

	def run(self):
		self._output.parsedFile = self.calcObj.parsedFile

	@property
	def output(self):
		return [self._output]

	@property
	def namespaceAttrs(self):
		return [ ["parsedFile"] ]

	@property
	def preRunShellComms(self):
		return [self.calcObj.runComm]
