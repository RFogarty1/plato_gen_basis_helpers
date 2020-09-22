
import plato_pylib.shared.custom_errors as custErrors
import types
from . import base_flow as baseFlow


class ParsedFileWorkflow(baseFlow.BaseLabelledWorkflow):
	""" Basic workflow to run a single job and get the parsed output file from it

	"""

	def __init__(self, calcObj, catchParserErrors=False):
		""" Initializer
		
		Args:
			calcObj (CalcMethod object): Object contains methods to parse the file and run the job
		"""
		self.calcObj = calcObj
		self.catchParserErrors = catchParserErrors
		self._writeInpFiles()
		self._output =  types.SimpleNamespace( **{k:None for k in self.namespaceAttrs[0]} )

	def _writeInpFiles(self):
		self.calcObj.writeFile()

	def run(self):
		if self.catchParserErrors:
			try:
				self._run()		
			except custErrors.PlatoPylibParseFileError:
				pass
		else:
			self._run()
			
	def _run(self):
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
