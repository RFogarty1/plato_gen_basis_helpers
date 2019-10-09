import abc

#TODO: These actually have an optDict property, which is pretty important (though doesnt have to be a property, usually just attr)This should be clearer on the interface definition (i.e. here)
class CalcMethod(abc.ABC):

	@abc.abstractmethod
	def writeFile(self):
		""" Should write an input file. Path should be set such that its consistent with outFilePath """ 
		pass

	@property
	@abc.abstractmethod
	def outFilePath(self):
		""" Path to the output file """
		pass

	@property
	@abc.abstractmethod
	def nCores(self):
		""" Number of cores to use for THIS calculation. """
		pass


	@property
	@abc.abstractmethod
	def runComm(self):
		""" Command to run the job. Should be possible to pass to subprocess.call """
		pass

	@property
	@abc.abstractmethod
	def parsedFile(self):
		""" Minimal implementation is a Namespace containing information parsed from the output file.
			The simplest way to do this starting from a dict(inputDict) is:

			from types import SimpleNameSpace
			nameSpace = SimpleNameSpace(**inputDict)

			This means the keys will accesble as nameSpace.attrName. The return value of this function can
			therefore be any class. I may add a few standard attribute names later (e.g. unitCell, energies)
		"""
		pass


