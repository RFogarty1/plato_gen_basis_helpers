import abc

from . import creator_resetable_kwargs as baseCreator

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

			The output should implement the method_objs.ParsedFile interface as much as possible
		"""
		pass



class ParsedFile():

	@property
	def energies(self):
		""" Should return a plato_pylib Energies object """
		raise NotImplementedError("")

	@property
	def numbAtoms(self):
		""" The number of atoms in the calculation """
		raise NotImplementedError("")

	@property
	def unitCell(self):
		""" A plato_pylib UnitCell object; units should generally be in bohr """
		raise NotImplementedError("")


class StandardParsedOutputFile(ParsedFile):

	def __init__(self, energies=None, numbAtoms=None, unitCell=None):
		""" Initializer
		
		Args:
			energies: plato_pylib Energies object
			numbAtoms: The number of atoms in the calculation
			unitCell: plato_pylib UnitCell object; units should generally be in bohr

		"""
		self._energies = energies
		self._numbAtoms = numbAtoms
		self._unitCell = unitCell

	@classmethod
	def fromKwargDict(cls, **kwargs):
		allowedKwargs = ["energies", "numbAtoms", "unitCell"]
		outKwargs = {k:kwargs[k] for k in kwargs if k in allowedKwargs}
		return cls(**outKwargs)

	@property
	def energies(self):
		return self._energies

	@property
	def numbAtoms(self):
		return self._numbAtoms

	@property
	def unitCell(self):
		return self._unitCell



class CalcMethodFactoryBase(baseCreator.CreatorWithResetableKwargsTemplate):

	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs) #Will probably always be empty
	registeredKwargs.add("workFolder")
	registeredKwargs.add("fileName")
	registeredKwargs.add("kPts")
	registeredKwargs.add("geom")



class StubCalcMethodFromParsedFileObject(CalcMethod):

	def __init__(self, parsedFile):
		self._parsedFile = parsedFile

	def writeFile(self):
		pass

	@property
	def outFilePath(self):
		return ""

	@property
	def nCores(self):
		return 1

	@property
	def runComm(self):
		return list()

	@property
	def parsedFile(self):
		return self._parsedFile

