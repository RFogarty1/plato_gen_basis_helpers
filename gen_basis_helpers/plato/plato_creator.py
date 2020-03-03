
import os
import types

import plato_pylib.plato.mod_plato_inp_files as modPlatoInp
import plato_pylib.plato.parse_plato_out_files as parsePlatoOut

from ..shared import method_objs as baseObjs



class PlatoCalcObjFactoryStandard(baseObjs.CalcMethodFactoryBase):

	registeredKwargs = set(baseObjs.CalcMethodFactoryBase.registeredKwargs)

	def _createFromSelf(self):
		raise NotImplementedError("")


#TypeError: Can't instantiate abstract class PlatoCalcObj with abstract methods nCores, outFilePath, parsedFile, runComm, writeFile

class PlatoCalcObj(baseObjs.CalcMethod):
	""" Plato version of the CalcMethod interface. This object is used to encapsulate the writing/running/parsing of plato jobs (similar to the Calculator in ASE)

	"""

	def __init__(self, basePath, strDict, pathToRunCommFunction):
		""" Initializer
		
		Args:
			basePath: (str) Path to the input file without file extension [though if you include the extension it will just be removed at initiation time]
			strDict: (dict) Dictionary of Token:Value for plato input file e.g. ["blochstates":"-1\n10 10 10"]. Modifying this directly is messy and not recommended; this dictionary is expected to be fixed for the lifetime of this object
			pathToRunCommFunction: (function, f(str)) Takes the plato filePath and returns a bash command to run plato on it
	
		"""
		self.basePath = os.path.splitext( basePath )[0]
		self.strDict = strDict
		self.pathToRunCommFunction = pathToRunCommFunction

	@property
	def nCores(self):
		raise NotImplementedError("")

	@property
	def outFilePath(self):
		return self.basePath + ".out"

	@property
	def parsedFile(self):
		outDict = parsePlatoOut.parsePlatoOutFile_energiesInEv(self.outFilePath)
		return types.SimpleNamespace(**outDict)

	@property
	def runComm(self):
		return self.pathToRunCommFunction(self._getInpFilePath())

	def writeFile(self):
		inpFilePath = self._getInpFilePath()
		modPlatoInp.writePlatoOutFileFromDict(inpFilePath,self.strDict)

	def _getInpFilePath(self):
		return self.basePath + ".in"



	


