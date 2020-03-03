
import os
import types

import plato_pylib.plato.mod_plato_inp_files as modPlatoInp
import plato_pylib.plato.parse_plato_out_files as parsePlatoOut

from ..shared import method_objs as baseObjs
from ..shared import calc_methods as calcMeth


class PlatoCalcObjFactoryStandard(baseObjs.CalcMethodFactoryBase):

	registeredKwargs = set(baseObjs.CalcMethodFactoryBase.registeredKwargs)
	registeredKwargs.add("methodStr")
	registeredKwargs.add("dataSet") #Plato relative path
	registeredKwargs.add("gridVals") #Format varies based on methodStr; for dft fft grid is used, fot dft2 atom centred is used

	#Key function
	def _createFromSelf(self):
		startObj = self._getMethodObj()
		self._modifyMethodObj(startObj)
		outStrDict = self._getStrDictFromMethodObj(startObj)
		runCommFunct = startObj.runCommFunction
		outBasePath = os.path.join(self.workFolder, self.fileName)
		return PlatoCalcObj(outBasePath, outStrDict, runCommFunct)

	#If not set here they default to None; These defaults will always be overriden if the explicit kwarg is passed to init (or create) to set them
	def _setDefaultInitAttrs(self):
		self.fileName = "plato_file.in"


	#Not the type of object we actually want in the end; wrong interface and incomplete
	def _getMethodObj(self):
		return calcMeth.createPlatoMethodObj(self.methodStr)

	def _modifyMethodObj(self, methObj):
		if self.kPts is not None:
			methObj.kpts = self.kPts
		if self.gridVals is not None:
			methObj.integGrid = self.gridVals
		if self.dataSet is not None:
			methObj.dataSet = self.dataSet

	def _getStrDictFromMethodObj(self, methObj):
		return methObj.getStrDictWithStruct(self.geom)


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



	


