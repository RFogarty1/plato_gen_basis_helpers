
import os
import types

from ..shared import method_objs as baseObjs

import plato_pylib.parseOther.parse_castep_files as parseCastep

class CastepCalcObj(baseObjs.CalcMethod):

	def __init__(self, basePath, paramFileDict, cellFileDict):
		""" Initializer
		
		Args:
			basePath: (str) Path, without extension, to write files to. If extension it will simply be removed (so not a problem)
			paramFileDict: (dict) Keys/vals are both strings; keys are castep *.param file keywords while vals are the string representation of the vlaues
			cellFileDict: (dict) Keys/vals both strings. Keys are castep *.cell file keywords while vals are their str representations. Preceding the key with %block will ensure its treated as a multiline block; but this generally should NOT be neccesary (exceptions are certain single line blocks)
				
		"""
		self.basePath = basePath
		self.paramFileDict = paramFileDict
		self.cellFileDict = cellFileDict

	def writeFile(self):
		cellPath, paramPath = self.basePath+".cell", self.basePath+".param"
		parseCastep.writeCastepParamFileFromDict(paramPath,self.paramFileDict)
		parseCastep.writeCastepCellFileFromTokens(cellPath, self.cellFileDict)

	@property
	def outFilePath(self):
		raise NotImplementedError("")

	@property
	def nCores(self):
		raise NotImplementedError("")

	@property
	def runComm(self):
		return list() #Havnt got castep on local  + generally youd want to run on a hpc

	@property
	def parsedFile(self):
		outFilePath = self.basePath + ".castep" #Only parse the basic output file for now
		parsedDict = parseCastep.parseCastepOutfile(outFilePath)
		outObj = types.SimpleNamespace(**parsedDict)
		outObj.unitCell.convAngToBohr()
		return outObj

	@property
	def basePath(self):
		return self._basePath

	@basePath.setter
	def basePath(self,val):
		self._basePath = os.path.splitext(val)[0]

