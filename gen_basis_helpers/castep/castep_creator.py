
import os
import types

from ..shared import method_objs as baseObjs
from . import method_register as methodReg

import plato_pylib.parseOther.parse_castep_files as parseCastep



class CastepCalcObjFactoryStandard(baseObjs.CalcMethodFactoryBase):

	registeredKwargs = set(baseObjs.CalcMethodFactoryBase.registeredKwargs)

	registeredKwargs.add("methodStr")
	registeredKwargs.add("pseudoPotDict") #{eleKey: name_of_pseudopotential}
	registeredKwargs.add("symmetryGenerate")
	registeredKwargs.add("cutoffEnergy")

	#Key function
	def _createFromSelf(self):
		paramDict = self.paramFileDict
		cellDict = self.cellFileDict
		basePath = self.basePath
		outObj = CastepCalcObj(basePath, paramDict, cellDict)
		return outObj

	def _createParamFileDict(self):
		outDict = methodReg.createParamDictFromMethodStr(self.methodStr)
		outDict["cut_off_energy"] = str(self.cutoffEnergy)
		return outDict

	def _createCellFileDict(self):
		outDict = dict()
		outDict.update( parseCastep.getCellGeomDictSectionFromUCell(self.geom) )
		self._updateDictWithPseudoPotStr(outDict)
		outDict["kpoint_mp_grid"] = " ".join([str(x) for x in self.kPts])
		if self.symmetryGenerate:
			outDict["symmetry_generate"] = ""
		return outDict

	def _updateDictWithPseudoPotStr(self,inpDict):
		if self.pseudoPotDict is None:
			return None #Castep can generate OTF pseudopotentials for all atoms in this case
		outList = list()
		for key,val in self.pseudoPotDict.items():
			currStr = "{} {}".format(key.capitalize(),val)
			outList.append(currStr)
		outList = sorted(outList)
		specPotVal = "\n".join(outList) + "\n"
		inpDict.update({"species_pot":specPotVal})


	@property
	def paramFileDict(self):
		""" Read-only dict that gets written to the .param file. Note modiying the returned dict has no effect on what gets written out
		"""
		return self._createParamFileDict()

	@property
	def cellFileDict(self):
		""" Read-only dict that gets written to the .cell file. Note that modifying the returned dict has no effect on what gets written out
		"""
		return self._createCellFileDict()

	@property
	def basePath(self):
		""" Read-only path to the base file (can be with or without extension; doesnt matter)
		"""
		return os.path.join(self.workFolder,self.fileName)


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

