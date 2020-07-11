
import os
import pathlib
import types

from ..shared import method_objs as baseObjs
from ..shared import geom_constraints as geoConstrainModule
from . import method_register as methodReg

import plato_pylib.parseOther.parse_castep_files as parseCastep



class CastepCalcObjFactoryStandard(baseObjs.CalcMethodFactoryBase):

	registeredKwargs = set(baseObjs.CalcMethodFactoryBase.registeredKwargs)

	registeredKwargs.add("methodStr")
	registeredKwargs.add("pseudoPotDict") #{eleKey: name_of_pseudopotential}
	registeredKwargs.add("symmetryGenerate")
	registeredKwargs.add("cutoffEnergy")

	#TODO: These should probably be moved up to the CalcMethodFactoryBase as soon as
	# they needed implementing on ANY other class
	registeredKwargs.add("runType")
	registeredKwargs.add("geomConstraints")


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
		outDict.update( self._getParamModDictBasedOnRunType() )
		return outDict

	def _createCellFileDict(self):
		outDict = dict()
		outDict.update( parseCastep.getCellGeomDictSectionFromUCell(self.geom) )
		self._updateDictWithPseudoPotStr(outDict)
		outDict["kpoint_mp_grid"] = " ".join([str(x) for x in self.kPts])
		if self.symmetryGenerate:
			outDict["symmetry_generate"] = ""
		outDict.update( self._getCellModDictBasedOnGeomConstraints() )
		return outDict

	def _updateDictWithPseudoPotStr(self,inpDict):
		if self.pseudoPotDict is None:
			return None #Castep can generate OTF pseudopotentials for all atoms in this case
		outList = list()
		for key,val in self.pseudoPotDict.items():
			currStr = "{} {}".format(key.capitalize(),val)
			outList.append(currStr)
		outList = sorted(outList)
		specPotVal = "\n".join(outList)
		inpDict.update({"species_pot":specPotVal})

	def _getParamModDictBasedOnRunType(self):
		outDict = dict()
		runStr = self.runType if self.runType is not None else "None" #Need to be able to apply .lower() to it

		if runStr.lower() == "geomOpt".lower():
			outDict["task"] = "GeometryOptimization".lower()

		return outDict

	def _getCellModDictBasedOnGeomConstraints(self):
		outDict = dict()
		if self.geomConstraints is None:
			geomConstraints = geoConstrainModule.GeomConstraints.initWithNoConstraints()
		else:
			geomConstraints = self.geomConstraints
		outDict["cell_constraints"] = getCellConstraintsStrFromGeomConstraintsObj(geomConstraints)
		eleList = [x[-1] for x in self.geom.cartCoords]
		ionicConstraints = getIonicConstraintsStrFromAtomicPosConstraintObj(geomConstraints.atomicPositionConstraints,eleList)
		if ionicConstraints is not None:
			outDict["%block ionic_constraints"] = ionicConstraints

		return outDict


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




def getCellConstraintsStrFromGeomConstraintsObj(geomConstraintsObj):
	""" Returns the cell_constraints section str for constraining cell angles/parameters
	
	Args:
		geomConstraintsObj (GeomConstraints object): Contains information on what (if anything) to constrain
			 
	Returns
		outStr (str): Contains the cell constraints in a format for castep *.cell files
 
	"""
	#Figure out the integers we need
	counter = 1
	outInts = list()
	cellConstraints = geomConstraintsObj.cellConstraints
	allConstraintVals = cellConstraints.lattParamsToFix + cellConstraints.anglesToFix
	for constrainVal in allConstraintVals:
		if constrainVal is True:
			outInts.append(0)
		else:
			outInts.append(counter)
			counter += 1

	#Convert the integers to a string
	outStrFmt = "{} {} {}\n{} {} {}" #First 3 for latt params, Second 3 for angles
	outStr = outStrFmt.format(*outInts)
	return outStr

def getIonicConstraintsStrFromAtomicPosConstraintObj(atomicPosConstraintObj, eleList):
	""" Returns the ionic_constraints section str for fixing atomic cartesian positions
	
	Args:
		atomicPosConstraintObj (AtomicPositionConstraints obj): Contains information on what (if anything) to constrain
		eleList (iter of strs): A list of the elements in the structure in order; needed since castep requires both the element symbol and its position within those elements (e.g. for Mg,Zr,Mg the atomic indices are 1,1,2)

	Returns
		 outStr (str): Contains any ionic constraints in a format for castep *.cell files
 
	"""
	if atomicPosConstraintObj.constraintsPresent is False:
		return None

	#Get all strings in a list
	idxList = getCastepIndicesFromEleList(eleList)
	constraintsList = list()
	currConstraintIdx = 1
	for constraint in atomicPosConstraintObj.atomicCartConstraints:
		if constraint.constraintsPresent:
			atomIdx = idxList[constraint.atomIdx]
			ele = eleList[constraint.atomIdx]
			if constraint.fixX:
				constraintsList.append("{} {} {} 1.0 0.0 0.0".format(currConstraintIdx, ele, atomIdx))
				currConstraintIdx += 1

			if constraint.fixY:
				constraintsList.append("{} {} {} 0.0 1.0 0.0".format(currConstraintIdx, ele, atomIdx))
				currConstraintIdx += 1

			if constraint.fixZ:
				constraintsList.append("{} {} {} 0.0 0.0 1.0".format(currConstraintIdx, ele, atomIdx))
				currConstraintIdx += 1


	#Combine the list and return
	outStr = "\n".join(constraintsList)

	return outStr

def getCastepIndicesFromEleList(eleList):
	""" Get atomic indices as used by castep
	
	Args:
		eleList (iter of Str): Each entry is one element symbol
			 
	Returns
		 idxList (iter of ints): The element indices castep-style. This means the indices start at one and refer to indices within an atomic speciies. E.g. ["X","Y","X"] would lead to [1,1,2] returned
 
	"""
	idxList = list()
	elesFound = list()
	numbInstances = list()
	for ele in eleList:
		if ele in elesFound:
			currEleIdx = elesFound.index(ele)
			numbInstances[currEleIdx] += 1
			currIdx = numbInstances[currEleIdx]
		else:
			elesFound.append(ele)
			numbInstances.append(1)
			currIdx = 1

		idxList.append(currIdx)

	return idxList


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
		outDir = os.path.dirname(self.basePath)
		pathlib.Path(outDir).mkdir(parents=True,exist_ok=True)
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
		return getParsedFileObjFromCastepOutputFile(outFilePath)

	@property
	def basePath(self):
		return self._basePath

	@basePath.setter
	def basePath(self,val):
		self._basePath = os.path.splitext(val)[0]


def getParsedFileObjFromCastepOutputFile(outFilePath):
	parsedDict = parseCastep.parseCastepOutfile(outFilePath)
	outObj = baseObjs.StandardParsedOutputFile.fromKwargDict(**parsedDict)
	outObj.unitCell.convAngToBohr()
	return outObj


