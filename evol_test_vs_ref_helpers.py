#!/usr/bin/python3

import itertools
import os
import shutil
import sys

import plato_pylib.plato.mod_plato_inp_files as platoInp
import plato_pylib.parseOther.parse_castep_files as parseCastep
import plato_pylib.plato.parse_plato_out_files as platoOut
import plato_pylib.utils.job_running_functs as jobRun

import dft_file_helpers as dftHelpers

''' Functions to help run test of energy vs volume data vs a reference () for plato basis sets '''

DFT_ONLY_OPT_KEYS = [x.lower() for x in ["model", "optimizemesh", "fftGridSpacing", "ACGGridSpacing", "paralleliseFFT", "SpinMixLevels", "diagonalisationMethod",
 "DiagonalisationWorkSize",  "MaxBond",  "DensityRWFlag",  "writeAtomDensityFlag",  "WavefunctionFlag", 
 "HyperfineFlag", "ACGEwaldParameters",  "ACGMeshReduction",  "ACGMinPartitionWt",  "ACGAngularMeshType",  "DiameterNLV",
 "IntegralPartitionFlag",  "DensityFit",  "NumericVna",  "SpinMixScheme", "MixThresHold", 
    "MixMetric",  "SplitDensityFlag"] ]

def writeTb2OrTb1FileFromDict(outPath, inpDict):
	#Get default options
	optDict = loadDefaultTb2OptDict()
	optDict = {k.lower():v for k,v in optDict.items()}
	#Update any options requested by user
	dCopy = {k.lower():v for k,v in inpDict.items()}
	optDict.update(dCopy)

	optDict = getPlatoStrDictFromOptDict_tb1OrTb2(optDict)

	#Write the file
	platoInp.writePlatoOutFileFromDict(outPath, optDict)


def getPlatoStrDictFromOptDict_tb1OrTb2(inpDict):
	return platoInp.getStrDictFromOptDict(inpDict, "dft2")


def loadDefaultTb1OptDict():
	return platoInp.getDefOptDict("tb1")

def loadDefaultTb2OptDict():
	return platoInp.getDefOptDict("dft2")

def getPlatoOptDictFromCastepOutFile(casFile, mode=None):
	modeToDefDict = {"tb2": loadDefaultTb2OptDict, "tb1": loadDefaultTb1OptDict}
	if mode is None:
		mode == "tb2"

	defaultDict = modeToDefDict[mode]()
	return _getUpdatedOptDictFromCastepOutFile(casFile,defaultDict)



def getPlatoTb2OptDictFromCastepOutFile(casFile):
	outDict = loadDefaultTb2OptDict()
	return _getUpdatedOptDictFromCastepOutFile(casFile, outDict)



def _getUpdatedOptDictFromCastepOutFile(casFile,outDict):
	parsedFile = parseCastep.parseCastepOutfile(casFile)

	outDict["natom"] = parsedFile["numbAtoms"]

	#Get parameters from unitcell
	uCell = parsedFile["unitCell"]
	uCell.convAngToBohr()
	outDict["cellsize"] = uCell.getLattParamsList()
	lattVects = uCell.getLattVects()
	for idx,lVect in enumerate(lattVects):
		lattVects[idx] = [x/outDict["cellsize"][idx] for x in lVect]
	outDict["cellvec"] = lattVects
	outDict["format"] = 0

	#The fractional co-ordinations - Not really checked up on this yet
	outDict["atoms"] = list()
	for fCoord in uCell.fractCoords:
		outDict["atoms"].append( " ".join([str(x) for x in fCoord]) )

	return outDict


#At first this will only be for angular grids, but will define an interface
#that can easily be extended to radial grids
def getGridSpacingVsEnergyFolder(folder, mode=None, angIdx=None):
	#Could feasibly add a step to check whether grid is angular or radial here, and delegate
	#each to a separate function
	return getAtomGridSpacingVsEnergyFolder(folder, mode, angIdx)

def getAtomGridSpacingVsEnergyFolder(folder,mode:"str-either angular or radial", angIdx:"two angular values-this says whether to take 0th or 1st"= 0):
	if mode!="angular" and mode!="radial":
		raise ValueError("{} is an invalid value for mode".format(mode))	

	#Get list of all out files
	outFileList = [ os.path.join(folder,x) for x in os.listdir(folder) if x.endswith('.out') ]
	baseFilePathList =  [ os.path.splitext(x)[0] for x in outFileList ]

	#Grab spacing from the input files
	gridSpacing = [_grabAtomGridSpacingFromInpFile(x+'.in', mode) for x in baseFilePathList]

	#Grab cohesive energies from the output files
	parsedOutFiles = [platoOut.parsePlatoOutFile(x) for x in outFileList]
	cohEnergies = [x["energies"].electronicCohesiveE for x in parsedOutFiles]

	#combine + return
	outData = list()
	for spacing,energy in itertools.zip_longest(gridSpacing, cohEnergies):
		outData.append( (spacing,energy) )

	return outData

def _grabAtomGridSpacingFromInpFile(filePath, mode, angIdx=0):
	parsedFile = platoInp.tokenizePlatoInpFile(filePath)
	gridPart = parsedFile["integralmesh"].split("\n")

	if (mode == "angular"):
		gridVal = int( gridPart[1].split()[angIdx+1] )
	elif (mode == "radial"):
		gridVal = int( gridPart[1].split()[0] )
	else:
		raise ValueError("{} is an invalid mode for _grabAtomGridSpacingFromInpFile".format(mode))
	return gridVal




class IntegGridConverger:
	def __init__(self, templateFile, mode:"angular, radial or kpt", convVariables, **kwargs):
		kwargs = {k.lower():v for k,v in kwargs.items()}
		self.templateFile = templateFile
		self.mode = mode
		self.convVariables = list(convVariables)
		self.atomGridFixedVal = kwargs.get("atomGridFixedVal".lower(),None)

		self.runComm = "dft" if self.mode=="fftgrid" else "dft2" 

	def runConvCalcs(self, nCores=1):
		allInpFiles = self._createConvFiles()
		allRunComms = jobRun.pathListToPlatoRunComms(allInpFiles,self.runComm)
		jobRun.executeRunCommsParralel(allRunComms, nCores)

	def getVarVsCohesive(self):
		parsedOutFiles = [platoOut.parsePlatoOutFile(x.replace(".in",".out")) for x in self._createConvFiles(createFiles=False)]
		if self.mode=="angular" or self.mode=="radial" or self.mode=="kpt":
			cohEnergies = [x["energies"].electronicCohesiveE for x in parsedOutFiles]
		elif self.mode=="fftgrid":
			cohEnergies = [x["energies"].electronicTotalE for x in parsedOutFiles]

		if self.mode=="angular" or self.mode == "radial": 
			convVars = self.convVariables
		elif self.mode=="kpt":
			convVars = [min(x) for x in self.convVariables]
		elif self.mode=="fftgrid":
			convVars = self.convVariables
		else:
			raise ValueError("{} is an invalid option for self.mode".format(self.mode))

		outData = list()
		for spacing,energy in itertools.zip_longest(convVars, cohEnergies):
			outData.append( (spacing,energy) )
	
		return outData
			

	def _createConvFiles(self,createFiles=True):
		if self.mode=="angular" or self.mode=="radial":
			outFiles = _createTb2AtomCentGridConvFromTemplate(self.templateFile, self.mode, self.convVariables, self.atomGridFixedVal, createFiles=createFiles)
		elif self.mode=="kpt":
			outFiles = self._createKptConvFilesFromTemplate(createFiles=createFiles)
		elif self.mode=="fftgrid":
			outFiles = self._createFFTConvFilesFromTemplate(createFiles=createFiles)
		else:
			raise ValueError("{} is an invalid option for self.mode".format(self.mode))

		return outFiles

	def _createKptConvFilesFromTemplate(self,createFiles=True):
		outFileList = list()
		outFormStr = "-1\n{} {} {}"
		baseFolder = os.path.split(self.templateFile)[0]

		fileNameFormat = "{}_{}_{}.in"
		for varVal in self.convVariables:
			fName = fileNameFormat.format(*varVal)
			fPath = os.path.join(baseFolder,fName)
			outFileList.append(fPath)
			if createFiles:
				shutil.copy2(self.templateFile, fPath)
				platoInp.modPlatoInpOptions(fPath, {"blochstates": outFormStr.format(*varVal)})
		return outFileList

	def _createFFTConvFilesFromTemplate(self,createFiles=True):
		outFileList = list()
		fileNameFormat = "grid_{}.in"
		baseFolder = os.path.split(self.templateFile)[0]

		for varVal in self.convVariables:
			fName = fileNameFormat.format( str(varVal).replace(".","pt") )
			fPath = os.path.join(baseFolder,fName)
			outFileList.append(fPath)
			if createFiles:
				shutil.copy2(self.templateFile,fPath)
				platoInp.modPlatoInpOptions(fPath, {"fftGridSpacing": str(varVal)})
		return outFileList


def _createTb2AtomCentGridConvFromTemplate(templatePath, mode, variableVals, fixedVal, createFiles=True):
	if mode != "angular" and mode!="radial":
		raise ValueError("{} is an invalid mode for createGridConvFilesAtomCentreGrid".format(mode))

	#Set the string for grid spacing in the inp file
	if mode == "radial":
		gridStr = "1\n{var} {fix} {fix}"
	else:
		gridStr = "1\n{fix} {var} {var}"

	#Create the files
	outFileList = list()
	fileNameFormat = "{}_{}.in"
	baseFolder = os.path.split(templatePath)[0]
	for varVal in variableVals:
		fName = fileNameFormat.format(varVal,mode)
		fPath = os.path.join(baseFolder,fName)
		outFileList.append(fPath)
		if createFiles:
			shutil.copy2(templatePath, fPath)
			platoInp.modPlatoInpOptions( fPath, {"integralmesh":gridStr.format(var=varVal,fix=fixedVal)} )
	return outFileList




