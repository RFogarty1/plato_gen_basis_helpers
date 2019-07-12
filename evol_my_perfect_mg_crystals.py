

''' Very specific code with the purpose of helping to run/analyse calculations to get bulk modulii/eqm volumes '''
''' etc. for Mg perfect crystals '''

import copy
import itertools as it
import os
import sys

import numpy as np

import calc_methods as calcMethods
import ref_data as refData

sys.path.append('/media/ssd1/rf614/Work/usr_scripts/coding/Plato_Analysis_Lib_Functions')
import fit_bulk_mod as fitBMod



def getStandardMultiCrystsDictFromMethodStrsAndDataSetFolder(methodStrs:list, dataSetFolder:"str to tell plato where to find model"):
	startFolder = getStandardStartFolderBulkModCalcsBasedOnCurrDir()
	kPts = getStandardKPtDict()
	structDict = getStandardStructDictHcpFccBcc()
	gridDict = getStandardGridDictsForMethods(methodStrs)
	startDict = {"kpts": kPts,
	             "startFolder": startFolder,
	             "structDict": structDict,
	             "dataFolder": dataSetFolder}
	optsDict = createAllEvolOptsDictFromStartDictAndMethodStrs(startDict, methodStrs, gridDict)
	allMultiCrystsDict = {k:v.getMultiCrystalObj() for k,v in optsDict.items()}

	return allMultiCrystsDict


#TODO: Remove all the duplicated code from getStandardMultiCrystsDictFromMethodStrsAndDataSetFolder
def getLargeVolMultiCrystsDictsFromMethodStrsAndDataSetFolder(methodStrs:list, dataSetFolder:"str to tell plato where to find model"):
	startFolder = getStandardStartFolderBulkModCalcsBasedOnCurrDir()
	startFolder = os.path.join(startFolder,"larger_volumes")
	kPts = getStandardKPtDict()

	#Get structures then set volumes to be larger
	maxVolume, volStep = 250,5
	structDict = getStandardStructDictHcpFccBcc()
	for key in structDict:
		volVals = _getVolumeValues(maxVolume,volStep,len(structDict[key]))
		for x,volVal in it.zip_longest( structDict[key],volVals  ):
			nAtoms = len(x.fractCoords)
			x.volume = volVal*nAtoms



	gridDict = getStandardGridDictsForMethods(methodStrs)
	startDict = {"kpts": kPts,
	             "startFolder": startFolder,
	             "structDict": structDict,
	             "dataFolder": dataSetFolder}


	optsDict = createAllEvolOptsDictFromStartDictAndMethodStrs(startDict, methodStrs, gridDict)
	allMultiCrystsDict = {k:v.getMultiCrystalObj() for k,v in optsDict.items()}

	return allMultiCrystsDict

def _getVolumeValues(maxVolume,volStep,numbStructs):
	outVols = list()
	for sIdx in range(numbStructs):
		currVol = maxVolume - (volStep*sIdx)
		outVols.append(currVol)
	return outVols

def getStandardGridDictsForMethods(methodStrs):
	outDict = dict()
	for x in methodStrs:
		if x.startswith("dft2_"):
			outDict[x] = getStandardGridValDictDft2()
		elif x.startswith("dft_"):
			outDict[x] = getStandardFFTGridVals()
		else:
			outDict[x] = None
	return outDict

def getStandardStartFolderBulkModCalcsBasedOnCurrDir():
	return os.path.abspath( os.path.join("work_folder", "bmod_calcs","att1") )


def getStandardStructDictHcpFccBcc():
	structDict  = {"hcp": refData.getUCellsForBulkModCalcs("hcp"),
	               "bcc": refData.getUCellsForBulkModCalcs("bcc"),
	               "fcc": refData.getUCellsForBulkModCalcs("fcc")}
	return structDict


def getStandardGridValDictDft2():
	return {"fcc": [50,50,50], "bcc": [60,50,50], "hcp": [50,50,50] } 


def getStandardFFTGridVals():
	return {key:0.15 for key in ["fcc","bcc","hcp"]} 

def getStandardKPtDict():
	return {"fcc": [20,20,20], "bcc": [20,20,20], "hcp": [20,20,12] }

def getRefBModMockedMultiCrystal():
    bmodDict = _getRefBModDict()
    mockObj = EvolMultiCrystType(None)
    mockObj.fittedEos = bmodDict
    return mockObj

def _getRefBModDict():
    refEosDict = {"hcp": refData.getPlaneWaveEosFitDict("hcp"),
                    "bcc": refData.getPlaneWaveEosFitDict("bcc"),
                    "fcc": refData.getPlaneWaveEosFitDict("fcc")}
    return refEosDict


def createAllEvolOptsDictFromStartDictAndMethodStrs(startDict, methodStrs, gridValsDict):
	""" Creates a set of options(held in EvolCalcOptions obj) for each method in methodStrs
	
	Args:
		startDict: dict with keys "kptsDict", "startFolder", "structDict", "datafolder" (case insensitive).
		           startFolder is the root dir to make jobs.
		           structDict links a set of keys (e.g. hcp) to lists of UnitCell structs (plato_pylib objs)
		           kpts links each structure to a list of kpts to use
		           dataFolder is a string telling plato where to find the model (dataset in plato inp files)
		methodStrs: A list of strings, each indicating a method defined in calc_methods.py (or externally with 
		            same interface)
		gridValsDict: Each key is a method string, and it links to its own dict. This dict contains the grid values
		              to use for every structure (e.g. "fcc":[50,40,40] may be a key/val)

	Returns:
		outDict: Keys are methodStrs, values are EvolCalcOptions objects. These encapsulate the calculation method/structures to use.
		         Though really their main purpose is probably just to have getMultiCrystalObj() called on them
	
	Raises:
		None
	"""
	optDict = {k.lower():v for k,v in startDict.items()}
	outDict = dict()
	for key in methodStrs:
		outDict[key] = EvolCalcOptions(key, optDict["kpts"], optDict["startfolder"], optDict["structdict"], optDict["datafolder"], gridValsDict[key])
	
	return outDict


class EvolCalcOptions():    
    def __init__(self, methodStr, kpts, startFolder, structDict, dataFolder:"str telling plato where to find the model", gridVals):
        self.kpts = kpts
        self.methodStr = methodStr
        self.startFolder = startFolder
        self.structDict = structDict
        self.gridVals = gridVals
        self.dataFolder = dataFolder
        #Create the method obj
        self.platoMethodDict = self._createPlatoMethodDict()
    
    def _createPlatoMethodDict(self):
        outDict = dict()
        for key in self.structDict.keys():
            currMethodObj = calcMethods.createPlatoMethodObj(self.methodStr)
            currMethodObj.dataSet = self.dataFolder
            currMethodObj.kpts = self.kpts[key]
            if self.gridVals is not None:
                currMethodObj.integGrid = self.gridVals[key]
            outDict[key] = currMethodObj
        return outDict
    
    def getMultiCrystalObj(self):
        multiCrystDict = dict()
        for key in self.structDict.keys():
            outFolder = os.path.join(self.startFolder, self.methodStr, key)
            currObjs = _getCalcObjsListFromStructListAndFolderAndMethodObj(self.structDict[key], outFolder, self.platoMethodDict[key])
            multiCrystDict[key] = EvolSingleCrystType(currObjs)
        return EvolMultiCrystType(multiCrystDict)
    
class EvolMultiCrystType():
    def __init__(self, calcObjDict):
        self.calcObjDict = calcObjDict
    
    def writeFiles(self):
        for key in self.calcObjDict:
            self.calcObjDict[key].writeFiles()
    
    def getRunComms(self):
        allRunComms = list()
        for key in self.calcObjDict:
            currRunComms = self.calcObjDict[key].getRunComms()
            allRunComms.extend(currRunComms)
        return allRunComms
    
    
    def fitEos(self, eos):
        bmodDict = dict()
        for key in self.calcObjDict.keys():
            outPaths = self.calcObjDict[key].outPaths
            currBModDict = fitBMod.getBulkModFromOutFilesAseWrapper(outPaths,eos=eos)
            bmodDict[key] = currBModDict
        self.fittedEos = bmodDict
    
    def appendBlankEos(self):
        blankDict = dict()
        for key in self.calcObjDict.keys():
            blankDict[key] = {k:np.nan for k in ["v0", "b0", "e0"]}
            outPaths = self.calcObjDict[key].outPaths
            energies,vols = fitBMod.getVolAndEnergiesForASEFromOutFileList(outPaths)
            blankDict[key]["data"] = np.array( (energies,vols) )
        self.fittedEos = blankDict
    
    @property
    def bulkModDict(self):
        try:
            outDict = {k:v["b0"] for k,v in self.fittedEos.items()}
            return outDict
        except (KeyError, AttributeError):
            return None
    
    @property
    def volDict(self):
        try:
            outDict = {k:v["v0"] for k,v in self.fittedEos.items()}
            return outDict
        except (KeyError, AttributeError):
            return None
    
    @property
    def deltaE0Dict(self):
        absE0 = self._getAbsE0Dict()
        minE0 = min(absE0.values())
        relE0 = {k:(v["e0"]-minE0) for k,v in self.fittedEos.items()}
        return relE0
    
    @property
    def deltaE0EvolPlotDict(self):
        absE0Dict = self._getAbsE0Dict()
        outDict = dict()
        for key in self.fittedEos.keys():
            outDict[key] = copy.deepcopy(self.fittedEos[key]["data"])
            outDict[key][:,1] = outDict[key][:,1] - min(absE0Dict.values())
        return outDict

    
    def _getAbsE0Dict(self):
        return {k:v["e0"] for k,v in self.fittedEos.items()}
   






#Creating this will need all structures in essence. We need to just hold a list of calcObjs,
#which i must generate in some other way
class EvolSingleCrystType():
    def __init__(self, calcObjList):
        self.calcObjList = calcObjList
    
    def writeFiles(self):
        for x in self.calcObjList:
            x.writeFile()
    
    def getRunComms(self):
        outList = list()
        for x in self.calcObjList:
            outList.append( x.getRunComm() )

        return outList
            
    @property
    def outPaths(self):
        return [x.filePath + x.outFileExt for x in self.calcObjList]
    
    @property
    def parsedFiles(self):
        return [x.parseOutFile() for x in self.calcObjList]
        
def _getCalcObjsListFromStructListAndFolderAndMethodObj(structList, saveFolder, methodObj):
    outList = list()
    fNameFormat = "aval_{}.in"
    factFunct = calcMethods.getPlatoCalcObjFromInpPathAndStrDictAndRunCommFunction
    for currStruct in structList:
        aVal = currStruct.lattParams["a"]
        fName = fNameFormat.format( "{:.3f}".format(aVal).replace(".","pt") )
        outPath = os.path.abspath (  os.path.join( saveFolder,fName )  )
        currStrDict = methodObj.getStrDictWithStruct(currStruct)
        currObj = factFunct(outPath, currStrDict, methodObj.runCommFunction)
        outList.append(currObj)
    return outList

