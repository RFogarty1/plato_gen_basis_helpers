

''' Very specific code with the purpose of helping to run/analyse calculations to get bulk modulii/eqm volumes '''
''' etc. for Mg perfect crystals '''

from collections import OrderedDict
import copy
import itertools as it
import os
import sys

import numpy as np

import calc_methods as calcMethods
import ref_data_mg as refData

sys.path.append('/media/ssd1/rf614/Work/usr_scripts/coding/Plato_Analysis_Lib_Functions')
import fit_bulk_mod as fitBMod




def getStandardMultiCrystsDictFromMethodStrsAndDataSetFolder(methodStrs:list, dataSetFolder:"str to tell plato where to find model", refStr=None):
	startFolder = getStandardStartFolderBulkModCalcsBasedOnCurrDir()
	kPts = getStandardKPtDict()
	structDict = getStandardStructDictHcpFccBcc()
	gridDict = getStandardGridDictsForMethods(methodStrs)
	startDict = {"kpts": kPts,
	             "startFolder": startFolder,
	             "structDict": structDict,
	             "dataFolder": dataSetFolder}
	optsDict = createAllEvolOptsDictFromStartDictAndMethodStrs(startDict, methodStrs, gridDict)
	allMultiCrystsDict = copy.deepcopy(optsDict) #Overly cautious really
	for key,val in allMultiCrystsDict.items():
		allMultiCrystsDict[key] = val.getMultiCrystalObj() 
#	allMultiCrystsDict = {k:v.getMultiCrystalObj() for k,v in optsDict.items()}

	#Create a mock object for the reference string if it isnt present
	if refStr is not None:
		if refStr in methodStrs:
			pass
		else:
			allMultiCrystsDict[refStr] = getRefBModMockedMultiCrystal()

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
	mockObj = EvolMultiCrystType({k:None for k in bmodDict},label="reference")
	mockObj.writeFiles = lambda: None
	mockObj.getRunComms = lambda: list()
	mockObj.fitEos = _decorateFunctToDoNothing(mockObj.fitEos)
	mockObj.fittedEos = bmodDict
	return mockObj

def _getRefBModDict():
    refEosDict = {"hcp": refData.getPlaneWaveEosFitDict("hcp"),
                    "bcc": refData.getPlaneWaveEosFitDict("bcc"),
                    "fcc": refData.getPlaneWaveEosFitDict("fcc")}
    return refEosDict

def _decorateFunctToDoNothing(funct):
	def outFunct(*args,**kwargs):
		return None
	return outFunct


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
	outDict = OrderedDict()
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
        return EvolMultiCrystType(multiCrystDict, label=self.methodStr)
 


class EvolMultiCrystsAndMethods():

	def __init__(self, multiCrystObjList, label=""):
		self.objList = multiCrystObjList
		self.label = label

	def writeFiles(self):
		for obj in self.objList:
			obj.writeFiles()

	def getRunComms(self):
		outList = list()
		for obj in self.objList:
			outList.extend( obj.getRunComms() )
		return outList

	def fitEos(self,eos,fitBlankIfError=None):
		for obj in self.objList:
			if fitBlankIfError is None:
				obj.fitEos(eos)
			else:
				obj.fitEos(eos,fitBlankIfError=fitBlankIfError)

	@property
	def _tableHeadings(self):
		return ["Method", "v0/ bohr cubed", "b0 / GPa", "Delta e0 / eV"]


	def getTablesForTabulate(self):
		tabDictList = list()		
		for obj in self.objList:
			tabDictList.append(obj.getTableDictForTabulate())

	
		uniqueKeys = self._getUniqueKeysFromDictList(tabDictList)	

		outTabList = list()
		print("uniqueKeys = {}".format(uniqueKeys))
		for key in uniqueKeys:
			currKeyTab = [[key]]
			currKeyTab.append(self._tableHeadings)
			for currDict in tabDictList:
				try:
					currTab = currDict[key]
				except KeyError:
					pass
				else: #This runs only when the exception doesnt occur
					currTab.pop(0)
					currTab.pop(0)
					currKeyTab.extend(currTab)
			outTabList.append(currKeyTab)

		return outTabList

	def _getUniqueKeysFromDictList(self, tabDictList):
		allKeys = list()
		for currDict in tabDictList:
			allKeys.extend(currDict.keys())	
		return list(set(allKeys))	


	def getPlotDataDict(self, keyOrder=None):
		outDict = dict()
		for currObj in self.objList:
			currDict = currObj.deltaE0EvolPlotDict
			if keyOrder is None:
				currDict = [val for k,v in currDict.items()]
			else:
				currDict = [currDict[k] for k in currDict.keys()]
			outDict[currObj.label] = currDict #RHS is actually now a list			
		return outDict





class EvolMultiCrystType():
	def __init__(self, calcObjDict, label=""):
		self.calcObjDict = calcObjDict
		self.label = label
	
	def writeFiles(self):
		for key in self.calcObjDict:
			self.calcObjDict[key].writeFiles()
	
	def getRunComms(self):
		allRunComms = list()
		for key in self.calcObjDict:
			currRunComms = self.calcObjDict[key].getRunComms()
			allRunComms.extend(currRunComms)
		return allRunComms
	
	#Eos is essentially the 2nd function of this class, and makes it tricky	
	def fitEos(self, eos, fitBlankIfError=False):
		bmodDict = dict()
		for key in self.calcObjDict.keys():
			outPaths = self.calcObjDict[key].outPaths
			try:
				currBModDict = fitBMod.getBulkModFromOutFilesAseWrapper(outPaths,eos=eos)
				bmodDict[key] = currBModDict
			except RuntimeError as e:
				if fitBlankIfError:
					bmodDict[key] = self.appendBlankEos()
				else:
					raise(e)
		self.fittedEos = bmodDict
	
	def appendBlankEos(self):
		blankDict = dict()
		for key in self.calcObjDict.keys():
			blankDict[key] = {k:np.nan for k in ["v0", "b0", "e0"]}
			outPaths = self.calcObjDict[key].outPaths
			energies,vols = fitBMod.getVolAndEnergiesForASEFromOutFileList(outPaths)
			blankDict[key]["data"] = np.array( (energies,vols) )
		self.fittedEos = blankDict


	def getTableDictForTabulate(self):
		tabDict = dict()
		for key in self.calcObjDict.keys():
			tabList = list()
			tabList.append(self.label)
			colHeadings = self._tableHeadings
			v0,b0,e0 = self.volDict[key], self.bulkModDict[key], self.deltaE0Dict[key]
			strList = [ "{:.3f}".format(x) for x in [v0,b0,e0] ]
			values = [self.label] + strList
			tabList.append(colHeadings)
			tabList.append(values)
			tabDict[key] = tabList
		return tabDict


	@property
	def _tableHeadings(self):
		return ["Method", "v0/ bohr cubed", "b0 / GPa", "Delta e0 / eV"]
	
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

