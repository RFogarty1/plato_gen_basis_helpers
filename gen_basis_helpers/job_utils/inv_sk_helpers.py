
import itertools
import os
import shutil
import subprocess


import numpy as np
from concurrent import futures

import plato_pylib.plato.mod_plato_inp_files as platoInp
from ..shared import ch_dir as chDir

import sys

import plato_pylib.plato.parse_inv_sk as parseInvSK
import plato_pylib.plato.parse_tbint_files as parseTbint

'''Purpose of this code is to make it easier to run/analyse tb2 calculations (mainly inv-sk) '''


#----------------->Functions for running inv-sk calculations <-----------------------

def runInvSkParralel(inpFilePaths, nCores):

	actStartDir = os.getcwd()
	startDirs = list()
	absFilePaths = [os.path.abspath(x) for x in inpFilePaths]
	for inpPath in absFilePaths:
		currStart = os.path.split(inpPath)[0]
		startDirs.append(currStart)

	jobsLeft = len(absFilePaths)
	print("A total of {} jobs will be run".format(jobsLeft))
	numbThreads = min( jobsLeft, nCores )

	with futures.ThreadPoolExecutor(numbThreads) as executor:
		for i in executor.map(runInvSk, absFilePaths, startDirs):
			jobsLeft -= 1
			print("jobsLeft = {}".format(jobsLeft))

	os.chdir(actStartDir)

def runInvSk(inpPath, startDir=None):
	if startDir is None:
		startDir = os.getcwd()

	#Step 1 = get base file name
	folder, fullFName = os.path.split(inpPath)
	baseFName, unused = os.path.splitext(fullFName)
	if folder=="":
		folder = startDir

	#Step 2 = create a temporary folder + copy file over
	tempDir = os.path.join(folder, baseFName)
	runFilePath = os.path.join(tempDir,fullFName)
	os.mkdir(tempDir) #It SHOULDNT ever already exist; hence not supressing the file exists error
	shutil.copy2(inpPath, runFilePath)

	#Step 3 = modify inversesk flag in file [to make sure its set] + figure out which elements are present
	_modInvSkFlagToOn(runFilePath)

	#Step 4 = Figure out what the output inv-sk filepaths will be[use a single function ldo]
	outputInvSkFileNames = _getOutputInvSkFileNames(runFilePath)

	#Step 5 = run the actual inv-sk calculation [Mod the file to use inverse-SK if not already present]
	with chDir.ChDir(startDir, tempDir):
		subprocess.check_call("cd {};tb2 {}".format(tempDir,baseFName),shell=True) #Often breaks when the cd {} is missing

	#Step 6 = rename the relevant inv-sk output files to something unique
	for x in outputInvSkFileNames:
		filePath = os.path.join(tempDir,x)
		outPath = os.path.join(tempDir, baseFName + "_" + x)
		os.rename( filePath, outPath )

	#Step 7 = move back all files except the input file (NOT copy)
	for fName in os.listdir(tempDir):
		if fName.endswith(".in"):
			os.remove( os.path.join(tempDir,fName) )
		else:
			fPath = os.path.join(tempDir,fName)
			os.rename(fPath, os.path.join(folder,fName))

	#Step 8 = delete temporary directory
	os.rmdir(tempDir)


def _modInvSkFlagToOn(inpFilePath):
	tokenizedFile = platoInp.tokenizePlatoInpFile(inpFilePath)
	outOptDict = {k.lower():v for k,v in tokenizedFile.items()}
	outOptDict["inversesk"] = "1"
	platoInp.writePlatoOutFileFromDict(inpFilePath,outOptDict)


def _getOutputInvSkFileNames(inpFilePath):
	tokenizedFile = {k.lower():v for k,v in platoInp.tokenizePlatoInpFile(inpFilePath).items()}
	if int(tokenizedFile["format"]) != 0:
		raise ValueError("Geometry format needs to be {}; but format={}".format(0,tokenizedFile["format"]))	

	#Get the elements present
	nAtoms = tokenizedFile["natom"]
	atoms = tokenizedFile["atoms"].split("\n")
	elementsPresent = list()
	for x in atoms:
		currElement = x.split()[-1]
		if currElement not in (elementsPresent):
			elementsPresent.append(currElement)

	#Get the file names from them
	outFileNames = [ "{}_{}_SK.csv".format(a,b) for a,b in itertools.product(elementsPresent,elementsPresent)   ]

	return outFileNames



#----------------->Functions for parsing inv-sk calculations <-----------------------

def getTbintHoppingListStructFormat(tbintPath,intType="hopping"):
	allInts = parseTbint.getIntegralsFromBdt(tbintPath)
	adtPaths = parseTbint.getAdtFilePathsFromBdt(tbintPath)
	shellsA, shellsB = parseTbint.parseAdtFile(adtPaths[0])["numbShells"], parseTbint.parseAdtFile(adtPaths[1])["numbShells"]

	#Put all the info in the required structure (TODO: extract this method somehow w/strat patt.)
	allBondTypes = ["sigma","pi","delta"]
	outLists = list()
	for shellAIdx in range(shellsA):
		currListA = list()
		for shellBIdx in range(shellsB):
			currDict = dict()
			for bondType in allBondTypes:
				 currDict[bondType] = getTbintIntegralTableFromObjDict(allInts, shellAIdx, shellBIdx, bondType, intType=intType)
			currListA.append(currDict)
		outLists.append(currListA)

	return outLists


def getXvsYStructsFromSKFileList(inpFileList, xProp="distance", yProp="hVal".lower(), tbintPath=None, inclXtal=True):
	individualListStructs = [ getXvsYStructsFromSKFile(currPath, xProp=xProp, yProp=yProp, tbintPath=tbintPath,inclXtal=inclXtal) for currPath in inpFileList ]
	mergedLists = individualListStructs[0]

	for x in individualListStructs[1:]:
		_appendListStruct(mergedLists,x)

	return mergedLists


def _appendListStruct(baseListStruct, appendStruct):
	for shellA in range(len(baseListStruct)):
		for shellB in range(len(baseListStruct[0])):
			for bondType in baseListStruct[shellA][shellB].keys():
				baseArray = baseListStruct[shellA][shellB][bondType]
				if baseArray is not None:
					baseListStruct[shellA][shellB][bondType] = np.concatenate( (baseArray,appendStruct[shellA][shellB][bondType]) )



def getXvsYStructsFromSKFile(inpPath, xProp="distance", yProp="hVal".lower(), tbintPath=None, inclXtal=True):
	''' Access np array you want by outLists[shellAIdx][shellBIdx][bondType], e.g. outLists[0][1]["pi"] '''
	''' Bond Type can be "sigma", "pi", "delta" '''

	parsedFileObj = parseInvSK.parseInvSK(inpPath)
	if inclXtal is False:
		parsedFileObj.removeXtalFieldTerms()

	return getXvsYStructsFromSkObj(parsedFileObj, xProp="distance", yProp="hVal".lower())



def getXvsYStructsFromSkObj(parsedFileObj, xProp="distance", yProp="hVal".lower(), tbintPath=None):
	allBondTypes = ["sigma","pi","delta"]

	#Need to figure out all shell combinations
	shellIndA = parsedFileObj.getShellIndices(atomIdx=0)
	shellIndB = parsedFileObj.getShellIndices(atomIdx=1)

	#For looking at errors vs tbint we need the tbint output file
	if tbintPath is not None:
		tbintObjs = parseTbint.getIntegralsFromBdt(tbintPath)
	else:
		tbintObjs = None

	#Put all the info in the required structure
	outLists = list()
	for shellAIdx in range(len(shellIndA)):
		currListA = list()
		for shellBIdx in range(len(shellIndB)):
			currDict = dict()
			for bondType in allBondTypes:
				 currDict[bondType] = getInvSKDataTwoShells(parsedFileObj, shellAIdx, shellBIdx, bondType, xProp, yProp, tbintDataObjs=tbintObjs)
			currListA.append(currDict)
		outLists.append(currListA)

	return outLists


#---------->Lower level inv-sk parsers<--------------


def getInvSKDataTwoShells(skDataObj, shellA, shellB, bondType, xProp, yProp, tbintDataObjs:dict=None):

	#Inefficient way of checking whether shell/bondtype combination is valid
	testData = skDataObj.getAllValsOrbPair("hval",shellA,shellB,bondType=bondType) 
	if (testData[0][1] is None):
		return None

	xData = _getPropFromInvSKData(skDataObj, xProp.lower(), shellA, shellB, bondType,tbintDataObjs)
	yData = _getPropFromInvSKData(skDataObj, yProp.lower(), shellA, shellB, bondType,tbintDataObjs)

	zippedXY = [(x,y) for x,y in itertools.zip_longest(xData,yData)]

	return np.array(zippedXY)


def getTbintIntegralTableFromObjDict( tbintDataObjs:dict, shellA, shellB, bondType, intType="hopping"):
	bondTypeToOrbSubIdx = {"sigma":1,"pi":2,"delta":3}
	orbSubIdx = bondTypeToOrbSubIdx[bondType]

	for integral in tbintDataObjs[intType]:
		if (integral.shellA==shellA) and (integral.shellB==shellB) and (integral.orbSubIdx==orbSubIdx):
			return integral.integrals
	return None


#NOTE: SOME VERY INEFFICIENT PARTS TO THIS
def _getPropFromInvSKData(skDataObj, prop, shellA, shellB, bondType=None,tbintDataObjs=None):
	if prop == "distance":
		outData = [x for x,y in skDataObj.getAllValsOrbPair("hval", shellA, shellB, bondType=bondType)]
	elif prop == "screenfunct":
		outData = [y for x,y in skDataObj.getAllValsOrbPair("screenFunct", shellA, shellB, bondType=bondType)]
	elif prop == "tanh_screenfunct":
		argVals = [y for x,y in skDataObj.getAllValsOrbPair("screenFunct", shellA, shellB, bondType=bondType)]
		outData = [math.tanh(x) for x in argVals]
	elif prop == "screenfunct_times_hr" or prop=="tanh_screenfunct_times_hr":
		currTbintTable = getTbintIntegralTableFromObjDict( tbintDataObjs, shellA, shellB, bondType )
		screenFData = [y for x,y in skDataObj.getAllValsOrbPair("screenFunct", shellA, shellB, bondType=bondType)]
		if prop=="tanh_screenfunct_times_hr":
			screenFData = [math.tanh(x) for x in screenFData]
		distanceData = [x for x,y in skDataObj.getAllValsOrbPair("screenFunct", shellA, shellB, bondType=bondType)]
		outData = [_getNearestTbintValue(dist,currTbintTable)*sFunct for dist,sFunct in itertools.zip_longest(distanceData,screenFData)]
	elif prop == "screenFunctAngDep_times_hr".lower():
		currTbintTable = getTbintIntegralTableFromObjDict( tbintDataObjs, shellA, shellB, bondType )
		screenFData = [y for x,y in skDataObj.getAllValsOrbPair("screenFunctAngDep".lower(), shellA, shellB, bondType=bondType)]
		distanceData = [x for x,y in skDataObj.getAllValsOrbPair("screenFunctAngDep".lower(), shellA, shellB, bondType=bondType)]
		outData = [_getNearestTbintValue(dist,currTbintTable)*sFunct for dist,sFunct in itertools.zip_longest(distanceData,screenFData)]
	elif prop == "hval":
		outData = [y for x,y in skDataObj.getAllValsOrbPair("hval", shellA, shellB, bondType=bondType)]
	elif prop == "hval_tbint":
		currTbintTable = getTbintIntegralTableFromObjDict( tbintDataObjs, shellA, shellB, bondType )
		distVals = [x for x,y in skDataObj.getAllValsOrbPair("hval", shellA, shellB, bondType=bondType)]
		outData = [_getNearestTbintValue(dist,currTbintTable) for dist in distVals]
	elif prop == "herror":
		distVsHVals = skDataObj.getAllValsOrbPair("hval", shellA, shellB, bondType=bondType)
		currTbintTable = getTbintIntegralTableFromObjDict( tbintDataObjs, shellA, shellB, bondType )
		outData = list()
		for x,y in distVsHVals:
			tbintVal = _getNearestTbintValue(x , currTbintTable)
			outData.append(tbintVal-y)
	else:
		raise ValueError("{} is an invalid property keyword".format(prop))

	return outData





def _getNearestTbintValue(dist, tbintDataTable:"iter; row 1=distances, row 2 = int vals"):
	#Convert to np first to get the differences faster
	npTable = np.array(tbintDataTable)
	idx = (np.abs(npTable[:,0] - dist)).argmin()
	return tbintDataTable[idx][1]

	#Can be a good 100x slower for real data-sets
#	diffs = [abs(x-dist) for x,unused in tbintDataTable] 
#	minDiff = min(diffs)
#	minDiffIdx = diffs.index(minDiff)
#	return tbintDataTable[minDiffIdx][1]




