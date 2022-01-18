import math
import copy
import itertools as it
import time

from ..analyse_md import thermo_data as thermoDataObjs
from ..analyse_md import traj_core as trajObjHelp
from concurrent import futures

import multiprocessing

import plato_pylib.shared.ucell_class as uCellHelp

def parseLammpsLogFile(inpPath):
	fileAsList = _getFileAsListFromInpPath(inpPath)

	outDict = dict()
	lIdx = 0

	while lIdx < len(fileAsList):
		line = fileAsList[lIdx]
		if "timestep " in line:
			timeStep = float( line.strip().split()[-1] )
		elif "Step " in line:
			lIdx,thermoObj = _parseThermoSection(fileAsList,lIdx)

		lIdx+=1

	#Get final dict
	time = [x*timeStep for x in thermoObj.dataDict["step"]]
	thermoObj.dataDict["time"] = time
	outDict["thermo_data"] = thermoObj
	outDict["timestep"] = timeStep

	return outDict

def _parseThermoSection(fileAsList, lineIdx):

	#1) Get headers and which columns they are
	outDict = dict()
	supportedHeaders = ["Step", "Temp", "TotEng", "Press", "PotEng", "KinEng"]
	headers = fileAsList[lineIdx].strip().split()
	headerToCol = dict()
	for header in headers:
		if header in supportedHeaders:
			headerToCol[header] = headers.index(header)
			outDict[header] = list()
	lineIdx+=1


	#2)Parse the data from the file
	while lineIdx<len(fileAsList):
		line = fileAsList[lineIdx]
		if "Loop" in line:
			break

		splitLine = line.strip().split()
		for key in outDict:
			currCol = headerToCol[key]
			currVal = float( splitLine[currCol] )
			outDict[key].append(currVal)
		lineIdx+=1


	#3)Translate this into the thermo object and return that 
	headersToKeys = {"Step":"step", "Temp":"temp", "TotEng":"eTotal","Press":"pressure",
	                 "PotEng":"ePot", "KinEng":"eKinetic"}
	outKwargDict = dict()
	for key in outDict.keys():
		outKwargDict[ headersToKeys[key] ] = outDict[key]

	outObj = thermoDataObjs.ThermoDataStandard.fromStdKwargs(**outKwargDict)

	return lineIdx, outObj


def getTrajectoryFromLammpsDumpFile(inpFile, timeStep=None, typeIdxToEle=None, nCores=1):
	""" Parses the trajectory from a lammps file into a TrajectoryInMemory object. MUST BE ATOMS format
	
	Args:
		inpFile: (str) Path to the file contaiing the trajectory (*.lammpstrj file)
		timeStep: (float, Optional) Time for one step. Needed to figure out the simulation time at any trajectory step, default is to just set that value to None
		typeIdxToEle: (dict) Keys are str versions of integers (e.g. str(4)) while values are the elements those integers correspond to in the dump file
		nCores: (int) Number of cores to parralelise over

	Returns
		outTraj: (TrajectoryInMemory obj) Contains trajectory info
 
	"""
	fileAsList = _getFileAsListFromInpPath(inpFile)
	typeIdxToEleToUse = {int(k):v for k,v in typeIdxToEle.items()} if typeIdxToEle is not None else None

	#Remove any blank lines at the end of fileAsList:
	nBlankLines = 0
	for idx,line in enumerate(fileAsList[-1:0]):
		strippedLine = line.strip()
		if line=="":
			break
		else:
			nBlankLines += 1
	fileAsList = fileAsList[:-1-nBlankLines]

	#TEMPORARY: Figure out all the lines containing cartesian co-ordinates and unit-cell params
	lineIdx = 0 
	startSectionIndices = list()
	while lineIdx <len(fileAsList):
		line = fileAsList[lineIdx]
		if "ITEM: TIMESTEP" in line:
			startSectionIndices.append(lineIdx) 
		lineIdx += 1

	#TEMPORARY: Parse the individual start indices
	nSteps = len(startSectionIndices)
	trajStepIndices = [x for x in range(nSteps)]
	endSectionIndices = [x for x in startSectionIndices[1:]] + [len(fileAsList)]
	fileSlices = [ copy.deepcopy(fileAsList[sIdx:eIdx]) for sIdx,eIdx in zip(startSectionIndices,endSectionIndices) ]
	outTrajSteps = [None]*nSteps

	#Parrelel attempt B: Chunk the work into 1/nCore chunks based on the number of fileSlices:
	chunkSize = int(nSteps/nCores) if nSteps>nCores else 1

	if nCores==1:
		outList = list()
		inpArg = zip(trajStepIndices,fileSlices)
		output = _parseNSections(inpArg)

		for (idx,val) in output:
			outTrajSteps[idx] = val
		outObj = trajObjHelp.TrajectoryInMemory(outTrajSteps)

#	if nCores==1:
#		outObjs = list()
#		lineIdx = 1 #Start +1 from the "ITEM: TIMESTEP" part
#		while lineIdx < len(fileAsList):
#			lineIdx, currObj = _parseNextSectionDumpFile(fileAsList, lineIdx) #Ends on ITEM:TIMESTEP
#			outObjs.append(currObj)
#			lineIdx += 1
#		outObj = trajObjHelp.TrajectoryInMemory(outObjs)


	if nCores>1:
		print("RUNNING CHUNKSIZE>1 CODE")
		allInpArgs = zip(trajStepIndices,fileSlices)
		batchedArgs = _splitIterableIntoBatchesOfSizeN(allInpArgs, chunkSize)
		outList, tempList = list(), list()
		with multiprocessing.Pool(nCores) as pool:
			output = pool.map(_parseNSections,batchedArgs)

		for x in output:
			outList.extend(x)

		for (idx,val) in outList:
			outTrajSteps[idx] = val
		outObj = trajObjHelp.TrajectoryInMemory(outTrajSteps)

	#Put the time on each trajectory step if we know it
	if timeStep is not None:	
		for step in outObj:
			step.time = step.step*timeStep

	if typeIdxToEleToUse is not None:
		for step in outObj:
			currCartCoords = step.unitCell.cartCoords
			for idx,coord in enumerate(currCartCoords):
				currCartCoords[idx][-1] = typeIdxToEleToUse[ int(currCartCoords[idx][-1]) ]
			step.unitCell.cartCoords = currCartCoords

	return outObj



#Function for splitting iterable into batches
def _splitIterableIntoBatchesOfSizeN(iterable, batchSize, pad=False, fillvalue=None):
	"grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
	args = [iter(iterable)] * batchSize
	if pad:
		return it.zip_longest(fillvalue=fillvalue, *args)
	else:
		return zip(*args)


#Function allows each processor to parse multiple-sections, which should be cheaper than constantly switching jobs
def _parseNSections(allInpArgs):
	#	stepIdx, fileSlice = allInpArgs
#	time.sleep(2)


	outList = list()
	for currArgs in allInpArgs:
		currRes = _parseSingleSection(currArgs)
		outList.append(currRes)

	#Temp: Testing why we dont get a speedup (suspect its largely a memory issue)
#	for i in range(100*len(outList)):
#		math.factorial(2000)


	return outList

#Defining outside function scope makes it picklable
def _parseSingleSection(inpArgs):
	sIdx, fileSlice = inpArgs
	lineIdx = 0
#		outTrajSteps[sIdx] = _getTrajStepFromFileAsListAndLineIdx(fileAsList, lIdx)
#		outTrajSteps[sIdx] = _getTrajStepFromFileAsListAndLineIdx(fileSlice, 0)
	outStep =  _getTrajStepFromFileAsListAndLineIdx(fileSlice, 0)
	return (sIdx, outStep)





def _parseNextSectionDumpFile(fileAsList, lineIdx):
	stepNumber = int(fileAsList[lineIdx].strip())
	while lineIdx<len(fileAsList):
		line = fileAsList[lineIdx]
		if "ITEM: TIMESTEP" in line:
			break
		elif "ITEM: BOX" in line:
			outCell = _getEmptyCellFromBoxSectionOfDumpFile(fileAsList, lineIdx)
			lineIdx += 1
		elif "ITEM: ATOMS" in line:
			lineIdx, cartCoords = _getCartCoordsFromDumpSection(fileAsList, lineIdx+1)
		else:
			lineIdx += 1

	outCell.cartCoords = cartCoords
	trajStepObj = trajObjHelp.TrajStepBase(unitCell=outCell, step=stepNumber)


	return lineIdx,trajStepObj

def _getTrajStepFromFileAsListAndLineIdx(fileAsList, lineIdx):
	lineIdx += 1 #input is the line containing "ITEM: TIMESTEP". Whereas we want the line one after that
	stepNumber = int(fileAsList[lineIdx].strip())
	while lineIdx<len(fileAsList):
		line = fileAsList[lineIdx]
		if "ITEM: TIMESTEP" in line:
			break
		elif "ITEM: BOX" in line:
			outCell = _getEmptyCellFromBoxSectionOfDumpFile(fileAsList, lineIdx)
			lineIdx += 1
		elif "ITEM: ATOMS" in line:
			lineIdx, cartCoords = _getCartCoordsFromDumpSection(fileAsList, lineIdx+1)
		else:
			lineIdx += 1

	outCell.cartCoords = cartCoords
	trajStepObj = trajObjHelp.TrajStepBase(unitCell=outCell, step=stepNumber)
	return trajStepObj


def _getEmptyCellFromBoxSectionOfDumpFile(fileAsList, lineIdx):
	if len(fileAsList[lineIdx].split()) != 9:
		raise ValueError("Seems like current file is not triclinic; cant parse orthogonal cells at the moment")

	#Get the bound values
	xLoBound, xHiBound, xy = [float(x) for x in fileAsList[lineIdx+1].strip().split()]
	yLoBound, yHiBound, xz = [float(x) for x in fileAsList[lineIdx+2].strip().split()]
	zLo, zHi, yz = [float(x) for x in fileAsList[lineIdx+3].strip().split()]

	#Get the values for the cell. Formulae from https://lammps.sandia.gov/doc/Howto_triclinic.html
	xLo = xLoBound - min(0, xy, xz, xy+xz)
	xHi = xHiBound - max(0, xy, xz, xy+xz)
	yLo = yLoBound - min(0, yz)
	yHi = yHiBound - max(0, yz)

	outVects = [ [xHi-xLo, 0      , 0      ],
	             [xy     , yHi-yLo, 0      ],
	             [xz     , yz     , zHi-zLo] ]

	return uCellHelp.UnitCell.fromLattVects(outVects)

def _getCartCoordsFromDumpSection(fileAsList, lineIdx):
	#1) Get the full cart-coords assuming they are ordered by atomId
	cartCoords, atomIds = list(), list()
	while lineIdx<len(fileAsList):
		if "ITEM" in fileAsList[lineIdx]:
			break

		splitLine = fileAsList[lineIdx].strip().split()
		atomIds.append( int(splitLine[0]) )
		eleVal = str( splitLine[1] )
		currXyz = [float(x) for x in splitLine[-3:]]
		currCoords = currXyz + [eleVal]
		cartCoords.append(currCoords)
		lineIdx += 1

	#2) Get output coordinates by 
	outCartCoords = [list() for x in cartCoords]
	for inpIdx,currId in enumerate(atomIds):
		outCartCoords[currId-1] = cartCoords[inpIdx]

	lineIdx -= 1
	return lineIdx, outCartCoords

def _getFileAsListFromInpPath(inpPath):
	with open(inpPath,"rt") as f:
		outStr = f.read()
	return outStr.split("\n")


