
from ..analyse_md import thermo_data as thermoDataObjs

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

	return outDict

def _parseThermoSection(fileAsList, lineIdx):

	#1) Get headers and which columns they are
	outDict = dict()
	supportedHeaders = ["Step", "Temp", "TotEng", "Press"]
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
	headersToKeys = {"Step":"step", "Temp":"temp", "TotEng":"eTotal","Press":"pressure"}
	outKwargDict = dict()
	for key in outDict.keys():
		outKwargDict[ headersToKeys[key] ] = outDict[key]

	outObj = thermoDataObjs.ThermoDataStandard.fromStdKwargs(**outKwargDict)

	return lineIdx, outObj


def _getFileAsListFromInpPath(inpPath):
	with open(inpPath,"rt") as f:
		outStr = f.read()
	return outStr.split("\n")


