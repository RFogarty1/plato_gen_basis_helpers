

import copy
import itertools as it
import os

import plato_pylib.parseOther.parse_cp2k_files as parseCP2KHelp
import plato_pylib.shared.unit_convs as uConvHelp

from ..misc import nudged_band_paths as nebPathHelp

def getNebPathFromParsedFileObj(parsedFileObj):
	""" Gets a NudgedBandPath object from a parsedFile
	
	Args:
		parsedFileObj: (types.SimpleNamespace) This has .final_neb_summary .neb_geoms and .neb_energies attributes
			 
	Returns
		nebPath: (NudgedBandPathStandard) Object containing the neb path. Has a toDict/fromDict method for easy serialization
 
	"""
	#Get the distance of each image along the pathway, rather than distance to previous image
	distsFromPrev = parsedFileObj.final_neb_summary["dists"]
	outDists = [0]
	for dist in distsFromPrev:
		outDists.append( dist+outDists[-1] )
	geoms = parsedFileObj.neb_geoms
	energies = parsedFileObj.neb_energies

	#Get the steps
	outSteps = list()
	for geo, energyObj, dist in it.zip_longest(geoms, energies, outDists):
		print(geo, energyObj, dist)
		currStep = nebPathHelp.NudgedBandStepStandard(geom=geo, energies=energyObj, dist=dist)
		outSteps.append(currStep)

	return nebPathHelp.NudgedBandPathStandard( outSteps )



def parseNudgedBandCalcStandard(cpoutPath, convAngToBohr=False):
	""" Parse info from the various Nudged elastic band output files, with the *.cpout file path as the sole input
	
	Args:
		cpoutPath: (str) The path to the output *.cpout file
			 
	Returns
		outDict: (dict) Contains all the parsed info
 
	"""
	outDict = parseNebSummaryFile(cpoutPath)
	numbReplicas = len(outDict["final_neb_summary"]["energies"])

	#Get the geometries
	nebGeoms = _getNebGeomsFromCpoutPathAndNumbReplicasAndEmptyCell(cpoutPath, numbReplicas, outDict["unitCell"])

	#Get the energies
	nebEnergies = _getEnergiesFromCpoutPathAndNumbReplicas(cpoutPath, numbReplicas)

	#Merge everything
	outDict["neb_geoms"] = nebGeoms
	outDict["neb_energies"] = nebEnergies

	#Convert ang to bohr if needed
	if convAngToBohr:
		[x.convAngToBohr() for x in outDict["neb_geoms"]]
		outDict["unitCell"].convAngToBohr()

	return outDict


def _getNebGeomsFromCpoutPathAndNumbReplicasAndEmptyCell(cpoutPath, nReplicas, inpCell):
	xyzPaths = _getXyzFilesFromCpoutPathAndNumbReplicas(cpoutPath, nReplicas)
	cartCoords = _parseFinalGeomFromIterOfXyzFiles(xyzPaths)
	outGeoms = list()
	for currCoords in cartCoords:
		currGeom = copy.deepcopy(inpCell)
		currGeom.cartCoords = currCoords
		outGeoms.append(currGeom)
	return outGeoms

def _getEnergiesFromCpoutPathAndNumbReplicas(cpoutPath, nReplicas):
	bandOutPaths = _getBandOutFilesFromCpoutPathAndNumbReplicas(cpoutPath, nReplicas)
	energyObjs = _parseFinalEnergyObjsFromIterOfOutFiles(bandOutPaths)
	return energyObjs

def _parseFinalGeomFromIterOfXyzFiles(xyzFiles):
	outGeoms = list()
	for xyzFile in xyzFiles:
		currGeoms = parseCP2KHelp.parseXyzFromGeomOpt(xyzFile, startGeomIdx=0)
		finalGeom = currGeoms["all_geoms"][-1].cartCoords
		outGeoms.append(finalGeom)
	return outGeoms

def _parseFinalEnergyObjsFromIterOfOutFiles(outFiles):
	outEnergies = list()
	for outFile in outFiles:
		currParsed = parseCP2KHelp.parseCpout(outFile,ThrowIfTerminateFlagMissing=False) #Only last file has terminate flag
		currEnergyObj = currParsed["energies"]
		outEnergies.append(currEnergyObj)
	return outEnergies

def _getXyzFilesFromCpoutPathAndNumbReplicas(cpoutPath, numbReplicas):
	folder, filename = os.path.split(cpoutPath)
	baseFileName = os.path.splitext(filename)[0]
	numbDigits = len(str(int(numbReplicas)))
	extFmt = "-pos-Replica_nr_{:0" + str(numbDigits) + "}-1.xyz"
	outPaths = [os.path.join(folder, baseFileName+extFmt.format(x)) for x in range(1,numbReplicas+1)] 
	return outPaths

def _getBandOutFilesFromCpoutPathAndNumbReplicas(cpoutPath, numbReplicas):
	folder, filename = os.path.split(cpoutPath)
	baseFileName = os.path.splitext(filename)[0]
	numbDigits = len(str(int(numbReplicas)))
	extFmt = "-BAND{:0" + str(numbDigits) + "}.out"
	outPaths = [os.path.join(folder, baseFileName+extFmt.format(x)) for x in range(1,numbReplicas+1)] 
	return outPaths

def parseNebSummaryFile(inpPath):
	""" Parses the summary *.cpout file for a CP2K Nudged elastic band calculation
	
	Args:
		inpPath: Path to *.cpout file
			 
	Returns
		outDict: (dict) Contains all info parsed from summary file
 
	"""
	fileAsList = _getFileAsListFromInpFile(inpPath)
	outParser = parseCP2KHelp.CpoutFileParser()
	parseCP2KHelp._addSearchWordAndFunctToParserObj("BAND TYPE", _parseNebSummarySection, outParser, handleParsedDictFunct=_handleParsedNebSummarySection)
	return outParser.getOutDictFromFileAsList(fileAsList)


def _getFileAsListFromInpFile(inpFile):
	with open(inpFile,"rt") as f:
		fileAsList = f.readlines()
	return fileAsList


def _parseNebSummarySection(fileAsList, lineIdx):

	outDict = dict()
	endStr = "*****"
	haToEv = uConvHelp.RYD_TO_EV*2

	while (lineIdx<len(fileAsList)) and (endStr not in fileAsList[lineIdx]):
		currLine = fileAsList[lineIdx]
		currSplitLine = currLine.strip().split()

		if "BAND TYPE" in currLine:
			if "OPTIMIZATION" not in currLine:
				outDict["bandType"] = currSplitLine[-1]
		elif "STEP NUMBER" in currLine:
			outDict["stepNumber"] = int(currSplitLine[-1])
		elif "BAND TOTAL ENERGY" in currLine:
			outDict["bandEnergy"] = haToEv*float(currSplitLine[-1])

		elif "DISTANCES REP" in currLine:
			startIdx = lineIdx
			dists = list()
			currLine = fileAsList[lineIdx].replace("DISTANCES REP =","")
			while "ENERGIES" not in currLine:
				currSplitLine = currLine.strip().split()
				currVals = [float(x) for x in currSplitLine]
				dists.extend(currVals)
				lineIdx += 1
				currLine = fileAsList[lineIdx].replace("DISTANCES REP =","")

			outDict["dists"] = dists
			lineIdx = startIdx

		elif "ENERGIES [au]" in currLine:
			startIdx = lineIdx
			energies = list()
			currLine = fileAsList[lineIdx].replace("ENERGIES [au] =","")
			while "BAND" not in currLine:
				currSplitLine = currLine.strip().split()
				currVals = [float(x)*haToEv for x in currSplitLine]
				energies.extend(currVals)
				lineIdx += 1
				currLine = fileAsList[lineIdx].replace("ENERGIES [au] =","")

			outDict["energies"] = energies
			lineIdx = startIdx

		lineIdx += 1

	return outDict, lineIdx

def _handleParsedNebSummarySection(parserInstance, outDict):
	parserInstance.outDict["final_neb_summary"] = outDict




