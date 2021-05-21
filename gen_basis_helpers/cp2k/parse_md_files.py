
import copy
import itertools as it
import math

import numpy as np

import plato_pylib.parseOther.parse_cp2k_files as parseCP2KHelp
import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.shared.unit_convs as uConvHelp

from ..analyse_md import analyse_metadyn_hills as metadynHillsHelp

#TODO: Relative paths please
import gen_basis_helpers.analyse_md.thermo_data as thermoDataHelp
import gen_basis_helpers.analyse_md.traj_core as trajHelp


def parseMdInfoFromMultipleCpoutAndXyzPaths(cpoutPaths, xyzPaths, tempKindPaths=None, velocityPaths=None, forcePaths=None):
	""" Parse an MD trajectory spread over multiple cp2k restarts
	
	Args:
		cpoutPaths: (iter of str) Paths to *.cpout files
		xyzPaths: (iter of str) Paths to *.xyz files
		tempKindPaths: (iter of str, Optional) Paths to files containing atomic temperatures
		velocityPaths: (iter of str, Optional) Paths to files containing velocities as dumped by motion section (basically an xyz format)
		forcePaths: (iter of str, Optional) Paths to files containing forces as dumped by motion section (basically an xyz format)

	Returns
		outDict: Contains trajectory/thermodynamic info
 
	"""
	outDict = dict()
	tempKindPaths = [None for x in range(len(cpoutPaths))] if tempKindPaths is None else tempKindPaths
	velocityPaths = [None for x in range(len(cpoutPaths))] if velocityPaths is None else velocityPaths
	forcePaths = [None for x in range(len(cpoutPaths))] if forcePaths is None else forcePaths

	#1) Get all the dicts
	parsedDicts = list()
	for cpoutPath, xyzPath, tKindPath, vPath, fPath in it.zip_longest(cpoutPaths, xyzPaths, tempKindPaths, velocityPaths, forcePaths):
		currDict = parseFullMdInfoFromCpoutAndXyzFilePaths(cpoutPath, xyzPath, tempKindPath=tKindPath, velocityPath=vPath, forcePath=fPath)
		parsedDicts.append(currDict)

	#2) Merge all the trajectories
	allTraj = [ currDict["trajectory"] for currDict in parsedDicts ]
	mergedTraj = trajHelp.getMergedTrajInMemory(allTraj)
	outDict["trajectory"] = mergedTraj

	#3) Merge the thermo data
	allThermoData = [currDict["thermo_data"] for currDict in parsedDicts]
	mergedThermo = thermoDataHelp.getMergedStandardThermoData(allThermoData)
	outDict["thermo_data"] = mergedThermo

	return outDict

def parseFullMdInfoFromCpoutAndXyzFilePaths(cpoutPath, xyzPath, tempKindPath=None, velocityPath=None, forcePath=None):
	""" Description of function
	
	Args:
		cpoutPath: (str) Path to *.cpout file
		xyzPath: (str) Path to *.xyz file
		tempKindPath: (str, optional) Path to a file containing atomic temperature for each kind of atom
 
	"""
	outDict = parseCpoutForMDJob(cpoutPath)
	outSteps = parseCp2kMdXyzFile(xyzPath)

	#If step 0 is in the xyz then we need to add the initial trajectory/thermal info to outDict
	if outSteps[0]["step"] == 0:
#		stepZero = trajHelp.TrajStepBase(unitCell=outDict["init_md_cell"], step=0, time=0)
		stepZero = trajHelp.TrajStepFlexible(unitCell=outDict["init_md_cell"], step=0, time=0)
		outDict["trajectory"].trajSteps.insert(0, stepZero)
		for key in outDict["thermo_data"].dataDict:
			outDict["thermo_data"].dataDict[key].insert(0, outDict["init_thermo_dict"][key])


	#Combine infor in cpout and xyz files to get full geometries
	outDict["trajectory"] = _getMergedTrajectoryFromParsedCpoutAndXyz(outDict, outSteps)

	#Optionally add velocities and forces
	if velocityPath is not None:
		_parseVelocitiesAndAddToTrajSteps(velocityPath, outDict["trajectory"].trajSteps)
	if forcePath is not None:
		_parseForcesXyzAndAddToTrajSteps(forcePath, outDict["trajectory"].trajSteps)


	#If tempKindPath is present, then parse it and add to the thermo data object
	if tempKindPath is not None:
		parsedTempKind = parseAtomTempFile(tempKindPath)
		firstGeom = [x for unused,x in zip([0],iter(outDict["trajectory"]))][0].unitCell 
		idxToLabelDict = _getKindIdxToSymbolDict( [ x[-1] for x in firstGeom.cartCoords] ) 
		_addAtomicTempsToThermoDataObj(parsedTempKind, outDict["thermo_data"], idxToLabelDict=idxToLabelDict)

	assert outDict["thermo_data"].dataListLengthsAllEqual

	outDict.pop("init_thermo_dict")
	outDict.pop("init_md_cell")

	return outDict

#Also modifies in place
#Relies on step numbers being in order for both parsedCpout and parsedXyz
def _getMergedTrajectoryFromParsedCpoutAndXyz(parsedCpout, parsedXyz):
	trajCpout, trajXyz = parsedCpout["trajectory"].trajSteps, parsedXyz
	outList = list()

#	assert len(trajXyz)>=len(trajCpout)

	idxCpout, idxXyz = 0, 0

	#Will throw an index error if we reach the end of cpout without finding the expected unit cell
	while idxXyz<len(parsedXyz):
		while idxCpout<len(trajCpout):
			if (trajCpout[idxCpout].step == parsedXyz[idxXyz]["step"]):
				trajCpout[idxCpout].unitCell.cartCoords = parsedXyz[idxXyz]["coords"]
				outList.append( trajCpout[idxCpout] )
				break
			else:
				idxCpout+=1

		idxXyz+=1

	return trajHelp.TrajectoryInMemory( outList )



def parseCpoutForMDJob(outFile, parser=None, raiseIfTerminateFlagMissing=False):
	parser = _getStandardCpoutMDParser() if parser is None else parser
	if raiseIfTerminateFlagMissing is False:
		def _setTerminateFlagToTrue(instance):
			instance.outDict["terminate_flag_found"] = True
		parser.finalStepsFunctions.append(_setTerminateFlagToTrue)

	fileAsList = parseCP2KHelp._getFileAsListFromInpFile(outFile)
	outDict = parser.getOutDictFromFileAsList(fileAsList)
	return outDict

def _getStandardCpoutMDParser():
	outObj = parseCP2KHelp.CpoutFileParser()
	outObj = parseCP2KHelp.CpoutFileParser( )
	parseCP2KHelp._addSearchWordAndFunctToParserObj("GO CP2K GO", _parseMDInitSection, outObj, handleParsedDictFunct=_mdInitSectionParseDictHandler)
	parseCP2KHelp._addSearchWordAndFunctToParserObj("MD_INI| MD initialization", _parseMDInitSection_secondFmt, outObj, handleParsedDictFunct=_mdInitSectionParseDictHandler)
	parseCP2KHelp._addSearchWordAndFunctToParserObj("STEP NUMBER", _parseMdStepInfo, outObj, handleParsedDictFunct=_mdStepsParseDictHandler)
	parseCP2KHelp._addSearchWordAndFunctToParserObj("Step number", _parseMdStepInfo, outObj, handleParsedDictFunct=_mdStepsParseDictHandler)
	outObj.finalStepsFunctions.append(_mdStepsFinalStepFunct)
	return outObj

def _parseMDInitSection(fileAsList, lineIdx):
	outDict = dict()
	lineIdx += 1
	endStr = "GO CP2K GO"
	haToEv = uConvHelp.RYD_TO_EV*2

	while lineIdx<len(fileAsList):
		currLine = fileAsList[lineIdx]
		if endStr in fileAsList[lineIdx]:
			break
		if "INITIAL CELL LNTHS[bohr]" in currLine:
			splitLine = currLine.strip().split()
			outDict["lattParams"] = [float(x)*uConvHelp.BOHR_TO_ANG for x in splitLine[-3:]]
		if "INITIAL CELL ANGLS[deg]" in currLine:
			splitLine = currLine.strip().split()
			outDict["lattAngles"] = [float(x) for x in splitLine[-3:]]
		if "INITIAL POTENTIAL ENERGY[hartree]" in currLine:
			outDict["ePot"] = float( currLine.strip().split()[-1] )*haToEv
		if "INITIAL KINETIC ENERGY" in currLine:
			outDict["eKinetic"] = float( currLine.strip().split()[-1] )*haToEv
		if "TEMPERATURE" in currLine:
			outDict["temp"] = float( currLine.strip().split()[-1] )
		if "INITIAL PRESSURE[bar]" in currLine:
			outDict["pressure"] = float( currLine.strip().split()[-1] )

		lineIdx+=1

	return outDict, lineIdx+1


def _parseMDInitSection_secondFmt(fileAsList, lineIdx):
	outDict = dict()
	lineIdx += 1
	haToEv = uConvHelp.RYD_TO_EV*2

	while lineIdx <len(fileAsList):
		currLine = fileAsList[lineIdx]
		if "MD_INI" not in currLine:
			break
		if "Cell lengths [bohr]" in currLine:
			splitLine = currLine.strip().split()
			outDict["lattParams"] = [float(x)*uConvHelp.BOHR_TO_ANG for x in splitLine[-3:]]
		if "Cell angles [deg]" in currLine:
			splitLine = currLine.strip().split()
			outDict["lattAngles"] = [float(x) for x in splitLine[-3:]]
		if "Potential energy [hartree]" in currLine:
			outDict["ePot"] = float( currLine.strip().split()[-1] )*haToEv
		if "Kinetic energy [hartree]" in currLine:
			outDict["eKinetic"] = float( currLine.strip().split()[-1] )*haToEv
		if "Temperature [K]" in currLine:
			outDict["temp"] = float( currLine.strip().split()[-1] )
		lineIdx += 1

	return outDict, lineIdx


def _mdInitSectionParseDictHandler(parserInstance, outDict):
	cellDict = {"lattParams": outDict["lattParams"], "lattAngles": outDict["lattAngles"]}
	initCell = uCellHelp.UnitCell(**cellDict)

	possibleKeysThermo = ["ePot", "eKinetic", "pressure", "step", "time", "temp"]

	parserInstance.outDict["init_md_cell"] = initCell
	parserInstance.outDict["init_thermo_dict"] = {k:v for k,v in outDict.items() if k in possibleKeysThermo}
	parserInstance.outDict["init_thermo_dict"]["step"] = 0 #We will only EVER use these values if its step 0 anyway so....
	parserInstance.outDict["init_thermo_dict"]["time"] = 0
	parserInstance.outDict["traj_geoms"] = list()

def _mdStepsParseDictHandler(parserInstance, outDict):

	#1) Deal with the thermal info stuff
	possibleKeysThermo = ["ePot", "eKinetic", "pressure", "step", "time", "temp"]

	#Initialise the arrays if needed
	if parserInstance.outDict.get("thermo_arrays", None) is None:
		parserInstance.outDict["thermo_arrays"] = dict()
		for key in possibleKeysThermo:
			if key in outDict.keys():
				parserInstance.outDict["thermo_arrays"][key] = list()

	#Add the current data to the arrays
	for key in outDict.keys():
		if key in possibleKeysThermo:
			parserInstance.outDict["thermo_arrays"][key].append( outDict[key] )

	#2) If needed, deal with the lattice parameters
	currTrajGeom = copy.deepcopy( parserInstance.outDict["init_md_cell"] )
	if "lattParams" in outDict.keys():
		currTrajGeom.setLattParams(outDict["lattParams"])


	parserInstance.outDict["traj_geoms"].append(currTrajGeom)


def _mdStepsFinalStepFunct(parserInstance):
	outThermoDict = parserInstance.outDict["thermo_arrays"]
	parserInstance.outDict["thermo_data"] = thermoDataHelp.ThermoDataStandard(outThermoDict)

	trajSteps = list()
	thermoInfo = parserInstance.outDict["thermo_arrays"]
	assert len( thermoInfo["step"] ) == len( parserInstance.outDict["traj_geoms"] )

	for idx, geom in enumerate(parserInstance.outDict["traj_geoms"]):
#		currStep = trajHelp.TrajStepBase(unitCell=geom, step=thermoInfo["step"][idx], time=thermoInfo["time"][idx])
		currStep = trajHelp.TrajStepFlexible(unitCell=geom, step=thermoInfo["step"][idx], time=thermoInfo["time"][idx])
		trajSteps.append(currStep)

	parserInstance.outDict["trajectory"] = trajHelp.TrajectoryInMemory(trajSteps)

	parserInstance.outDict.pop("thermo_arrays")
	parserInstance.outDict.pop("traj_geoms")

def _parseMdStepInfo(fileAsList, lineIdx):
	outDict = dict()
	endStr = "*********"

	haToEv = uConvHelp.RYD_TO_EV*2
	while lineIdx<len(fileAsList):
		currLine = fileAsList[lineIdx].upper()
		if endStr in currLine:
			break
		if "POTENTIAL ENERGY" in currLine:
			outDict["ePot"] = float(currLine.strip().split()[-2])*haToEv
		if "KINETIC ENERGY" in currLine:
			outDict["eKinetic"] = float(currLine.strip().split()[-2])*haToEv
		if "PRESSURE" in currLine:
			outDict["pressure"] = float(currLine.strip().split()[-2])
		if "STEP NUMBER" in currLine:
			outDict["step"] = int(currLine.strip().split()[-1])
		if "TIME [FS]" in currLine:
			outDict["time"] = float(currLine.strip().split()[-1])
		if "TEMPERATURE [K]" in currLine:
			outDict["temp"] = float(currLine.strip().split()[-2])
		if "CELL LNTHS[BOHR]             " in currLine:
			outDict["lattParams"] = [float(x) for x in currLine.strip().split()[-3:]]
		if "CELL ANGLS[DEG]" in currLine: 
			outDict["lattAngles"] = [float(x) for x in currLine.strip().split()[-3:]]
		lineIdx += 1

	return outDict, lineIdx+1



def _parseVelocitiesAndAddToTrajSteps(velXyzPath, trajSteps):
	""" Gets velocities from a *.xyz dump and adds them to trajSteps when step indices match
	
	Args:
		velXyzPath: (str, Path)
		trajSteps: (iter of TrajStepFlexible)
			 
	Returns
		Nothing; works in place
 
	"""
	parsedDicts = parseCp2kMdXyzFile(velXyzPath)
	valsAsCoords = _getValuesOfXyzParsedDictsForTrajSteps(parsedDicts, trajSteps)

	for tStep,val in it.zip_longest(trajSteps,valsAsCoords):

		if val is None:
			currVals = None
		else:
			currVals = [x[:3] for x in val]

		tStep.addExtraAttrDict({"velocities": {"value":currVals, "cmpType":"numericalArray"}})

def _parseForcesXyzAndAddToTrajSteps(forcesXyzPath, trajSteps):
	""" See _parseVelocitiesAndAddToTrajSteps; this is basically the same except with atomic forces rather than velocities """

	parsedDicts = parseCp2kMdXyzFile(forcesXyzPath)
	valsAsCoords = _getValuesOfXyzParsedDictsForTrajSteps(parsedDicts, trajSteps)

	for tStep,val in it.zip_longest(trajSteps,valsAsCoords):

		if val is None:
			currVals = None
		else:
			currVals = [x[:3] for x in val]

		tStep.addExtraAttrDict({"forces": {"value":currVals, "cmpType":"numericalArray"}})


def _getValuesOfXyzParsedDictsForTrajSteps(parsedDicts, trajSteps):
	pDictIdx, tStepIdx = 0,0

	outVals = list()
	while (pDictIdx<len(parsedDicts)) and (tStepIdx<len(trajSteps)):
		pDictStep, tStepDict = parsedDicts[pDictIdx]["step"], trajSteps[tStepIdx].step

		if pDictStep==tStepDict:
			outVals.append( parsedDicts[pDictIdx]["coords"] )
			pDictIdx += 1
			tStepIdx += 1

		elif pDictStep<tStepDict:
			pDictIdx += 1

		elif pDictStep>tStepDict:
			outVals.append( None )
			tStepIdx += 1

		else:
			raise ValueError("Shouldnt ever reach here")

	#TODO: Currently untested; so i should test it really
	#Fill in any remaining values
	if len(outVals) < len(trajSteps):
		for step in trajSteps[ len(outVals): ]:
			outVals.append(None)

	return outVals




def parseCp2kMdXyzFile(inpXyz):
	fileAsList = _readFileIntoList(inpXyz)
	lineIdx = 0
	outDicts = list()
	while lineIdx < len(fileAsList):
		if fileAsList[lineIdx].split() == "":
			lineIdx += 1
		else:
			lineIdx, currDict = _parseSingleStepXyz(fileAsList, lineIdx)
			outDicts.append(currDict)

	#Want to only get results from last md-run. Adjacent step numbers not increasing is generally the sign of two jobs appending
	startIdx = 0
	stepNumbs = [x["step"] for x in outDicts]

	if len(stepNumbs)>1:
		for idx,step in enumerate(stepNumbs[1:],start=1):
			if step<= stepNumbs[idx-1]:
				startIdx = idx

	return outDicts[startIdx:]


def _readFileIntoList(inpPath):
	with open(inpPath,"rt") as f:
		outList = f.readlines()
	return outList

def _parseSingleStepXyz(fileAsList, lineIdx):
	outDict = dict()
	nAtoms = int( fileAsList[lineIdx].strip() )
	lineIdx += 1
	splitLine = fileAsList[lineIdx].strip().replace(","," ").replace("="," ").split()
	lineIdx += 1
	outDict["step"] = int(splitLine[1])
	outDict["time"] = float(splitLine[3])
	outCoords = list()

	for idx in range(nAtoms):
		splitLine = fileAsList[lineIdx].strip().split()
		currCoords = [float(x) for x in splitLine[1:4]] + [splitLine[0]]
		outCoords.append(currCoords)
		lineIdx += 1

	outDict["coords"] = outCoords

	return lineIdx, outDict


def parseAtomTempFile(inpFile):
	fileAsList = parseCP2KHelp._getFileAsListFromInpFile(inpFile)
	outDict = {"step":list(), "time":list()}
	idx = 0

	while idx<len(fileAsList):
		if idx==0:
			nKinds = len(fileAsList[idx].strip().split()) - 2
			outDict["kindTemp"] = [list() for x in range(nKinds)]

		currSplitLine = fileAsList[idx].strip().split()

		#If this step isnt > previous step it means up to now came from a previous run
		#(result of CP2K appending to files rather than overwriting)
		if len(outDict["step"]) > 0:
			currStep = int(currSplitLine[0])
			if currStep <= outDict["step"][-1]:
				outDict["step"] = list()
				outDict["time"] = list()
				outDict["kindTemp"] = [list() for x in range(nKinds)]

		#Parse this line
		outDict["step"].append( int(currSplitLine[0]) )
		outDict["time"].append( float(currSplitLine[1]) )
		for kindIdx in range(nKinds):
			outDict["kindTemp"][kindIdx].append( float(currSplitLine[2+kindIdx]) )


		idx+=1

	return outDict


def _addAtomicTempsToThermoDataObj(atomicTempDict, thermoObj, idxToLabelDict=None):
	idxToLabelDict = dict() if idxToLabelDict is None else idxToLabelDict

	#Step 1 = check step indices are consistent
	stepsThermo, stepsAtomTemp = thermoObj.dataDict["step"], atomicTempDict["step"]
	assert stepsThermo==stepsAtomTemp
	
	#Step 2 = add to thermo obj
	for idx,atomTemps in enumerate(atomicTempDict["kindTemp"]):
		currLabel = idxToLabelDict.get(idx,idx)
		currKey = "kindTemp_{}".format(currLabel)
		thermoObj.dataDict[currKey] = atomTemps

def _getKindIdxToSymbolDict(eleList):
	foundEles = list()
	idx = 0
	outDict = dict()
	for ele in eleList:
		if ele not in foundEles:
			foundEles.append(ele)
			outDict[idx] = ele
			idx+=1
	return outDict



#Metadynamics files below
def parseIterOfMetadynamicsHillsLogFiles(inpPaths):
	""" Runs parseMetadynamicsHillsLogFile on each inpPath and merges the output
	
	Args:
		inpPaths: (iter of str) Iter of *HILLS*.metadynLog paths
			 
	Returns
		outDicts: (iter of dicts) 1 dict per collective variable (in order)

	outDicts:
		keys: "time", "position", "scale", "height"
		values: Each is an iter of values; they should all be the same length. Units are exactly as they are in the output file
 
	WARNING:
		Only tested for up to two collective variables at time of writing
 
	"""
	outDicts = [parseMetadynamicsHillsLogFile(x) for x in inpPaths]
	assert all([len(x)==len(outDicts[0]) for x in outDicts])

	output = list()

	for colVarIdx in range(len(outDicts[0])):
		currDict = outDicts[0][colVarIdx]
		for fIdx in range(1,len(outDicts)):
			thisFileDict = outDicts[fIdx][colVarIdx]
			currDict["time"].extend( thisFileDict["time"] )
			currDict["position"].extend( thisFileDict["position"] )
			currDict["scale"].extend( thisFileDict["scale"] )
			currDict["height"].extend( thisFileDict["height"] )
	
		output.append(currDict)

	return output

def parseMetadynamicsHillsLogFile(inpPath):
	""" Parses basic information from a *HILLS*.metadynLog file.
	
	Args:
		inpPath: (str, file path) A *HILLS*.metadynLog file. Contains information on the hills spawned in an MD simulation
			 
	Returns
		outDicts: (iter of dicts) 1 dict per collective variable (in order)

	outDicts:
		keys: "time", "position", "scale", "height"
		values: Each is an iter of values; they should all be the same length. Units are exactly as they are in the output file
 
	WARNING:
		Only tested for up to two collective variables at time of writing

	"""
	genericArray = np.array(_parseGenericMetadynLogFile(inpPath))
	
	#Figure out how many variables we have
	nCols = len( genericArray[0] )
	nColVars = 1 + int( (nCols-4)/2 )
	nExtraColVars = nColVars-1

	#Parse column by column
	#step 1 - get the shared values
	timesAll = [x[0] for x in genericArray]
	heightsAll = [x[3+2*nExtraColVars] for x in genericArray]

	#Step 2 = parse the specific values and put in the dicts
	outDicts = list()
	for varIdx in range(nColVars):
		colVarIdx = 1+varIdx
		scaleIdx = colVarIdx + nColVars
		currPositions = [x[colVarIdx] for x in genericArray]
		currScales = [x[scaleIdx] for x in genericArray]
		currTimes = [x for x in timesAll]
		currHeights = [x for x in heightsAll]
		currDict = {"time":currTimes, "position":currPositions, "scale":currScales, "height":currHeights}
		outDicts.append(currDict)

	return outDicts

def _parseGenericMetadynLogFile(inpPath):
	fileAsList = _readFileIntoList(inpPath)
	idx = 0
	outArray = list()
	while idx<len(fileAsList):
		splitLine = fileAsList[idx].strip().split()
		if len(splitLine) == 0:
			break
		else:
			outArray.append( [float(x) for x in splitLine] )
		idx += 1

	return outArray

def getGroupedMultiDimGaussHillsFromParsedHillsFileOutput(inpDicts):
	""" Gets a instance representing the combined Gaussian functions from the output of parseMetadynamicsHillsLogFile. Should let the potential be visualised
	
	Args:
		inpDicts: (iter of dicts) Each element corresponds to 1 collective variable. The dicts contain keys for scale/position/height for a set of Gaussian Hills
			 
	Returns
		outHillsFunct: (GroupedMultiDimGaussHills) Instance representing the combination of hills. Has functions to evaluate each individual component AND their rums
 
	"""
	#Step 1 = get MetadynHillsInfo class
	metaDynHillsInstance = getMetadynHillsInfoFromParsedHillsFileOutput(inpDicts)
	return metaDynHillsInstance.createGroupedHills()

#	#Step 1 = get all the 1 dimensional gaussian functions for each collective variable
#	oneDimGauFuncts = list()
#	for cVarIdx in range(len(inpDicts)):
#		heights, positions, scales = [inpDicts[cVarIdx][key] for key in ["height","position","scale"]] 
#		currFuncts = [metadynHillsHelp.OneDimGaussianHill(height=height,pos=pos,scale=scale) for height,pos,scale in it.zip_longest(heights, positions, scales)]
#		oneDimGauFuncts.append( currFuncts )
#
#	#Step 2 = merge the 1 dimensional gaussian functions
#	assert all([len(x)==len(oneDimGauFuncts[0]) for x in oneDimGauFuncts])
#	multiDimFuncts = list()
#	for hillIdx in range(len(oneDimGauFuncts[0])):
#		currOneDims = [ x[hillIdx] for x in oneDimGauFuncts ]
#		multiDimFuncts.append( metadynHillsHelp.MultiDimGaussHill(currOneDims) )
#
#	#Step 3 = combine all the multi dimensional gaussian functions
#	return metadynHillsHelp.GroupedMultiDimGaussHills(multiDimFuncts)


def getMetadynHillsInfoFromParsedHillsFileOutput(inpDicts):
	""" Gets an instance representing combined Gaussian hill functions from the output of parseMetadynamicsHillsLogFile. The class makes it simpler to combine hills and look at potential evolution over time
	
	Args:
		inpDicts: (iter of dicts) Each element corresponds to 1 collective variable. The dicts contain keys for scale/position/height for a set of Gaussian Hills
			 
	Returns
		outObj: (MetadynHillsInfo instance) Data storage class essentially

	"""
	times = inpDicts[0]["time"]
	heights, scales, positions = list(), list(), list()

	for idx in range(len(times)):
		heights.append( [x["height"][idx] for x in inpDicts] ) 
		scales.append( [x["scale"][idx] for x in inpDicts] )
		positions.append( [x["position"][idx] for x in inpDicts] )

	outKwargs = {"times":times, "positions":positions, "scales":scales, "heights":heights}

	return metadynHillsHelp.MetadynHillsInfo(**outKwargs)


