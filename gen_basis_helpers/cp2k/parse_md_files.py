
import copy
import itertools as it

import plato_pylib.parseOther.parse_cp2k_files as parseCP2KHelp
import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.shared.unit_convs as uConvHelp

#TODO: Relative paths please
import gen_basis_helpers.analyse_md.thermo_data as thermoDataHelp
import gen_basis_helpers.analyse_md.traj_core as trajHelp


def parseMdInfoFromMultipleCpoutAndXyzPaths(cpoutPaths, xyzPaths):
	outDict = dict()

	#1) Get all the dicts
	parsedDicts = list()
	for cpoutPath, xyzPath in it.zip_longest(cpoutPaths, xyzPaths):
		currDict = parseFullMdInfoFromCpoutAndXyzFilePaths(cpoutPath, xyzPath)
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

def parseFullMdInfoFromCpoutAndXyzFilePaths(cpoutPath, xyzPath):
	outDict = parseCpoutForMDJob(cpoutPath)
	outSteps = parseCp2kMdXyzFile(xyzPath)

	#If step 0 is in the xyz then we need to add the initial trajectory/thermal info to outDict
	if outSteps[0]["step"] == 0:
		stepZero = trajHelp.TrajStepBase(unitCell=outDict["init_md_cell"], step=0, time=0)
		outDict["trajectory"].trajSteps.insert(0, stepZero)
		for key in outDict["thermo_data"].dataDict:
			outDict["thermo_data"].dataDict[key].insert(0, outDict["init_thermo_dict"][key])


	outDict["trajectory"] = _getMergedTrajectoryFromParsedCpoutAndXyz(outDict, outSteps)

	assert outDict["thermo_data"].dataListLengthsAllEqual

	outDict.pop("init_thermo_dict")
	outDict.pop("init_md_cell")

	return outDict

#Also modifies in place
#Relies on step numbers being in order for both parsedCpout and parsedXyz
def _getMergedTrajectoryFromParsedCpoutAndXyz(parsedCpout, parsedXyz):
	trajCpout, trajXyz = parsedCpout["trajectory"].trajSteps, parsedXyz
	outList = list()

	assert len(trajXyz)<=len(trajCpout)

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
#	import pdb
#	pdb.set_trace()

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
		currStep = trajHelp.TrajStepBase(unitCell=geom, step=thermoInfo["step"][idx], time=thermoInfo["time"][idx])
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

