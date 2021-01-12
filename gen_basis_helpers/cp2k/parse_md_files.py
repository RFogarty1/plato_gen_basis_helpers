
import copy

import plato_pylib.parseOther.parse_cp2k_files as parseCP2KHelp
import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.shared.unit_convs as uConvHelp

#TODO: Relative paths please
import gen_basis_helpers.analyse_md.thermo_data as thermoDataHelp
import gen_basis_helpers.analyse_md.traj_core as trajHelp


def parseFullMdInfoFromCpoutAndXyzFilePaths(cpoutPath, xyzPath):
	outDict = parseCpoutForMDJob(cpoutPath)
	outSteps = parseCp2kMdXyzFile(xyzPath)

	#If step 0 is in the xyz then we need to add the initial trajectory/thermal info to outDict
	if outSteps[0]["step"] == 0:
		stepZero = trajHelp.TrajStepBase(unitCell=outDict["init_md_cell"], step=0, time=0)
		outDict["trajectory"].trajSteps.insert(0, stepZero)
		for key in outDict["thermo_data"].dataDict:
			outDict["thermo_data"].dataDict[key].insert(0, outDict["init_thermo_dict"][key])

	assert len(outSteps)==len(outDict["trajectory"].trajSteps)
	for tStep, xyzStep in zip(outDict["trajectory"].trajSteps, outSteps):
		if tStep.step == xyzStep["step"]:
			tStep.unitCell.cartCoords = xyzStep["coords"]
		else:
			raise ValueError("Step number mismatch between cpout and xyz files: tStep.step = {} but xyzStep[\"step\"] = {}".format(tStep.step,xyzStep["step"]))

	assert outDict["thermo_data"].dataListLengthsAllEqual

	outDict.pop("init_thermo_dict")
	outDict.pop("init_md_cell")

	return outDict

def parseCpoutForMDJob(outFile, parser=None):
	parser = _getStandardCpoutMDParser() if parser is None else parser
	fileAsList = parseCP2KHelp._getFileAsListFromInpFile(outFile)
	outDict = parser.getOutDictFromFileAsList(fileAsList)
	return outDict

def _getStandardCpoutMDParser():
	outObj = parseCP2KHelp.CpoutFileParser()
	outObj = parseCP2KHelp.CpoutFileParser( )
	parseCP2KHelp._addSearchWordAndFunctToParserObj("GO CP2K GO", _parseMDInitSection, outObj, handleParsedDictFunct=_mdInitSectionParseDictHandler)
	parseCP2KHelp._addSearchWordAndFunctToParserObj("ENSEMBLE TYPE                =", _parseMdStepInfo, outObj, handleParsedDictFunct=_mdStepsParseDictHandler)
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
			outDict["lattParams"] = [float(x) for x in splitLine[-3:]]
		if "INITIAL CELL ANGLS[deg]" in currLine:
			splitLine = currLine.strip().split()
			outDict["lattAngles"] = [float(x) for x in splitLine[-3:]]
		if "INITIAL POTENTIAL ENERGY[hartree]" in currLine:
			outDict["ePot"] = float( currLine.strip().split()[-1] )*haToEv
		if "INITIAL KINETIC ENERGY" in currLine:
			outDict["eKinetic"] = float( currLine.strip().split()[-1] )*haToEv
		if "TEMPERATURE" in currLine:
			outDict["temp"] = float( currLine.strip().split()[-1] )

		lineIdx+=1

	return outDict, lineIdx+1


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
		currLine = fileAsList[lineIdx]
		if endStr in currLine:
			break
		if "POTENTIAL ENERGY[hartree]" in currLine:
			outDict["ePot"] = float(currLine.strip().split()[-2])*haToEv
		if "KINETIC ENERGY [hartree]" in currLine:
			outDict["eKinetic"] = float(currLine.strip().split()[-2])*haToEv
		if "PRESSURE [bar]" in currLine:
			outDict["pressure"] = float(currLine.strip().split()[-2])
		if "STEP NUMBER" in currLine:
			outDict["step"] = int(currLine.strip().split()[-1])
		if "TIME [fs]" in currLine:
			outDict["time"] = float(currLine.strip().split()[-1])
		if "TEMPERATURE [K]" in currLine:
			outDict["temp"] = float(currLine.strip().split()[-2])
		if "CELL LNTHS[bohr]             " in currLine:
			outDict["lattParams"] = [float(x) for x in currLine.strip().split()[-3:]]
		if "CELL ANGLS[deg]" in currLine: 
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

