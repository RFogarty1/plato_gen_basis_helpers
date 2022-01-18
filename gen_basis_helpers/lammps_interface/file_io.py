
import collections

import plato_pylib.shared.ucell_class as uCellHelp

_DATA_PARSABLE_HEADERS = ["LAMMPS Atom File", "Masses", "Atoms", "Bonds", "Angles"]


def getUCellObjectFromDataFile(inpPath, atomStyle="full", massDict=None):
	massDict = massDict if massDict is not None else uCellHelp.getEleKeyToMassDictStandard()

	tokenizedFile = tokenizeDataFile(inpPath)
	#Step 1 = get lattice parameters
	outCell = _getUnitCellObjFromDataHeaderSection(tokenizedFile["LAMMPS Atom File"])

	#Step 2 = get the atomic co-ordinates; this depends on atomStyle
	if atomStyle=="full":
		typeIdxToEle = _getTypeIdxToEleFromMassesSectionAndMassDict(tokenizedFile["Masses"],massDict)
		cartCoords = _getCartCoordsFromAtomsSection_atomStyleFull(tokenizedFile["Atoms"], typeIdxToEle)
		outCell.cartCoords = cartCoords
	else:
		raise NotImplementedError("atomStyle={} is not implemented yet".format(atomStyle))

	return outCell

def _getUnitCellObjFromDataHeaderSection(headerStr):
	sectionAsList = [x.strip() for x in headerStr.split("\n")]
	tiltFactors = [0,0,0]
	diagParams = [0,0,0]

	for line in sectionAsList:
		#Get diag elements of lattice vectors
		if ("xlo" in line) or ("ylo" in line) or ("zlo" in line):
			splitLine = line.split() #Defaults to splitting all whitespace
			currParam = float(splitLine[1]) - float(splitLine[0])
			if "xlo" in line:
				currIdx = 0
			elif "ylo" in line:
				currIdx = 1
			elif "zlo" in line:
				currIdx = 2
			else:
				raise ValueError("No idea how loop ended up here - it shouldnt be possible")
			diagParams[currIdx] = currParam

		#Deal with tilt factors
		if "xy xz yz" in line:
			tiltFactors = [float(x) for x in line.split()[:3]]

	lattVectors = [ [diagParams[0], 0, 0],
	                [tiltFactors[0], diagParams[1], 0],
	                [tiltFactors[1], tiltFactors[2], diagParams[2]] ]
	return uCellHelp.UnitCell.fromLattVects(lattVectors)
	

def _getTypeIdxToEleFromMassesSectionAndMassDict(massesStr, massDict, massTol=1e-2):
	sectionAsList = [x.strip() for x in massesStr.split("\n")]
	outDict = dict()
	for line in sectionAsList:
		if line!="":
			splitLine = line.split()
			idx, mass = int(splitLine[0]), float(splitLine[1])
			currKeys = list()
			for k,v in massDict.items():
				if abs(float(v)-mass) < massTol:
					currKeys.append(k)
			assert len(currKeys)==1, "Found {} key for mass = {}".format(len(currKeys), mass)
			outDict[idx] = currKeys[0]

	return outDict

def _getCartCoordsFromAtomsSection_atomStyleFull(atomsStr, idxToEle):
	sectionAsList = [x.strip() for x in atomsStr.split("\n")]
	
	outCoords = list()
	for line in sectionAsList:
		if line!="":
			splitLine = line.split()
			currEle = idxToEle[ int(splitLine[2]) ]
			currXYZ = [float(x) for x in splitLine[-3:]]
			currCoord = currXYZ + [currEle]
			outCoords.append(currCoord)
	return outCoords

def tokenizeDataFile(inpPath):
	""" Takes a LAMMPS data file and converts it into an OrderedDict
	
	Args:
		inpPath: (str) Path to the input file
			 
	Returns
		outDict: (OrderedDict) Keys are headers, vals are strs for that section
 
	"""
	fileAsList = _getFileAsListFromInpPath(inpPath)
	outDict = collections.OrderedDict()
	currIdx = 0

	while currIdx<len(fileAsList):
		for pattern in _DATA_PARSABLE_HEADERS:
			if pattern in fileAsList[currIdx]:
				unused, currStr = _parseOneSection(fileAsList, currIdx)
				outDict[pattern] = currStr
		currIdx += 1

	return outDict


def writeDataFileFromTokens(outPath, tokens):
	""" Writes a set of tokens to an output file in format of LAMMPS data file
	
	Args:
		outPath: (str) Path to the output file we write to
		tokens: (OrderedDict) Keys are headers, vals are the string to write
 
	Returns
		Nothing but will create the file (assuming the directory exists)
	 
	"""
	outStr = ""
	for key in tokens.keys():
		outStr += key.rstrip() + "\n"
		outStr += "\n" #Blank line between header and body
		outStr += tokens[key].rstrip() + "\n"
		outStr += "\n" #Blank link between all body text regions
	_writeFileFromStr(outPath,outStr)
	


def writeScriptFileFromTokens(outPath, tokens):
	""" Writes a set of tokens to an output file in the format of a LAMMPS script file
	
	Args:
		outPath: (str) Path to the output file we write to
		tokens: (OrderedDict) Keys are commands, vals are the arguments passed to commands
 
	Returns
		 Noting but will create the file (assuming the directory exists)
 
	"""
	outStr = ""
	for key in tokens.keys():
		outStr += key.rstrip() + " "
		outStr += tokens[key].strip() + "\n"
	_writeFileFromStr(outPath,outStr)

#NOTE: This may not work for some files that i dont generate. I think blank lines are probably allowed within any section (their DEFO allowed within the header part)
def _parseOneSection(fileAsList, startIdx):
	outIdx = startIdx+2 #We dont want the final blank line
	outListEntries = list()

	while outIdx<len(fileAsList):
		currLine = fileAsList[outIdx]
		if currLine.strip()=="":
			break
		else:
			outListEntries.append(currLine)
			outIdx+=1

	return outIdx,"\n".join(outListEntries)



def _getFileAsListFromInpPath(inpPath):
	outStr = _getStrFromInpPath(inpPath)
	return outStr.split("\n")


def _getStrFromInpPath(inpPath):
	with open(inpPath,"rt") as f:
		outStr = f.read()
	return outStr

def _writeFileFromStr(inpPath, fileStr):
	with open(inpPath,"wt") as f:
		outStr = f.write(fileStr)
	return outStr
