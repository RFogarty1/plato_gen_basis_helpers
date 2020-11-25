
import collections

_DATA_PARSABLE_HEADERS = ["LAMMPS Atom File", "Masses", "Atoms", "Bonds", "Angles"]

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
