
import os
import json

from ..cp2k import parse_pdos_files as parsePdosHelp
from ..misc import dicts_to_db as dictsDbHelp

def dumpPdosFilesFromStdOut_simple(stdOutObj, startDir):
	""" Takes a parsedFile Standard output obj and dumps data into a sub-path based on eleKey/structKey/methodKey of stdOutObj 
	
	Args:
		stdOutObj: Restricted to the "parsedOutFile" type here
		startDir: (str) Path to a base folder for writing database related things. Files will be written to subpaths starting here
			 
	"""
	#Figure out where we dump it
	assert len(stdOutObj.label)==1
	outDict = dict()
	pathDict = getOutPathDictForPdosDump(stdOutObj, startDir)
	outPath = os.path.join(startDir, pathDict["pdos_path_ext"], pathDict["pdos_filename"])

	#Get the dict that we dump
	outDict = dict()
	pdosDict = stdOutObj.data[0][0].parsedFile.pdos
	
	if pdosDict.get("atomKinds", None) is not None:
		outDict["atomKinds"] = _turnListOfPdosObjsToListOfDicts( pdosDict["atomKinds"] )

	if pdosDict.get("atomLists", None) is not None:
		outDict["atomLists"] = _turnListOfPdosObjsToListOfDicts( pdosDict["atomLists"] )

	dictsDbHelp.dumpDictsToFilePath([outDict], outPath)


def _turnListOfPdosObjsToListOfDicts(pdosList):
	return [x.toDict() for x in pdosList]

def getOutPathDictForPdosDump(stdOutObj, startDir):
	""" 
	
	Args:
		stdOutObj: Restricted to the "parsedOutFile" type here
		startDir: (str) Path to a base folder for writing database related things. Files are written to subpaths starting here
			 
	Returns
		outDict: (dict) Contains pathways (relative to a startDir) for dumping nudged band calculation
 
	"""
	assert len(stdOutObj.label)==1
	eleKey = stdOutObj.label[0].eleKey
	structKey = stdOutObj.label[0].structKey
	methodKey = stdOutObj.label[0].methodKey

	pathExt = os.path.join(eleKey, structKey, methodKey)
	outPath = os.path.join(startDir, eleKey, structKey, methodKey, "pdos.json")
	outDict = {"pdos_path_ext":pathExt, "pdos_filename": "pdos.json"}
	return outDict



def getPdosDictFromRecord(startDir, inpRecord):
	""" Standard function to get a pdos dict object from a pdos calculation record
	
	Args:
		startDir: (str) Path to a base folder for writing database related things. Files will be written to subpaths starting here
		inpRecord: (dict) Standard record containing neb pathway, should have partly been generated with dumpPdosFilesFromStdOut_simple

	Returns
		outDict: keys are "atomLists" and "atomKinds". Values are iters of PdosFragmentStandard objects
 
	"""
	inpPath = os.path.join(startDir, inpRecord["pdos_path_ext"], inpRecord["pdos_filename"])
	inpDicts = 	_readSinglePdosIntoListFromJsonFile(inpPath)

	outDicts = dict()
	if inpDicts.get("atomLists",None) is not None:
		outDicts["atomLists"] = [parsePdosHelp.PdosFragmentStandard.fromDict(x) for x in inpDicts["atomLists"]]

	if inpDicts.get("atomKinds",None) is not None:
		outDicts["atomKinds"] = [parsePdosHelp.PdosFragmentStandard.fromDict(x) for x in inpDicts["atomKinds"]]

	return outDicts


def _readSinglePdosIntoListFromJsonFile(inpFile):
	with open(inpFile,"rt") as f:
		outList = json.load(f)
	assert len(outList)==1
	return outList[0]




