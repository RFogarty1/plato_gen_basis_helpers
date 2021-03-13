
import os


#import gen_basis_helpers.cp2k.parse_neb_files as parseNebHelp

from ..cp2k import parse_neb_files as parseNebHelp
from ..misc import shared_io as sharedIoHelp
from ..misc import nudged_band_paths as pathwayHelp

def dumpNebFilesFromStdOut_simple(stdOutObj, startDir):
	""" Takes a parsedFile Standard output obj and dumps nudged elastic band pathway (only file for now) into a sub-path based on eleKey/structKey/methodKey of stdOutObj 
	
	Args:
		stdOutObj: Restricted to the "parsedOutFile" type here
		startDir: (str) Path to a base folder for writing database related things. Files will be written to subpaths starting here
			 
	"""
	assert len(stdOutObj.label)==1

	pathDict = getOutPathDictForNebFileDump(stdOutObj, startDir)
	outPath = os.path.join(startDir, pathDict["neb_path_ext"], pathDict["neb_pathway_filename"])
	nebPathway = parseNebHelp.getNebPathFromParsedFileObj( stdOutObj.data[0][0].parsedFile )

	sharedIoHelp.dumpObjWithToDictToJson(nebPathway, outPath)

def getOutPathDictForNebFileDump(stdOutObj, startDir):
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
	outPath = os.path.join(startDir, eleKey, structKey, methodKey, "neb_path.json")
	outDict = {"neb_path_ext":pathExt, "neb_pathway_filename": "neb_path.json"}
	return outDict

def getNebPathwayFromRecord(startDir,inpRecord):
	""" Standard function to get a NudgedBandPathStandard object from a nudged band calculation record
	
	Args:
		startDir: (str) Path to a base folder for writing database related things. Files will be written to subpaths starting here
		inpRecord: (dict) Standard record containing neb pathway, should have partly been generated with dumpNebFilesFromStdOut_simple

	Returns
		pathway: (NudgedBandPathStandard) Contains object representing a nudged elastic band pathway
 
	"""
	outPath = _getNebPathwayFilePathFromRecord(startDir, inpRecord)
	return pathwayHelp.readNudgedBandPathwayFromJsonFileStandard(outPath)


def _getNebPathwayFilePathFromRecord(startDir, inpRecord):
	outPath = os.path.join(startDir, inpRecord["neb_path_ext"], inpRecord["neb_pathway_filename"])
	return outPath

