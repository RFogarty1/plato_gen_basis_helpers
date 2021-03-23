
import itertools as it
import os
import re

import plato_pylib.shared.unit_convs as uConvHelp


def parsePdosFromCpoutPath(cpoutPath):
	""" Parses partial/projected density of states (which cp2k dumps in special files) from a cpoutPath (we dont actually use that cpoutPath, just uses the folder to figure it out)
	
	Args:
		cpoutPath: (str) Path to the *.cpout file
			 
	Returns
		outDict: (dict) Keys are "atomKinds" and "listKinds". Values are iter of PdosFragmentStandard
 
	"""
	outDict = dict()
	outDict["atomKinds"] = _parsePdosKindFilesFromCpoutPath(cpoutPath)
	outDict["atomLists"] = _parsePdosListFilesFromCpoutPath(cpoutPath)
	return outDict


def _parsePdosKindFilesFromCpoutPath(cpoutPath):
	inpPaths = getPdosKindsPathsFromCpoutPath(cpoutPath)
	outPdos = [parsePdosFromFile(inpPath) for inpPath in inpPaths]
	return outPdos

def _parsePdosListFilesFromCpoutPath(cpoutPath):
	inpPaths = getPdosAtomicListsPathsFromCpoutPath(cpoutPath)
	outPdos = [parsePdosFromFile(inpPath) for inpPath in inpPaths]
	return outPdos

def getPdosKindsPathsFromCpoutPath(cpoutPath):
	""" Gets a list of paths containing pdos data for different atomic kinds from path to the main output (cpout) file
	
	Args:
		cpoutPath: (str) Path to the *.cpout file
			 
	Returns
		outPaths: (iter of str) Paths to any files contaninig pdos data for atomic kinds
 
	"""
	basePath = os.path.splitext(cpoutPath)[0]
	folderPath = os.path.split(basePath)[0]

	pattern = "k[0-9]+-1.pdos"
	relPaths = list()
	for outName in sorted( os.listdir( folderPath ) ):
		currMatches = re.findall(pattern, outName)
		if len(currMatches)==1:
			relPaths.append( os.path.join(folderPath,outName) )
	return sorted(relPaths)


def getPdosAtomicListsPathsFromCpoutPath(cpoutPath):
	""" Gets a list of paths containing pdos data for different atomic lists from path to the main output (cpout) file
	
	Args:
		cpoutPath: (str) Path to the *.cpout file
			 
	Returns
		outPaths: (iter of str) Paths to any files contaninig pdos data for atomic kinds
 
	"""
	basePath = os.path.splitext(cpoutPath)[0]
	folderPath = os.path.split(basePath)[0]
	
	pattern = "-list[0-9]+-1.pdos"
	relPaths = list()
	for outName in sorted( os.listdir( folderPath ) ):
		currMatches = re.findall(pattern, outName)
		if len(currMatches)==1:
			relPaths.append( os.path.join(folderPath,outName) )
	return sorted(relPaths)


class PdosFragmentStandard():

	def __init__(self, eigenValues=None, occs=None, fragName=None,
	             breakdowns=None, breakdownHeaders=None):
		""" Iniitalizer
		
		Args:
			eigenValues: (iter of floats) Each element is eigenvalue for one MO
			occs: (iter of floats) Each element is occupancy for one MO
			fragName: (str, or None) Name of fragment. e.g. "O" for the oxygen pdos
			breakdowns: (iter of iter of floats) Each element corresponds to 1 eigenvalue. Each entry (in 1 element) corresponds to the fractional contribution from a shell (or similar)
			breakdownHeaders: (iter of str) labels for the breakdowns in breakdownHeaders
 
		"""
		self._eqTol = 1e-5
		self.eigenValues = eigenValues
		self.occs = occs
		self.fragName = fragName
		self.breakdowns = breakdowns
		self.breakdownHeaders = breakdownHeaders

	@classmethod
	def fromDict(cls,inpDict):
		return cls(**inpDict)

	def toDict(self):
		outKeys = ["eigenValues", "occs", "fragName", "breakdowns", "breakdownHeaders"]
		return {k:getattr(self,k) for k in outKeys}


	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)

		directCmpAttrs = ["fragName", "breakdownHeaders"]
		num1dArrayAttrs = ["eigenValues", "occs"]
		num2dArrayAttrs = ["breakdowns"]

		#
		for attr in directCmpAttrs:
			if getattr(self, attr) != getattr(other, attr):
				return False

		#eigenvalues/occupancies
		for attr in num1dArrayAttrs:
			ours, theirs = getattr(self, attr), getattr(other, attr)
			if (ours is None) and (theirs is None):
				pass
			elif (ours is None) or (theirs is None):
				return False
			else:
				if len(ours) != len(theirs):
					return False
				for valA, valB in it.zip_longest(ours, theirs):
					if abs(valA-valB)>eqTol:
						return False

		#breakdowns
		for attr in num2dArrayAttrs:
			ours, theirs = getattr(self, attr), getattr(other, attr)

			if (ours is None) and (theirs is None):
				pass
			elif (ours is None) or (theirs is None):
				return False
			else:
				if len(ours) != len(theirs):
					return False

				for arrA, arrB in it.zip_longest(ours,theirs):
					if len(arrA) != len(arrB):
						return False
					for valA, valB in it.zip_longest(arrA, arrB):
						if abs(valA-valB) > eqTol:
							return False

		return True


def parsePdosFromFile(inpFile):
	""" Parses a single partial(projected) density of states from a cp2k format file (*.pdos)
	
	Args:
		inpFile: (str) Path to the file containing a pdos
			 
	Returns
		outPdos: (PdosFragmentStandard) Object contaninig all info on the pdos
 
	"""
	fileAsList = _readFileIntoList(inpFile)

	eigenVals, occs = list(), list()
	breakdowns = list()
	moIndices = list()
	lineIdx = 0

	#Check for an atomic kind
	atKindMatch = re.findall("atomic kind [A-Z,a-z]+ ", fileAsList[0]) 
	if len(atKindMatch)>0:
		assert len(atKindMatch)==1
		fragName = atKindMatch[0].strip().split()[-1]
	else:
		fragName = None

	#Get the breakdown headers
	numLine = fileAsList[2].strip().split()
	nBreakdown = len(numLine) - 3
	breakdownHeaders = fileAsList[1].strip().split()[-nBreakdown:]

	lineIdx = 2
	haToEv = uConvHelp.RYD_TO_EV*2
	while lineIdx < len(fileAsList):
		currLine = fileAsList[lineIdx].strip().split()
		moIndices.append(currLine[0])
		eigenVals.append( float(currLine[1])*haToEv )
		occs.append( float(currLine[2]) )
		breakdowns.append( [float(x) for x in currLine[3:]] )
		lineIdx += 1

	#Paranoid error check
	assert moIndices == sorted(moIndices)

	#Create the output object
	currKwargs = {"eigenValues":eigenVals, "occs":occs, "fragName":fragName,
	              "breakdowns":breakdowns, "breakdownHeaders":breakdownHeaders}

	return PdosFragmentStandard(**currKwargs)


def _readFileIntoList(inpFile):
	with open(inpFile, "rt") as f:
		outList = f.readlines()
	return outList











