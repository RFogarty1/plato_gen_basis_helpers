
import itertools as it
import os
import math

import plato_pylib.parseOther.parse_castep_files as parseCastep
import plato_pylib.shared.ucell_class as UCell
import plato_pylib.utils.fit_eos as fitBMod


def getPerfectHcpMinimalUCell(element):
	lattVects = [ [ 1.00000000, 0.00000000, 0.00000000],
	              [-0.50000000, 0.86602540, 0.00000000],
	              [ 0.00000000, 0.00000000, 1.00000000] ]
	
	fractCoords = [ [0.0,0.0,0.0,element], [ 1/3, 2/3, 0.5, element] ]
	idealCoverA = math.sqrt( 8/3 )

	perfectHcpCell = UCell.UnitCell.fromLattVects(lattVects, fractCoords=fractCoords)
	perfectHcpCell.setLattParams([1.0,1.0,idealCoverA])
	return perfectHcpCell

def getUCellsFromCastepBulkModFolder(refFolder):
	casOutFiles = getCastepOutPathsForFolder(refFolder)
	parsedUCells = [parseCastep.parseCastepOutfile(x)["unitCell"] for x in casOutFiles]
	[x.convAngToBohr() for x in parsedUCells]
	return parsedUCells


def getEosFitDictFromEosCastepFolder(refFolder,eos="murnaghan"):
	outPaths = getCastepOutPathsForFolder(refFolder)
	return _getEosDictFromFilePaths(outPaths,eos)


def getUCellInBohrFromCastepOutFile(outFilePath):
	uCell = parseCastep.parseCastepOutfile(outFilePath)["unitCell"]
	uCell.convAngToBohr()
	return uCell

def _getEosDictFromFilePaths(filePaths,eos):
	outDict = fitBMod.getBulkModFromOutFilesAseWrapper(filePaths, eos=eos)
	return outDict

def getCastepOutPathsForFolder(refFolder):
	return [os.path.join(refFolder,x) for x in os.listdir(refFolder) if x.endswith('.castep')]





def getDimerSepFromCastepOutFile(inpFile,inBohr=True):
	""" Gets separation
	
	Args:
	  inpFile(str): Path to *.castep output file
	  inBohr(opt, bool): If true convert distance from angstrom to bohr; if False use angstrom
 
	Returns
		 dimerSep: Separation between the two atoms in the castep file
 
	Raises:
		 AssertionError: If the last unit-cell in the input file doesnt have exactly 2 atoms
	"""
	parsedUCell = parseCastep.parseCastepOutfile(inpFile)["unitCell"]
	if inBohr:
		parsedUCell.convAngToBohr()
	return getDimerSepFromUCellObj(parsedUCell)


def getDimerSepFromUCellObj(uCellObj):
	""" Returns distance of a dimer in UnitCell format (see plato_pylib.shared.ucell_class)
	
	Args:
	 uCellObj: UnitCell object from plato_pylib representing geometry of a dimer.
			 
	Returns
		 dist: Separation (in same units as UnitCell obj) between the two atoms 
 
	Raises:
		 AssertionError: If UnitCell does not contain exactly 2 atoms
	"""

	cartCoords = uCellObj.cartCoords
	assert len(cartCoords)==2, "UnitCell must contain 2 atoms, input cell contained {}".format( len(cartCoords) )
	atomA, atomB = cartCoords[0][:3], cartCoords[1][:3]
	return _getDistTwoVectors(atomA,atomB)


def _getDistTwoVectors(vectA, vectB):
	diffVector = [a-b for a,b in it.zip_longest(vectA,vectB)]
	return math.sqrt( sum([x**2 for x in diffVector]) ) 




#Likely not gonna be used again, may remove if so
def _getUCellFromCrystalMakerCastepOutFile(refFile):
	tokenizedFile = parseCastep.tokenizeCastepCellFileAndRemoveBlockFromKeys(refFile)

	#Step 1 = get lattice vectors via the cell parameters. The fract co-ords should work fine for this (checked for octa supercell of hcp)
	cellParamStrList = tokenizedFile["lattice_abc"].strip().split()
	lattParams = [float(x) for x in cellParamStrList[1:4]]
	lattAngles = [float(x) for x in cellParamStrList[4:]]
	outUCell = UCell.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	outUCell.convAngToBohr()

	#Step 2 = get the fractional co-ords. Note I'm pretty sure this relies on me picking the "right" lattice vectors.
	outUCell.fractCoords = parseCastep._getFractCoordsFromTokenizedCellFile(tokenizedFile)

	return outUCell

