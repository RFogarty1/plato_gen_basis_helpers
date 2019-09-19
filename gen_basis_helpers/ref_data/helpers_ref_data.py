
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
	casOutFiles = _getCastepOutPathsForFolder(refFolder)
	parsedUCells = [parseCastep.parseCastepOutfile(x)["unitCell"] for x in casOutFiles]
	[x.convAngToBohr() for x in parsedUCells]
	return parsedUCells


def getEosFitDictFromEosCastepFolder(refFolder,eos="murnaghan"):
	outPaths = _getCastepOutPathsForFolder(refFolder)
	return _getEosDictFromFilePaths(outPaths,eos)

def _getEosDictFromFilePaths(filePaths,eos):
	outDict = fitBMod.getBulkModFromOutFilesAseWrapper(filePaths, eos=eos)
	return outDict

def _getCastepOutPathsForFolder(refFolder):
	return [os.path.join(refFolder,x) for x in os.listdir(refFolder) if x.endswith('.castep')]


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

