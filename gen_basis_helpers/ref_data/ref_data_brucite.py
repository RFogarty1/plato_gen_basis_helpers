
import copy
import os

from ..castep import castep_creator as castepCreator
from ..shared import config_vars as configVars


import plato_pylib.shared.ucell_class as uCell


BASE_FOLDER = os.path.join( configVars.CASTEP_DB_PATH, "brucite" )


def getExptBruciteStructure():
	""" Data taken from  https://doi.org/10.1007/BF00202300(Catti 1995, lattice paramters and z-coords) and https://doi.org/10.5006/1322 (remaning fractional co-ordinates) """

	zHydrogen = 0.413
	zOxygen   = 0.220
	fractCoords = [ [0  , 0   ,0          ,"Mg"],
	                [1/3, 2/3 ,    zOxygen, "O"],
	                [2/3, 1/3,   1-zOxygen, "O"],
	                [1/3, 2/3,   zHydrogen, "H"],
	                [2/3, 1/3, 1-zHydrogen, "H"] ]
	lattParamA = 3.15
	lattParamC = 4.77
	lattParams = [lattParamA, lattParamA, lattParamC]
	lattAngles = [90, 90, 120]

	outUCell = uCell.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	outUCell.fractCoords = fractCoords
	outUCell.convAngToBohr()
	return outUCell
	

def getStructsForEos(structType="hexagonal"):
	if structType!="hexagonal":
		raise ValueError("Only hexagonal cells implemented at the moment")

	#About 5% in each direction around the experimental value. A bit further for the high volume range though (since PW overestimate volume)
	aVals = [3.005, 3.029, 3.053, 3.077, 3.102, 3.126, 3.150,
	         3.174, 3.198, 3.223, 3.247, 3.271, 3.295,
	         3.319, 3.343, 3.367, 3.391] 

	startCell = getPlaneWaveGeom()
	cOverA = startCell.lattParams["c"]/startCell.lattParams["a"]

	allOutCells = list()
	for a in aVals:
		currLattParams = [a,a,cOverA*a]
		currCell = copy.deepcopy(startCell)
		currCell.setLattParams(currLattParams)
		currCell.convAngToBohr()
		allOutCells.append(currCell)

	return allOutCells


#Plane wave geometry
def getPlaneWaveGeomParsedFileObject(structType="hexagonal"):
	structTypeToFunct = {"hexagonal":_getPlaneWaveHexagonalParsedFileObj}
	return structTypeToFunct[structType.lower()]()

def _getPlaneWaveHexagonalParsedFileObj():
	refPath = os.path.join(BASE_FOLDER,"opt_geom","hexagonal","geom_opt.castep")
	return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)

def getPlaneWaveGeom(structType="hexagonal"):
	return getPlaneWaveGeomParsedFileObject(structType).unitCell



