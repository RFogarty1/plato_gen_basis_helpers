
import plato_pylib.shared.ucell_class as uCell

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
	



