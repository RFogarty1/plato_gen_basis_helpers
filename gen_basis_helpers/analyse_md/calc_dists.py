
import math
import numpy as np

import MDAnalysis.lib.distances as distLib

from . import mdanalysis_interface as mdAnalysisInter

import plato_pylib.shared.ucell_class as uCellHelp


def calcNearestDistanceBetweenTwoSetsOfIndices(inpCell, indicesA=None, indicesB=None, minDist=1e-2):
	""" Calculates the nearest distance between atoms in indicesA and indicesB (uses nearest image convention)
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		indicesA: (Optional, iter of ints) Indices of the atoms to include for the first dimension; Default is to include ALL atoms
		indicesB: (Optional, iter of ints) Indices of the atoms to include for the second dimension; Default is indicesA
		minDist: (float,Optional) Minimum distance. We set this to avoid returning a value of 0 when there is overlap between indicesA and indicesB; if you want it to return 0 in that case then just set minDist to be negative
			 
	Returns
		outDist: (float) The minimum distance between atoms in indicesA and indicesB
 
	"""
	distMatrix = calcDistanceMatrixForCell_minImageConv(inpCell, indicesA=indicesA, indicesB=indicesB)
	for rIdx, unused in enumerate(distMatrix):
		for cIdx, unused in enumerate(distMatrix[rIdx]):
			if distMatrix[rIdx][cIdx] < minDist:
				distMatrix[rIdx][cIdx] = np.inf

	return np.min(distMatrix)


def calcDistanceMatrixForCell_minImageConv(inpCell, indicesA=None, indicesB=None):
	""" Calculates the distance matrix for coords in inpCell using the nearest image convention
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		indicesA: (Optional, iter of ints) Indices of the atoms to include for the first dimension; Default is to include ALL atoms
		indicesB: (Optional, iter of ints) Indices of the atoms to include for the second dimension; Default is indicesA
 
	Returns
		 distMatrix: (NxN numpy array) distMatrix[n][m] gives the distance between atom n and m
 
	"""
	#Sort out default args
	cartCoords = inpCell.cartCoords
	indicesA = [x for x in range(len(cartCoords))] if indicesA is None else indicesA
	indicesB = indicesA if indicesB is None else indicesB

	#Get the coords
	coordsA = np.array([ x[:3] for idx,x in enumerate(inpCell.cartCoords) if idx in indicesA ])
	coordsB = np.array([ x[:3] for idx,x in enumerate(inpCell.cartCoords) if idx in indicesB ])

	#Calculate the relevant distance matrix
	dims = mdAnalysisInter.getMDAnalysisDimsFromUCellObj(inpCell)
	distMatrix = distLib.distance_array(coordsA, coordsB, box=dims)
	return distMatrix

def calcSingleDistBetweenCoords_minImageConv(inpCell, coordA, coordB):
	""" Gets the distance between two co-ordinates in inpCell using the minimum image convention; thus, this will return the smallest distance possible with the PBCs
	
	Args:
		inpCell: (plato_pylib UnitCell object) We use this to get the lattice vectors needed to account for PBCs
		coordA: (len 3 iter) Co-ordinates of first atom
		coordB: (len 3 iter) Co-ordinates of second atom
 
	Returns
		dist: (float) The distance between the co-ordinates after taking PBCs into account
 
	"""
	dims = mdAnalysisInter.getMDAnalysisDimsFromUCellObj(inpCell)
	dists = distLib.calc_bonds(np.array([coordA[:3]]),np.array([coordB[:3]]),box=dims)
	assert len(dists)==1
	return dists[0]


def calcSingleAngleBetweenCoords_minImageConv(inpCell, coordA, coordB, coordC):
	""" Gets the angle (ABC) between three co-ordinates in inpCell using the minimum image convention; thus, this will use the closest possible A,B,C versions
	
	Args:
		inpCell: (plato_pylib UnitCell object) We use this to get the lattice vectors needed to account for PBCs
		coordA: (len 3 iter) Co-ordinates of first atom
		coordB: (len 3 iter) Co-ordinates of second atom
		coordC: (len 3 iter) Co-ordinates of third atom
 
	Returns
		angle: (float) The angle (in degrees) between the co-ordinates after taking PBCs into account
 
	"""
	dims = mdAnalysisInter.getMDAnalysisDimsFromUCellObj(inpCell)
	args = [ np.array([coordA[:3]]), np.array([coordB[:3]]), np.array([coordC[:3]]) ]
	angles = distLib.calc_angles(*args, box=dims)
	assert len(angles)==1
	return math.degrees(angles[0])


#TODO: Not sure if this works for dists >L/2. If not i maybe need to test and raise for it
def getNearestImageNebCoordsBasic(inpCell, coordA, coordB):
	""" Given cartCoordA, find the nearest image co-ordinate for coordB (useful, for example, when a bond crosses periodic boundaries)
	
	Args:
		inpCell: (plato_pylib UnitCell obj) Used to get the latt params/angles
		coordA: (len 3 iter) [x,y,z]. Cartesian co-ordinates expected
		coordB: (len 3 iter) [x,y,z]. Cartesian co-ordinates expected

	Returns
		outCoord: The image of coordB which is nearest to coordA
	
	"""
	lattVects = inpCell.lattVects
	startCartCoords = [coordA[:3], coordB[:3]]
	startFractCoords = [x[:3] for x in uCellHelp.getFractCoordsFromCartCoords(lattVects, startCartCoords)]

	shiftFCoordA = [0.5,0.5,0.5]
	fShiftVector = [x-y for x,y in zip(shiftFCoordA, startFractCoords[0])]
	shiftFCoordB = [x+y for x,y in zip(startFractCoords[1],fShiftVector)]

	#shiftFCoordB = _getFractCoordsFoldedIntoCell(shiftFCoordB)
	shiftFCoordB = _getImageInCellForInpFractCoords(shiftFCoordB)
	diffShiftFractVect = [b-a for b,a in zip(shiftFCoordB,shiftFCoordA)]
	coordBFractVector = [x+y for x,y in zip(diffShiftFractVect, startFractCoords[0])]
	nearestImageCoords = uCellHelp._getCartCoordsFromFract_NoElement(lattVects,coordBFractVector)

	return nearestImageCoords


def _getImageInCellForInpFractCoords(inpCoord, tolerance=1e-2):

	outCoords = list()
	for idx,val in enumerate(inpCoord[:3]):
		currVal = val
		if (val<-1*tolerance):
			shiftVal = math.ceil(abs(val))
			currVal += shiftVal
		elif (val>1+tolerance):
			shiftVal = -1*math.floor(val)
			currVal += shiftVal
		outCoords.append(currVal)

	return outCoords

