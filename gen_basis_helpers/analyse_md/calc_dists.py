
import copy
import math
import numpy as np

import MDAnalysis.lib.distances as distLib

from . import mdanalysis_interface as mdAnalysisInter

from ..shared import plane_equations as planeEqnHelp
from ..shared import cart_coord_utils as cartHelp

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
	coordsA = np.array( [cartCoords[idx][:3] for idx in indicesA] )
	coordsB = np.array( [cartCoords[idx][:3] for idx in indicesB] )

	#Calculate the relevant distance matrix
	dims = mdAnalysisInter.getMDAnalysisDimsFromUCellObj(inpCell)
	distMatrix = distLib.distance_array(coordsA, coordsB, box=dims)
	return distMatrix


#NOTE: I could probably extend this to a different plane; but suspect it would need to contain at least one cell vector
#(and maybe even two)
def calcHozDistMatrixForCell_minImageConv(inpCell, indicesA=None, indicesB=None, minTotInterPlaneDist=1e-5):
	""" Calculates matrix of horizontal distances (i.e. distance along surface plane) for inpCell
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		indicesA: (Optional, iter of ints) Indices of the atoms to include for the first dimension; Default is to include ALL atoms
		indicesB: (Optional, iter of ints) Indices of the atoms to include for the second dimension; Default is indicesA
		minTotInterPlaneDist: (float) We calculate hoz-distance by using totalDist-interPlaneDist. If these values are the same (e.g. hozDist=0 then float errors may make totalDist-interPlaneDist negative which leads to a domain error when square-rooting. minTotInterPlaneDist means to set hozDist to zero in this case
 
	Returns
		hozDistMatrix:  (NxM numpy array) distMatrix[n][m] gives the distance between atom n and m
 
	"""

	#Sort out default args
	cartCoords = inpCell.cartCoords
	indicesA = [x for x in range(len(cartCoords))] if indicesA is None else indicesA
	indicesB = indicesA if indicesB is None else indicesB

	#1) Get dist matrix
	#2) Get interplanar (PBCs wont affect this if we use the surface plane specifically)
	#3) hozDit = math.sqrt(totalDist**2 - interPlanarDist**2)
	totalDistMatrix = calcDistanceMatrixForCell_minImageConv(inpCell, indicesA=indicesA, indicesB=indicesB) #Somehow this seems to come out wrong
	surfPlane = cartHelp.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(inpCell)
	nearestImageNebCoordsMatrix = getNearestImageNebCoordsMatrixBasic(inpCell, indicesA=indicesA, indicesB=indicesB)



	outMatrix = np.zeros( [len(indicesA), len(indicesB)] )
	for dMatrixIdxA, atomIdxA in enumerate(indicesA):
		for dMatrixIdxB,atomIdxB in enumerate(indicesB):
			currTotalDist = totalDistMatrix[dMatrixIdxA][dMatrixIdxB]

			posA, posB = [x[:3] for x in nearestImageNebCoordsMatrix[dMatrixIdxA][dMatrixIdxB]] #Taking from here effectively handles PBCs
			currInterPlaneDist = planeEqnHelp.getInterPlaneDistTwoPoints(posA, posB, surfPlane)

			#Figure out the horizontal distance
			if abs(currTotalDist-currInterPlaneDist) < minTotInterPlaneDist:
				outMatrix[dMatrixIdxA][dMatrixIdxB] = 0
			else:
				outMatrix[dMatrixIdxA][dMatrixIdxB] = math.sqrt(currTotalDist**2 - currInterPlaneDist**2)

	return outMatrix

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



def getInterSurfPlaneSeparationTwoPositions(posA, posB, inpCell):
	""" Gets the minimum separation along the surface normal for two (cartesian) positios when given inpCell. The surface plane is defined by axb
	
	Args:
		posA: (len-3 float iter) [x,y,z]
		posB: (len-3 float iter) [x,y,z]
		inpCell: (plato_pylib UnitCell object). Used to get the a and b lattice vectors
 
	Returns
		outDist: (float) Absolute distance between the surface planes containing posA and posB [these planes are parralel to axb but have different d values when a plane is defined as ax + by +cz = d
 
	"""
	surfPlaneEqn = cartHelp.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(inpCell)
	nearestImageCoords = getNearestImageNebCoordsBasic(inpCell, posA, posB)
	return planeEqnHelp.getInterPlaneDistTwoPoints(posA, nearestImageCoords, surfPlaneEqn)


def getNearestImageNebCoordsMatrixBasic(inpCell, indicesA=None, indicesB=None):
	""" Gets a matrix where each element contains [coordA,coordB] where A,B are row/column indices. coordA is the same co-ordinate found in inpCell; coordB is the nearest neighbour for it
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		indicesA: (Optional, iter of ints) Indices of the atoms to include for the first dimension; Default is to include ALL atoms
		indicesB: (Optional, iter of ints) Indices of the atoms to include for the second dimension; Default is indicesA
			
	Returns
		outMatrix: (NxM matrix) outMatrix[n][m] = [coordN, coordM] where N and M are lengths of indicesA and indicesB
 
	"""
	#Sort out default args
	fractCoords = inpCell.fractCoords
	indicesA = [x for x in range(len(cartCoords))] if indicesA is None else indicesA
	indicesB = indicesA if indicesB is None else indicesB

	#Figure out the value for each row; we need useCell to exploit some relevant function in the ucell_class thing
	useCell = uCellHelp.UnitCell(lattParams=inpCell.getLattParamsList(), lattAngles=inpCell.getLattAnglesList())
	outMatrix = list()
	for idxA in indicesA:
		currRow = _getSingleRowOfNearestImageNebCoordsMatrix(fractCoords, useCell, idxA, indicesB)
		outMatrix.append(currRow)

	return outMatrix


def _getSingleRowOfNearestImageNebCoordsMatrix(startFractCoords, useCell, idxA, indicesB):
	#Setup the cell to have the fractCoords we want
#	startFractCoords = inpCell.fractCoords
	centralFractUnshifted = startFractCoords[idxA]
	fractCoords = [centralFractUnshifted] + [startFractCoords[idx] for idx in indicesB]
	useCell.fractCoords = fractCoords

	#Shift the cell such that the first index is right at the centre + fold all atoms into the cell
	centralFractShifted = [0.5,0.5,0.5]
	fractShiftVector = [x-y for x,y in zip(centralFractShifted, centralFractUnshifted)]
	reverseFractShiftVector = [-1*x for x in fractShiftVector]
	uCellHelp.applyTranslationVectorToFractionalCoords(useCell, fractShiftVector, foldInAfter=False)
	uCellHelp.foldAtomicPositionsIntoCell(useCell, tolerance=1e-5)

	#Now shift it back; this now contains only the nearest image co-ords for the index 0 val
	uCellHelp.applyTranslationVectorToFractionalCoords(useCell, reverseFractShiftVector, foldInAfter=False)

	outRow = list()
#	outFractCoords = useCell.fractCoords
	outCartCoords = useCell.cartCoords
	for coord in outCartCoords[1:]:
		partA = [x for x in outCartCoords[0]]
		partB = coord
		outRow.append( [partA,partB] )

	return outRow

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
	indicesA, indicesB = [0],[1]
	useCell = copy.deepcopy(inpCell)
	useCell.cartCoords = [ coordA[:3] + ["X"], coordB[:3] + ["X"] ] 
	outMatrix = getNearestImageNebCoordsMatrixBasic(useCell, indicesA=indicesA, indicesB=indicesB)
	return outMatrix[0][0][1][:3]


#TODO: Remove; it was previously used in getNearestImageNebCoordsBasic but isnt anymore
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

