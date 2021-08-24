
import copy
import collections
import itertools as it
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







def calcDistanceMatrixForCell_minImageConv(inpCell, indicesA=None, indicesB=None, sparseMatrix=False):
	""" Calculates the distance matrix for coords in inpCell using the nearest image convention
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		indicesA: (Optional, iter of ints) Indices of the atoms to include for the first dimension; Default is to include ALL atoms
		indicesB: (Optional, iter of ints) Indices of the atoms to include for the second dimension; Default is indicesA
		sparseMatrix: (Optional, Bool) If True the the matrix returned will be NxN (N being number of atoms in inpCell) even when indicesA or indicesB set. In that case we just set the unwanted indices to np.nan

	Returns
		 distMatrix: (NxN numpy array) distMatrix[n][m] gives the distance between atom n and m
 
	"""
	#Sort out default args
	cartCoords = inpCell.cartCoords
	indicesA = [x for x in range(len(cartCoords))] if indicesA is None else indicesA
	indicesB = indicesA if indicesB is None else indicesB

	#Edge case: Empty list passed for one of the indices
	if (indicesA==list()) or (indicesB==list()):
		if sparseMatrix:
			outMatrix = np.empty( (len(cartCoords),len(cartCoords)) )
			outMatrix[:] = np.nan
			return outMatrix

	#Get the coords
	coordsA = np.array( [cartCoords[idx][:3] for idx in indicesA] )
	coordsB = np.array( [cartCoords[idx][:3] for idx in indicesB] )

	#Calculate the relevant distance matrix
	dims = mdAnalysisInter.getMDAnalysisDimsFromUCellObj(inpCell)
	distMatrix = distLib.distance_array(coordsA, coordsB, box=dims)
	if sparseMatrix:
		outDim = len(cartCoords)
		distMatrix = _getTwoDimSparseMatrix(distMatrix, outDim, indicesA, indicesB)

	return distMatrix


 

#NOTE: I could probably extend this to a different plane; but suspect it would need to contain at least one cell vector
#(and maybe even two)
def calcHozDistMatrixForCell_minImageConv(inpCell, indicesA=None, indicesB=None, minTotInterPlaneDist=1e-5, sparseMatrix=False):
	""" Calculates matrix of horizontal distances (i.e. distance along surface plane) for inpCell
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		indicesA: (Optional, iter of ints) Indices of the atoms to include for the first dimension; Default is to include ALL atoms
		indicesB: (Optional, iter of ints) Indices of the atoms to include for the second dimension; Default is indicesA
		minTotInterPlaneDist: (float) We calculate hoz-distance by using totalDist-interPlaneDist. If these values are the same (e.g. hozDist=0 then float errors may make totalDist-interPlaneDist negative which leads to a domain error when square-rooting. minTotInterPlaneDist means to set hozDist to zero in this case
		sparseMatrix: (Optional, Bool) If True the the matrix returned will be NxN (N being number of atoms in inpCell) even when indicesA or indicesB set. In that case we just set the unwanted indices to np.nan
 
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

	#Return sparse matrix if requested
	if sparseMatrix:
		outDim = len(cartCoords)
		outMatrix = _getTwoDimSparseMatrix(outMatrix, outDim, indicesA, indicesB)

	return outMatrix


#Only sensible for dist matrices function really
#NOTE: Theres more copy call than i need...
class _SpecialMemoizationFunct():

	def __init__(self, funct, maxEntries=20):
		self.argCache = collections.deque(list(), maxlen=maxEntries)
		self.funct = funct

	def __call__(self, *args, **kwargs):

		#1) check if args already in cache. 
		currArgs = tuple([x for x in args])
		sortedKwargs = sorted(kwargs.items())
		for prevArgs,prevKwargs,output in self.argCache:
			if (prevArgs==currArgs) and (prevKwargs==sortedKwargs):
				return np.copy(output) #NEED TO RETURN A COPY, so theres no way to corrupt the memoized

		#2)If not, run the function and add to the cache
		output = self.funct(*args,**kwargs)
		toAppend =  [tuple([x for x in args]),sorted(kwargs.items()),np.copy(output)]
		self.argCache.appendleft( toAppend )

		return np.copy(output)

_calcDistanceMatrixForCell_minImageConv_memoized = _SpecialMemoizationFunct(calcDistanceMatrixForCell_minImageConv)
_calcHozDistMatrixForCell_minImageConv_memoized = _SpecialMemoizationFunct(calcHozDistMatrixForCell_minImageConv)

#This lets me insert a pdb at each call; which lets me track the cache more easily
def _calcHozDistMatrixForCell_minImageConv_memoized_debuggable(*args, **kwargs):
	return _calcHozDistMatrixForCell_minImageConv_memoized(*args, **kwargs)


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


def calcDistancesFromSurfPlaneForCell(inpCell, indices=None, planeEqn=None, sparseMatrix=False):
	""" Calculates distances of atoms from a surface plane for inpCell
	
	Args:
		inpCell: (plato_pylib UnitCell object) 
		indices: (iter of ints) Indices of atoms to calculate this distance for. Default is to use all indices in the cell
		planeEqn: (ThreeDimPlaneEquation object). MUST be parralel to axb; where axb are cell vector
		sparseMatrix: (Optional, Bool) If True the the matrix returned will be len-N (N being number of atoms in inpCell) even when indices is set. In that case we just set the unwanted indices to np.nan
			 
	Returns
		outDists: (iter of floats; same length as indices) Each is the distance from the plane
 
	"""
	#Create a cell to use for this
	useCell = copy.deepcopy(inpCell)
	useCartCoords = useCell.cartCoords

	#Deal with default args	
	indices = [idx for idx in range(len(useCartCoords))] if indices is None else indices
	planeEqn = cartHelp.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(useCell) if planeEqn is None else planeEqn

	#Add an atom in the input planeEqn
	newCoord = planeEqn.getPointClosestToOrigin() + ["X"]
	useCartCoords.append(newCoord)
	useCell.cartCoords = useCartCoords

	#Get the nearest coords to the surface plane
	indicesForAllAtoms = [idx for idx in range(len(useCartCoords))]
	coordMatrix = getNearestImageNebCoordsMatrixBasic(useCell, indicesA=[indicesForAllAtoms[-1]], indicesB=indicesForAllAtoms)

	#Get the distances
	outDists = list()
	for idx in indices:
		coordA, coordB = coordMatrix[0][idx]
		currDist = planeEqnHelp.getInterPlaneDistTwoPoints( coordA[:3], coordB[:3], planeEqn )
		outDists.append(currDist)

	#Convert output to sparse form if requested
	if sparseMatrix:
		outDists = _getOneDimSparseMatrix(outDists, len(inpCell.cartCoords), indices)
		outDists = outDists.tolist()

	return outDists


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



def getNearestImageVectorsForIdxPairs(inpCell, idxPairs, sparseMatrix=False):
	""" Gets the vectors [coordB - coordA] (i.e. the position vector to go from A to B) for each idx pair [A,B] using the minimum image convention. This can be MUCH faster than using getNearestImageVectorMatrixBasic if you only want a subset of vectors calculated
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		idxPairs: (iter of ints) The indices we want the position vectors for. If sparseMatrix output is used then the reverse vectors should automatically be calculated
		sparseMatrix: (Bool) If True return a sparsely populated matrix of values (explained more below)
 
	Returns
		outVects: one vector per idxPair if sparseMatrix is False, else NxNx3 matrix containing np.nan for most entries but [coordB-coordA] for those corresponding to values in idxPairs. Here "N" refers to the number of atoms in inpCell
 
	"""
	#Step 1: Get the MINIMAL nearest image matrix by doing it row by row
	useCell = uCellHelp.UnitCell(lattParams=inpCell.getLattParamsList(), lattAngles=inpCell.getLattAnglesList())
	fractCoords = inpCell.fractCoords
	fCoordsNoEles = np.array([x[:3] for x in fractCoords])
	eleList = [x[-1] for x in fractCoords]

	outDim = len(fractCoords)
#	sparseNearestNebMatrix = [ [list() for x in range(outDim)]  for unused in range(outDim) ]

	sparseNearestNebMatrix = np.empty( (outDim,outDim,2,3) )
	sparseNearestNebMatrix[:] = np.nan

	sortedIdxPairs = sorted(idxPairs)

	#Step 1.X: Loop over idx pairs, calculating the relevant row for sparseNearestNebMatrix
	idx = 0
	while idx < len(sortedIdxPairs):
		#Get all relevant idxB values for this row
		startIdx, rowIdx = idx, sortedIdxPairs[idx][0]
		endIdx = startIdx + 1

		while endIdx < len(sortedIdxPairs):
			if sortedIdxPairs[endIdx][0] != rowIdx:
				break
			else:
				endIdx += 1

		indicesB = [x[1] for x in idxPairs[startIdx:endIdx]]

		#Append to the sparse matrix and Update our idx value
		currRowVals = _getSingleRowOfNearestImageNebCoordsMatrix(fCoordsNoEles, eleList, useCell, rowIdx, indicesB)
		for idxB, currRowVal in it.zip_longest(indicesB, currRowVals):
			sparseNearestNebMatrix[rowIdx][idxB][0] = currRowVal[0][:3] #Dont need the element type
			sparseNearestNebMatrix[rowIdx][idxB][1] = currRowVal[1][:3]

		idx = endIdx

	#Step 2: Build sparse matrix with the vectors
	outSparseMatrix = np.empty( (outDim,outDim,3) )
	outSparseMatrix[:] = np.nan

	for idxPair in idxPairs:
		coordA, coordB = sparseNearestNebMatrix[idxPair[0]][idxPair[1]]
		outSparseMatrix[ tuple(idxPair)] = np.array([b-a for b,a in it.zip_longest(coordB[:3],coordA[:3])])
		outSparseMatrix[ tuple(reversed(idxPair)) ] = np.array([a-b for b,a in it.zip_longest(coordB[:3],coordA[:3])])

	#Step 3: Also build a sparse output matrix as we go; and optionally return it
	if sparseMatrix:
		outVals = outSparseMatrix
	else:
		outVals =  [ outSparseMatrix[tuple(idxPair)] for idxPair in idxPairs ] 

	return outVals


def getNearestImageVectorMatrixBasic(inpCell, indicesA=None, indicesB=None, sparseMatrix=False):
	""" Gets a matrix where each element contains [coordB-coordA] where A,B are row/column indices. coordA is that found in inpCell, coordB is the nearest neighbour image. Thus function gets the vectors between two atoms using minimum image convention
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		indicesA: (Optional, iter of ints) Indices of the atoms to include for the first dimension; Default is to include ALL atoms
		indicesB: (Optional, iter of ints) Indices of the atoms to include for the second dimension; Default is indicesA
		sparseMatrix: (Optional, Bool) If True the the matrix returned will be NxN (N being number of atoms in inpCell) even when indicesA or indicesB set. In that case we just set the unwanted indices to np.nan
			 
	Returns
		outMatrix: (NxM Matrix) outMatrix[n][m] = [coordM-coordN] where N and M are lengths of indicesA and indicesB
 
	"""
	outMatrix = getNearestImageNebCoordsMatrixBasic(inpCell, indicesA=indicesA, indicesB=indicesB)

	for rIdx in range(len(outMatrix)):
		for cIdx in range(len(outMatrix[rIdx])):
			coordA, coordB = outMatrix[rIdx][cIdx]
			outMatrix[rIdx][cIdx] = [b-a for b,a in it.zip_longest(coordB[:3],coordA[:3])]
#			outMatrix[cIdx][rIdx] = [a-b for b,a in it.zip_longest(coordB[:3],coordA[:3])]

	if sparseMatrix:
		cartCoords = inpCell.cartCoords
		outDim = len(cartCoords)
		outMatrix = _getTwoDimSparsePosVectorMatrix(np.array(outMatrix), outDim, indicesA, indicesB)

	return outMatrix

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
	fCoordsNoEles = np.array([x[:3] for x in fractCoords])
	eleList = [x[-1] for x in fractCoords]
	indicesA = [x for x in range(len(fractCoords))] if indicesA is None else indicesA
	indicesB = indicesA if indicesB is None else indicesB

	#Figure out the value for each row; we need useCell to exploit some relevant function in the ucell_class thing
	useCell = uCellHelp.UnitCell(lattParams=inpCell.getLattParamsList(), lattAngles=inpCell.getLattAnglesList())
	outMatrix = list()
	for idxA in indicesA:
		currRow = _getSingleRowOfNearestImageNebCoordsMatrix(fCoordsNoEles, eleList, useCell, idxA, indicesB)
		outMatrix.append(currRow)

	return outMatrix


def _getSingleRowOfNearestImageNebCoordsMatrix(startFractCoords, eleList, useCell, idxA, indicesB):
	#Setup the cell to have the fractCoords we want
#	startFractCoords = inpCell.fractCoords
	centralFractUnshifted = startFractCoords[idxA]
	fractCoords = np.array( [centralFractUnshifted.tolist()] + [startFractCoords[idx].tolist() for idx in indicesB] )

	#Shift the cell such that the first index is right at the centre + fold all atoms into the cell
	centralFractShifted = [0.5,0.5,0.5]
	fractShiftVector = [x-y for x,y in zip(centralFractShifted, centralFractUnshifted)]
	reverseFractShiftVector = [-1*x for x in fractShiftVector]

	useFractCoords = np.add(fractCoords, fractShiftVector)
	uCellHelp.foldFractCoordArrayToValsBetweenZeroAndOne(useFractCoords, tolerance=1e-5)

	#Now shift it back; this now contains only the nearest image co-ords for the index 0 val
	useFractCoords = np.add(useFractCoords, reverseFractShiftVector)

	#Convert to cartesian co-ordinates and 
	outCartCoords = uCellHelp._getCartCoordsFromFract_NoElement(useCell.lattVects, useFractCoords)

	outRow = list()
	for idxB,coord in enumerate(outCartCoords[1:]):
		partA = [x for x in outCartCoords[0]] + [eleList[idxA]] #Can reach pretty significant runtime portion
		partB = coord + [eleList[indicesB[idxB]]]
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
#	useCell = copy.deepcopy(inpCell) 
	useCell = uCellHelp.UnitCell(lattParams=inpCell.getLattParamsList(), lattAngles=inpCell.getLattAnglesList())
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


def getInterAtomicAnglesForInpGeom(inpCell, angleIndices, degrees=True):
	""" Gets a list of requested angles between atoms for inpCell taking PBCs into account.
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		angleIndices: (iter of len-3 iters) Each element contains [idxA,idxB,idxC] which means calculate the angle between idA,idxB,idxC for indices in inpCell.cartCoords
		degrees: (Bool) If True then return angles in degrees; else use radians

	Returns
		outAngles: (iter of floats) Length is the same as len(angleIndices); each is an angle calculated for a value in angleIndices
 
	"""
	#1) Check indices not empty; if empty just return empty list
	if len(angleIndices)==0:
		return list()

	cartCoords = inpCell.cartCoords

	#Convert to format needed for mdAnalysis function
	numbAngles = len(angleIndices)
	cartCoordsA = np.array([ np.array(cartCoords[angleIndices[idx][0]][:3]) for idx in range(numbAngles) ])
	cartCoordsB = np.array([ np.array(cartCoords[angleIndices[idx][1]][:3]) for idx in range(numbAngles) ])
	cartCoordsC = np.array([ np.array(cartCoords[angleIndices[idx][2]][:3]) for idx in range(numbAngles) ])

	#Calculate using md analysis
	dims = mdAnalysisInter.getMDAnalysisDimsFromUCellObj(inpCell)
	outAnglesRadians = distLib.calc_angles(cartCoordsA, cartCoordsB, cartCoordsC, box=dims)
	outAnglesRadians = [x for x in outAnglesRadians]

	#Convert output to the form i want
	if degrees:
		outAngles = [math.degrees(x) for x in outAnglesRadians]
	else:
		outAngles = outAnglesRadians

	return outAngles



def _getTwoDimSparsePosVectorMatrix(inpMatrix, outDim, indicesA, indicesB):
	""" Gets a sparsely populated 2-D vector matrix from a given inpMatrix. This allows us to not worry about how atom indices map to indicesA and indicesB, while still not calculating more values than we need. This differs SLIGHTLY from getting a sparse 3-d matrix in that we assume here that the vectors are len-3 by default (may extend later)
	
	Args:
		inpMatrix: (NxMx3 matrix) N and M are dimensions of indicesA and indicesB. inpMatrix[1][2] contains value between indicesA[1] and indicesB[2]
		outDim: (int) The length of one side of the output (square) matrix
		indicesA: (iter of ints) Indices used to build inpMatrix rows
		indicesB: (iter of ints) Indices used to build inpMatrix columns
 
	Returns
		outMatrix: (NxNx3 matrix) Where N is outDim. outMatrix[idxA][idxB] gets the value between idxA and idxB; if the information wasnt in inpMatrix then all elements of the vector will be np.nan.
 
	"""
	outMatrix = np.empty( (outDim,outDim,3) )
	outMatrix[:] = np.nan

	#Should probably factor this out; almost identical to the normal 2-dim case
	allInpIndices = ( tuple([idxA,idxB,idxC]) for idxA,idxB,idxC in it.product( range(len(indicesA)), range(len(indicesB)), range(3) ) )
	allOutIndices = ( tuple([idxA,idxB,idxC]) for idxA,idxB,idxC in it.product( indicesA, indicesB, range(3) ) )

	for inpIdx, outIdx in it.zip_longest(allInpIndices, allOutIndices):
		outMatrix[ tuple(outIdx) ] = inpMatrix[ tuple(inpIdx) ]
		outMatrix[ tuple( [x for x in reversed(outIdx[:2])] + [outIdx[2]]) ] = -1*inpMatrix[ tuple(inpIdx) ]

	return outMatrix

def _getTwoDimSparseMatrix(inpMatrix, outDim, indicesA, indicesB):
	""" Gets a sparsely populated 2-D matrix from a given inpMatrix. This allows us to not worry about how atom indices map to indicesA and indicesB, while still not calculating more values than we need
	
	Args:
		inpMatrix: (NxM matrix) N and M are dimensions of indicesA and indicesB. inpMatrix[1][2] contains value between indicesA[1] and indicesB[2]
		outDim: (int) The length of one side of the output (square) matrix
		indicesA: (iter of ints) Indices used to build inpMatrix rows
		indicesB: (iter of ints) Indices used to build inpMatrix columns
 
	Returns
		outMatrix: (NxN matrix) Where N is outDim. outMatrix[idxA][idxB] gets the value between idxA and idxB; if the information wasnt in inpMatrix then the value will be np.nan.
 
	"""
	outMatrix = np.empty( (outDim,outDim) )
	outMatrix[:] = np.nan

	#this is a smallish bottleneck; but unsure if theres any faster way to do it
	allInpIndices = ( tuple([idxA,idxB]) for idxA,idxB in it.product( range(len(indicesA)), range(len(indicesB)) ) )
	allOutIndices = ( tuple([idxA,idxB]) for idxA,idxB in it.product( indicesA, indicesB ) )

	for inpIdx,outIdx in it.zip_longest(allInpIndices,allOutIndices):
		outMatrix[ tuple(outIdx) ] = inpMatrix[ tuple(inpIdx) ]
		outMatrix[ tuple(reversed(outIdx)) ] = inpMatrix[ tuple(inpIdx) ]

	return outMatrix


def _getOneDimSparseMatrix(inpMatrix, outDim, indices):
	""" Get a sparsely populated 1-d matrix from a given 1-d matrix. 
	
	Args:
		inpMatrix: (len-N matrix, or iter of floats) N is len(indices) in this case
		outDim: (int) The length of one side of the output matrix; should be number of atoms in general
		indices: (iter of ints) Indices used to build inpMatrix
 
	Returns
		outMatrix: (len-N numpy array) Where N is outDim. outMatrix[idx] gets the value for a given atom index; if the information wasnt in inpMatrix then the value will be np.nan 
 
	"""
	outMatrix = np.empty( (outDim) )
	outMatrix[:] = np.nan

	inpIndices = [ idx for idx in range(len(indices)) ]
	outIndices = [ idx for idx in indices ]

	for inpIdx,outIdx in it.zip_longest(inpIndices, outIndices):
		outMatrix[outIdx] = inpMatrix[inpIdx]

	return outMatrix




