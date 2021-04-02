
import copy
import math
import itertools as it

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.utils.supercell as supCellHelp

from . import plane_equations as planeEqnHelp
from . import simple_vector_maths as vectHelp


def getAveragePlaneEqnForAtomIndices(inpGeom, indices, planeEqn):
	""" Returns a planeEqn which uses the mean d-value for the atom indices supplied.
	
	Args:
		inpGeom: (plato_pylib UnitCell object)
		indices: (iter of ints) Indices of the atoms to use
		planeEqn: (ThreeDimPlaneEquation) Defines the direction of the plane (only a,b,c need to be set really)
			 
	Returns
		outPlaneEqn: (ThreeDimPlaneEquation) 
 
	Raises:
		 Errors
	"""
	#Get all the d values
	dVals = list()
	for idx,coord in enumerate(inpGeom.cartCoords):
		if idx in indices:
			currDVal = planeEqn.calcDForInpXyz(coord[:3])
			dVals.append(currDVal)

	#Get the average value + return in form of a new plane equation
	outDVal = sum(dVals) / len(dVals)
	outCoeffs = planeEqn.coeffs[:3] + [outDVal]
	return planeEqnHelp.ThreeDimPlaneEquation(*outCoeffs)


def getMostCentralIdxFromList(inpGeom, indices):
	""" Gets the atom index which is closest to the centre of inpGeom when given a list of atom indices
	
	Args:
		inpGeom: (plato_pylib UnitCell object)
		indices: (iter of ints) Indices of atoms; we restrict the search to these indices

	Returns
		outIdx: (int) Index corresponding to the most central atom
	
	"""
	centralFractCoord = [0.5,0.5,0.5,"X"]
	centralCartCoord = uCellHelp.getCartCoordsFromFractCoords( inpGeom.lattVects,[centralFractCoord])[0][:3]
	otherPoints = [ inpGeom.cartCoords[idx][:3] for idx in indices ] 
	return indices[getIdxOfNearestPointToInputPoint(centralCartCoord, otherPoints)]


def getDistancesOfAtomsFromPlaneEquation_nearestImageAware(inpGeom, planeEqn, atomIndices, signed=False):
	""" Gets distances of atoms from an input plane
	
	Args:
		inpGeom: (plato_pylib UnitCell object)
		planeEqn: (ThreeDimPlaneEquation) The plane we want distances for
		atomIndices: (iter of int) 
		signed: (Bool) If True get signed distances, if False get unsigned distances
			 
	Returns
		outDists: (iter of float) Distances for each atomIndex from the plane
 
	"""
	
	startPlaneDVal = planeEqn.coeffs[-1]
	#0) Get the initial cartestian co-ordinates
	cartCoords = inpGeom.cartCoords
	tempCell = copy.deepcopy(inpGeom)

	lattParams, lattAngles = inpGeom.getLattParamsList() , inpGeom.getLattAnglesList()
	tempCellToGetNewPlaneEqn = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)

	#1) Figure out the shift we need to get the planeEqn that passes through [0.5,0.5,0.5]
	tempCellToGetNewPlaneEqn.fractCoords = [ [0.5,0.5,0.5,"X"] ]
	dValForCentreOfCell = planeEqn.calcDForInpXyz( tempCellToGetNewPlaneEqn.cartCoords[0][:3] )
	planeEqnForCentreOfCell = copy.deepcopy(planeEqn)
	planeEqnForCentreOfCell = planeEqnHelp.ThreeDimPlaneEquation( *(planeEqn.coeffs[:3] + [dValForCentreOfCell]) )
	shiftPlaneToCentreVector = planeEqnHelp.getVectorToMoveFromParallelPlanesAToB(planeEqn,planeEqnForCentreOfCell)

	#2) Translate the cell to get our plane equation at the centre of it
	tempCell = copy.deepcopy(inpGeom)
	newCellCartCoords = list()
	for coord in cartCoords:
		newCoord = [x+t for x,t in it.zip_longest(coord[:3], shiftPlaneToCentreVector)]
		newCellCartCoords.append( newCoord + [coord[-1]] )

	tempCell.cartCoords = newCellCartCoords
	uCellHelp.foldAtomicPositionsIntoCell(tempCell)

	#3) Get distances
	useCoords = tempCell.cartCoords
	outDists = list()
	for idx in atomIndices:
		if signed:
			currDist = planeEqnForCentreOfCell.getSignedDistanceOfPointFromPlane( useCoords[idx][:3] )
		else:
			currDist = planeEqnForCentreOfCell.getDistanceOfPointFromPlane( useCoords[idx][:3] )
		outDists.append(currDist)

	return outDists


def getHeightOfCell_abSurface(inpCell):
	""" Gets the height of a cell with the surface defined in the ab plane. For an orthogonal cell this is just c, else its the distance between top/bottom planes in the cell (<c) for other cells
	
	Args:
		inpCell: (UnitCell object) Contains lattice parameters and angles
			 
	Returns
		 height: (float) Height of the cell
 
	"""

	aVect,bVect,cVect = inpCell.lattVects
	pointOnTopPlane = cVect
	bottomPlane = planeEqnHelp.ThreeDimPlaneEquation.fromTwoPositionVectors(aVect,bVect)
	return bottomPlane.getDistanceOfPointFromPlane(pointOnTopPlane)


def shiftCoordsToLeaveEqualVacAboveAndBelowSurface(inpCell):
	""" Shift cartesian co-ordinates such that top/bottom surface atoms are equal distances from top/bottom of the cell. Note that the surface is defined as the ab plane
	
	Args:
		inpCell: (UnitCell object)
 
	"""

	#Get points on the top/bottom planes of the cell
	aVect, bVect, cVect = inpCell.lattVects
	pointOnBotPlane = aVect
	pointOnTopPlane = [a+c for a,c in zip(aVect,cVect)]

	#Get plane equations for top/bottom atom planes
	planeEqnTopAtoms = getPlaneEqnForOuterSurfaceAtoms(inpCell)
	planeEqnBotAtoms = getPlaneEqnForOuterSurfaceAtoms(inpCell, top=False)

	#Check that no atoms are above/below the cell
	signedDistTopAtomToTopCellPlane = planeEqnTopAtoms.getSignedDistanceOfPointFromPlane(pointOnTopPlane)
	signedDistBotAtomToBotCellPlane = planeEqnBotAtoms.getSignedDistanceOfPointFromPlane(pointOnBotPlane)
	errorTol = 1e-2
	if signedDistTopAtomToTopCellPlane<( -1*errorTol ):
		raise ValueError("Atoms in inpCell must be contained within cell boundaries")
	if signedDistBotAtomToBotCellPlane<( -1*errorTol ):
		raise ValueError("Atoms in inpCell must be contained within cell boundaries")

	#Get the amount of vacuum at the top/bottom of the cell to start with
	vacTop = planeEqnTopAtoms.getDistanceOfPointFromPlane(pointOnTopPlane)
	vacBot = planeEqnBotAtoms.getDistanceOfPointFromPlane(pointOnBotPlane)
	averageVac = (vacTop+vacBot) / 2
	vacDiff = averageVac - vacTop

	#Figure out the translation vector that gets that diff
	#-1*vacDiff is because vacDiff works on logic of shiting vacuum, while this function
	#works on logic of moving surface effectively
	topSurfNormal = planeEqnTopAtoms.coeffs[:3]
	transVector = getCVectorCorrespondingToSurfHeightShift(-1*vacDiff, topSurfNormal, cVect)

	#Apply the translation to all cart coords
	cartCoords = inpCell.cartCoords
	outCartCoords = list()
	for coord in cartCoords:
		currCoord = [x+t for x,t in it.zip_longest(coord[:3],transVector)] + [coord[-1]]
		outCartCoords.append(currCoord)

	inpCell.cartCoords = outCartCoords


#Effectively a backend function, but does get called outside this module so no _
def getCVectorCorrespondingToSurfHeightShift(heightShift, surfNormal, cVector):
	theta = vectHelp.getAngleTwoVectors(surfNormal, cVector)
	unitCVector = vectHelp.getUnitVectorFromInpVector(cVector)
	lenC = heightShift / math.cos(math.radians(theta))
	outVector = [lenC*u for u in unitCVector]
	return outVector


def getDistancesFromPointAlongPlane(inpPoint, otherPoints, planeEqn):
	""" Get "horizontal" distances (distance along a plane) of a set of points from an input point. By along plane i mean the second term in tVector = tVect(between-planes) + tVect(within-plane). If your plane equation happens to be [0,0,1] = N then this will correspond to the separation along the xy directions

	Args:
		inpPoint: [x,y,z]
		otherPoints: iter of [x,y,z]
		planeEqn: (ThreeDimPlaneEquation)

	Returns
		outDists: (iter of floats) Distances along the plane
 
	"""
	return [planeEqnHelp.getOutOfPlaneDistTwoPoints(inpPoint, x, planeEqn) for x in otherPoints]

#Functions including interfaces to unitCell class
def getNearestInPlaneDistanceGivenInpCellAndAtomIdx(inpCell, atomIdx, planeEqn, includeImages=True, planeTolerance=1e-2):
	""" Returns the nearest in-plane neighbour distance for atomIdx.
	
	Args:
		inpCell: plato_pylib UnitCell object
		atomIdx: index of atom in the unitCell
		planeEqn: ThreeDimPlaneEquation object, defines the plane to look in. Doesnt need to pass through atomIdx; this function will generate a shifted plane which DOES pass through it regardless
		includeImages: Bool (Optional), If True then periodic images (in all directions) will be considered
			 
	Returns
		 nebDist: The distance to the nearest in plane neighbour
 
	"""
	if includeImages:
		startCell = supCellHelp.getUnitCellSurroundedByNeighbourCells(inpCell)
	else:
		startCell = inpCell

	cartCoords = copy.deepcopy([x[:3] for x in startCell.cartCoords])
	inpPoint = cartCoords[atomIdx]
	cartCoords.pop(atomIdx)
	planeCoeffs = planeEqn.coeffs
	planeDValue = planeEqn.calcDForInpXyz(inpPoint)
	planeCoeffs[-1] = planeDValue

	newPlane = planeEqnHelp.ThreeDimPlaneEquation(*planeCoeffs)

	return getDistanceToNearestInPlanePointToInpPoint(inpPoint, cartCoords, newPlane, planeTolerance)


def getPlaneEqnForOuterSurfaceAtoms(inpCell, top=True):
	""" Gets plane equation for the outer atoms in the xy plane, the plane eqn should point towards vacuum/out from the surface at least
	
	Args:
		inpCell: (plato_pylib UnitCell object) 
		top: (Bool, default=True) If true the surface will be that furthest along c (the "top" surface), else it will be the surface most -ve along c (the bottom surface)
 
	Returns
		 planeEqn: (ThreeDimPlaneEquation object) outermost surface plane containing atoms
 
	"""
	abPlaneEqn = getABPlaneEqnWithNormVectorSameDirAsC(inpCell.lattVects)

	#Find the plane-equations for the top or bottom plane [ASSUMES PERFECTLY FLAT SURFACE for now]
	#(Stolen from the generic surface code for now; TODO: Factor this out)
	allDVals = list()
	for x in inpCell.cartCoords:
		allDVals.append( abPlaneEqn.calcDForInpXyz(x[:3]) )
	maxD, minD = max(allDVals), min(allDVals)
	
	if top:
		outPlaneEquation = planeEqnHelp.ThreeDimPlaneEquation( *(abPlaneEqn.coeffs[:3] + [maxD]) )
	else:
		outPlaneEquation = planeEqnHelp.ThreeDimPlaneEquation( *(abPlaneEqn.coeffs[:3] + [minD]) )
		outPlaneEquation.coeffs = [x*-1 for x in outPlaneEquation.coeffs]

	return outPlaneEquation


def getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(inpCell):
	""" See getABPlaneEqnWithNormVectorSameDirAsC docstring """
	return getABPlaneEqnWithNormVectorSameDirAsC(inpCell.lattVects)

def getABPlaneEqnWithNormVectorSameDirAsC(lattVects):
	""" Gets a plane equation for the a,b plane with the normal vector pointing along the direction of the c vector. Useful for defining "above" or "below" in a consistent way
	
	Args:
		lattVects: len 3 iter of len 3 iters, lattice vectors a,b and c.
			 
	Returns
		 planeEqn: ThreeDimPlaneEquation object, represents the plane equation for the ab plane whereby the normal vector points along the same rough direction as the c-vector (dot product is positive)
 
	"""
	aVect, bVect, cVect = lattVects
	startPlaneEqn = planeEqnHelp.ThreeDimPlaneEquation.fromTwoPositionVectors(aVect,bVect)
	normVector = startPlaneEqn.coeffs[:3]

	if (vectHelp.getDotProductTwoVectors(normVector,cVect) < 0):
		outCoeffs = [-1*x for x in startPlaneEqn.coeffs]
	else:
		outCoeffs = startPlaneEqn.coeffs

	return planeEqnHelp.ThreeDimPlaneEquation(*outCoeffs)


#Functions for getting inter-planar separation
def getUniquePlaneDistsAndAtomIndicesFromSurfacePlaneAndCartCoords(inpPlaneEqn, cartCoords, planeTolerance):
	idxVsDist = _getAtomIndicesVsSignedDistsFromPlane(inpPlaneEqn, cartCoords)
	return _getUniquePlaneDistsAndAtomIndicesFromIdxVsDistList(idxVsDist, planeTolerance)

def _getAtomIndicesVsSignedDistsFromPlane(inpPlaneEqn, cartCoords):
	idxVsDist = list()
	#Get the distances of each atom from the bottom/top plane of the surface
	for idx,coord in enumerate(cartCoords):
		currDist = inpPlaneEqn.getSignedDistanceOfPointFromPlane(coord[:3])
		idxVsDist.append( (idx,currDist) )

	return idxVsDist

def _getUniquePlaneDistsAndAtomIndicesFromIdxVsDistList(idxVsDist, planeTolerance):
	allDists = [ idxVsDist[0][1] ]
	uniquePlaneDists = [ idxVsDist[0][1] ]
	uniquePlaneAtomIndices = [ [idxVsDist[0][0]] ]
	for idx,dist in idxVsDist[1:]:
		sortedIntoPlane = False
		#Check if this plane has already been found
		for pIdx, pDist in enumerate(uniquePlaneDists):
			if abs(pDist-dist) < planeTolerance:
				uniquePlaneAtomIndices[pIdx].append(idx)
				sortedIntoPlane = True

		if not sortedIntoPlane:
			uniquePlaneDists.append( dist )
			uniquePlaneAtomIndices.append( [idx] ) 

	return uniquePlaneDists, uniquePlaneAtomIndices

#Functions for getting in-plane OR out-of plane nearest neighbours
def getDistanceToNearestInPlanePointToInpPoint(inpPoint, otherPoints, planeEqn, planeTolerance=1e-2):
	outCoords = getCoordsOfNearestInPlanePointToInpPoint(inpPoint, otherPoints, planeEqn, planeTolerance)
	outDist = vectHelp.getDistTwoVectors(inpPoint,outCoords)
	return outDist

def getDistanceToNearestOutOfPlanePointToInpPoint(inpPoint, otherPoints, planeEqn, planeTolerance=1e-2):
	outCoords = getCoordsOfNearestOutOfPlanePointToInpPoint(inpPoint, otherPoints, planeEqn, planeTolerance)
	outDist = vectHelp.getDistTwoVectors(inpPoint,outCoords)
	return outDist

def getCoordsOfNearestInPlanePointToInpPoint(inpPoint, otherPoints, planeEqn, planeTolerance=1e-2):
	outIdx = getIdxOfNearestInPlanePointToInpPoint(inpPoint, otherPoints, planeEqn, planeTolerance)
	return otherPoints[outIdx]

def getCoordsOfNearestOutOfPlanePointToInpPoint(inpPoint, otherPoints, planeEqn, planeTolerance=1e-2):
	outIdx = getIdxOfNearestOutOfPlanePointToInpPoint(inpPoint, otherPoints, planeEqn, planeTolerance)
	return otherPoints[outIdx]

def getIdxOfNearestOutOfPlanePointToInpPoint(inpPoint, otherPoints, planeEqn, planeTolerance=1e-2):
	filterFunct = getFilteredIndicesForCoordsOutOfInputPlane
	outIdx = _getIdxOfNearestPointToInpPointAfterFiltering(inpPoint, otherPoints, planeEqn, filterFunct, planeTolerance)
	return outIdx

def getIdxOfNearestInPlanePointToInpPoint(inpPoint, otherPoints, planeEqn, planeTolerance=1e-2):
	filterFunct = getFilteredIndicesForCoordsInInputPlane
	outIdx = _getIdxOfNearestPointToInpPointAfterFiltering(inpPoint, otherPoints, planeEqn, filterFunct, planeTolerance)
	return outIdx

def _getIdxOfNearestPointToInpPointAfterFiltering(inpPoint, otherPoints, planeEqn, filterFunct, planeTolerance=1e-2):
	relevantIndices = filterFunct(otherPoints, planeEqn, planeTolerance)
	relevantCoords = [otherPoints[idx] for idx in relevantIndices]
	idxInRelevantCoords = getIdxOfNearestPointToInputPoint(inpPoint, relevantCoords)
	outIdx = relevantIndices[idxInRelevantCoords]
	return outIdx

def getIndicesNearestInPlanePointsUpToN(nPoints, inpPoint, otherPoints, planeEqn, planeTolerance=1e-2):
	inPlaneIndices = getFilteredIndicesForCoordsInInputPlane(otherPoints, planeEqn, planeTolerance)
	inPlaneCoords = [otherPoints[idx] for idx in inPlaneIndices]
	dists = [vectHelp.getDistTwoVectors(inpPoint, x) for x in inPlaneCoords]
	idxVsDist = [ [idx,x] for idx,x in it.zip_longest(inPlaneIndices,dists)]
	sortedIndices = [ x[0] for x in sorted( idxVsDist, key=lambda x:x[1]) ]
	if nPoints<len(sortedIndices):
		outIndices = sortedIndices[:nPoints]
	else:
		outIndices = sortedIndices
	return outIndices


#Functions for filtering co-ordinates based on plane equations
def getFilteredIndicesForCoordsOutOfInputPlane(inpCoords, planeEqn, planeTolerance=1e-2):
	return _getFilteredIndicesBasedOnWhetherTheyAreInPlane(inpCoords, planeEqn, False, planeTolerance)

def getFilteredIndicesForCoordsInInputPlane(inpCoords, planeEqn, planeTolerance=1e-2):
	return _getFilteredIndicesBasedOnWhetherTheyAreInPlane(inpCoords, planeEqn, True, planeTolerance)

def _getFilteredIndicesBasedOnWhetherTheyAreInPlane(inpCoords, planeEqn, keepInPlane, planeTolerance=1e-2):
	""" Get a set of co-ordinates from inpCoords based on whether they lie on the input plane
	
	Args:
		inpCoords: iter of iters(len-3), x,y,z co-ordinates for a set of points
		planeEqn: ThreeDimPlaneEquation object, defines the plane
		keepInPlane: Bool, If True we keep co-ords in the same plane; else we keep co-ords NOT in the plane
		planeTolerance: float, Points further from the plane than this value are considered to be in different planes
			 
	Returns
		 outIndices: iter of ints, Indices in inpCoords for the filtered co-ordinates
 
	"""
	outOfPlaneIndices = list()
	inPlaneIndices = list()

	for idx,x in enumerate(inpCoords):
		currDist = planeEqn.getDistanceOfPointFromPlane(x)
		if (currDist < planeTolerance):
			inPlaneIndices.append(idx)
		else:
			outOfPlaneIndices.append(idx)

	if keepInPlane:
		outIndices = inPlaneIndices
	else:
		outIndices = outOfPlaneIndices

	return outIndices


#General function for getting nearest distances/coords to points
def getNearestDistanceToInputPoint(inpPoint, otherPoints):
	""" Return distance to point closest to inpPoint. For Example, inpPoint may be co-ordinates of an atom while otherPoints are all the neighbour co-ordinates; in which case this function will give the nearest neighbour distance (or one of them at least)
	
	Args:
		inpPoint: iter (usually len-3), co-ordinates of input point
		otherPoints: iter of iters, co-ordinates of other points
			 
	Returns
		 distance: float, Distance to nearest point
 
	"""
	nearestCoord = getNearestCoordToInputPoint(inpPoint, otherPoints)
	return vectHelp.getDistTwoVectors(inpPoint, nearestCoord)

def getNearestCoordToInputPoint(inpPoint, otherPoints):
	""" Return co-ordinates of point closest to inpPoint. For Example, inpPoint may be co-ordinates of an atom while otherPoints are all the neighbour co-ordinates; in which case this function will give the nearest neighbour (or one of them at least)
	
	Args:
		inpPoint: iter (usually len-3), co-ordinates of input point
		otherPoints: iter of iters, co-ordinates of other points
			 
	Returns
		 outCoord: iter, Co-ordinates of the nearest point
 
	"""
	outIdx = getIdxOfNearestPointToInputPoint(inpPoint, otherPoints)
	return [x for x in otherPoints[outIdx]]

def getIdxOfNearestPointToInputPoint(inpPoint, otherPoints):
	outIdx = 0
	minDistance = vectHelp.getDistTwoVectors(inpPoint, otherPoints[0])
	for idx,coord in enumerate(otherPoints[1:],1):
		currDist = vectHelp.getDistTwoVectors(inpPoint,coord)
		if (currDist < minDistance):
			minDistance = currDist
			outIdx = idx
	return outIdx

def getClosestDistanceBetweenTwoPoints(inpPoints):
	""" When given a list of input points, return the shortest distance between any two points
	
	Args:
		iter (usually len-3), co-ordinates of input point
			 
	Returns
		 distance: (float) Shortest distance between two points
 
	"""
	#Get a distance matrix with diag terms set to infinity so they dont interfere with np.min
	distMatrix = _getDistMatrixForSetOfCoords(inpPoints)
	for idx,unused in enumerate(distMatrix):
		distMatrix[idx][idx] = np.inf

	return np.min(distMatrix) 



def getClosestDistanceBetweenTwoElementsForInpCell(inpCell, eleA, eleB, inclImages=True, inclImageDims=None):
	""" Gets the closest distance between two elements in a list of cartesian co-ordinates
	
	Args:
		inpCell: (plato_pylib UnitCell object) 
		eleA: (str) Symbol for first element of interest
		eleB: (str) Symbol for second element of interest
		inclImages: (Bool, Default=True) Whether to include images
		inclImageDims: (len-3 bool-iter, Default=[True,True,True]) If inclImages is true then this determines the dimensions for which we include images. Setting some to False can lead to significantly faster runtimes

	Returns
		 outDist: (float) The smallest distance between eleA and eleB
 
	"""
	inclImageDims = [True,True,True] if inclImageDims is None else inclImageDims
	if inclImages is False:
		inclImageDims = [False, False, False]

	#Get distance matrices for central-central and central-image cells
	centralCoords, imageCoords = _getCentralAndImageCoordsFromInpCell(inpCell, imageDims=inclImageDims)
	distMatrixCentral = _getDistMatrixForSetOfCoords(centralCoords)
	distMatrixCentralAndImages = _getDistMatrixBetweenTwoSetsOfSeparateCoords(centralCoords, imageCoords)

	#figure out the relevant indices in central/image cells
	indicesACentral = [idx for idx,coord in enumerate(centralCoords) if coord[-1].upper()==eleA.upper()]
	indicesBCentral = [idx for idx,coord in enumerate(centralCoords) if coord[-1].upper()==eleB.upper()]

	if any(inclImageDims):
		indicesAImages  = [idx for idx,coord in enumerate(imageCoords)   if coord[-1].upper()==eleA.upper()]
		indicesBImages  = [idx for idx,coord in enumerate(imageCoords)   if coord[-1].upper()==eleB.upper()]
	else:
		indicesAImages = list()
		indicesBImages = list()

	#TODO: Refactor these two loops into functions (taking indices and distance matrices as args)
	#Get the minimum distance in the central cell
	minDistCentral = np.inf
	for idxA in indicesACentral:
		for idxB in indicesBCentral:
			if idxA!=idxB:
				currDist = distMatrixCentral[idxA][idxB]
				if currDist < minDistCentral:
					minDistCentral = currDist

	#TODO: Sort out the image cells next
	minDistWithImages = np.inf
	for idxA in indicesACentral:
		for idxB in indicesBImages:
			currDist = distMatrixCentralAndImages[idxA][idxB]
			if currDist <minDistWithImages:
				minDistWithImages = currDist

	minDist = min( [minDistCentral,minDistWithImages] )

	return minDist



def _getCentralAndImageCoordsFromInpCell(inpCell, imageDims=None):
	imageDims = [True, True, True] if imageDims is None else imageDims

	kwargs = {"alongA":imageDims[0], "alongB":imageDims[1], "alongC":imageDims[2], "removeStartCoords":True}
	cellWithImages = supCellHelp.getUnitCellSurroundedByNeighbourCells(inpCell, **kwargs)
	centralCoords = copy.deepcopy(inpCell.cartCoords)
	imageCoords = copy.deepcopy(cellWithImages.cartCoords)
	return centralCoords, imageCoords

#Works using distances
#NOTE: NOW HAS NO USE (since i changed implementation for getting central and image coords)
def _getIndicesOfDuplicatedAtomsInCoordsB(coordsA, coordsB, distTol=1e-2):
	duplicatedIndices = list()
	for idxA,coordA in enumerate(coordsA):
		for idxB,coordB in enumerate(coordsB):
			currDist = vectHelp.getDistTwoVectors(coordA[:3], coordB[:3])
			if (currDist<distTol):
				duplicatedIndices.append(idxB)
	return duplicatedIndices

def _getDistMatrixForSetOfCoords(inpCoords):
	return _getDistMatrixBetweenTwoSetsOfSeparateCoords(inpCoords, inpCoords)

def _getDistMatrixBetweenTwoSetsOfSeparateCoords(coordsA, coordsB):
	outMatrix = np.zeros( (len(coordsA), len(coordsB)) )
	for rIdx,coordA in enumerate(coordsA):
		for cIdx,coordB in enumerate(coordsB):
			outMatrix[rIdx][cIdx] = vectHelp.getDistTwoVectors(coordA[:3], coordB[:3])
	return outMatrix




