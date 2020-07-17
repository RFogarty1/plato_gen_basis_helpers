
import copy
import plato_pylib.utils.supercell as supCellHelp

from . import plane_equations as planeEqnHelp
from . import simple_vector_maths as vectHelp


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
