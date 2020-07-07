
from . import simple_vector_maths as vectHelp


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

