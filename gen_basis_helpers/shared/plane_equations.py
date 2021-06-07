
import itertools as it
import math

import numpy as np

from . import simple_vector_maths as vectHelp

class ThreeDimPlaneEquation():
	""" Class representing a 3-dimension plane equation ax + by +cz = d

	"""
	def __init__(self, a, b, c, d):
		""" Initializer
		
		Args:
			a,b,c,d (floats): Parameters for the plane equation ax + by + cz = d
				 
		"""
		self._eqTol = 1e-6
		self.a = a
		self.b = b
		self.c = c
		self.d = d

	@classmethod
	def fromTwoPositionVectors(cls, inpVectA, inpVectB, normaliseCoeffs=True):
		""" Alternative initializer. Creates the object using two input POSITION vectors (i.e. both MUST pass through origin)
		
		Args:
			inpVectA: (len 3 float iter) Position vector
			inpVectB: (len 3 float iter) Position vector
			normaliseCoeffs: (Bool, Optional) If True then always return coefficients for the normalised (i.e. unit) vector normal to the plane. Default=True
				 
		"""
		vectA, vectB = np.array(inpVectA), np.array(inpVectB)
		normVect = np.cross(vectA,vectB)
		if normaliseCoeffs:
			lenNormVect = math.sqrt( sum([x**2 for x in normVect]) )
			normVect = np.array( [x/lenNormVect for x in normVect] )

		#D is always zero, since the position vectors are also points on the plane
		outCoeffs = [x for x in normVect] + [0] 

		return cls(*outCoeffs)

	def __eq__(self,other):
		
		eqTol = min(self._eqTol, other._eqTol)

		for cA,cB in it.zip_longest(self.coeffs, other.coeffs):
			if (abs(cA-cB) > eqTol):
				return False

		return True


	def getSignedDistanceOfPointFromPlane(self, inpXyz):
		""" Calculates the signed fistance of a point from the plane. If the normal vector points towards the point, the distance is +ve, if it points in the opposite direction it is negative 

		Args:
			inpXyz: (len 3 float iter) [x,y,z] co-ordinates

		Returns:
			outDist: (float) The signed distance between the input point and the nearest point on this plane
		"""
		#Step 1 = Find the d value for the parralel plane this lies on; the displacement vector between point and plane is then the vector normal to this plane
		#Step 2  = use the fact the normal vector points the same way (or the exact opposite way) as the displacement vector to get a signed distance
		dValForThisPoint = self.calcDForInpXyz(inpXyz)
		diffInDVals = dValForThisPoint - self.d #Using the absolute value here would give the unsigned distance
		lenNormalVectorToThisPlane = math.sqrt( (self.a**2) + (self.b**2) + (self.c**2) )
		outDist = diffInDVals / lenNormalVectorToThisPlane
		return outDist


	def getDistanceOfPointFromPlane(self, inpXyz):
		""" Calculates the distance of a point from the plane
		
		Args:
			inpXyz: (len 3 float iter) [x,y,z] co-ordinates
				 
		Returns
			outDist: (float) The distance between input point and the nearest point on this plane
	 
		"""
		return abs( self.getSignedDistanceOfPointFromPlane(inpXyz) )


	def calcDForInpXyz(self, inpXyz):
		""" For a given xyz calculate d. If d=self.d then the point lies on this plane; else it lies on a parralel plane with d being the output to this function
		
		Args:
			inpXyz: (len 3 float iter) co-ordinates for a point
				 
		"""
		assert len(inpXyz) == 3
		return sum( [param*coord for param,coord in it.zip_longest([self.a,self.b,self.c],inpXyz)] )


	
	def getPointClosestToOrigin(self):
		""" Returns the point on this plane closest to the origin """
		coeffs = self.coeffs
		normVal = 1/( sum([x**2 for x in self.coeffs[:3]]) )
		outPoint = [x*coeffs[-1]*normVal for x in self.coeffs[:3]]

		#Internal checks; can probably remove later
		#Firstly check point is expected distance from plane
		errorTol = 1e-4 #
		expDist = self.getDistanceOfPointFromPlane([0,0,0])
		actDist = vectHelp.getLenOneVector(outPoint)
		if abs(expDist-actDist)>errorTol:
			raise ValueError("Some mistake in my coding here")
		#Secondly check point is distance of zero from the plane (i.e. it lies on the plane)
		if abs( self.getDistanceOfPointFromPlane(outPoint) ) > errorTol:
			raise ValueError("Some mistake in my coding here")

		return outPoint

	@property
	def coeffs(self):
		""" Returns a,b,c,d coeffs (in ax + by +cz = d) as a len-4 iter
		"""
		return [self.a,self.b,self.c,self.d]

	@coeffs.setter
	def coeffs(self,val):
		self.a, self.b, self.c, self.d = val



def getOutOfPlaneDistTwoPoints(posA, posB, planeEqn):
	""" Description of function
	
	Args:
		posA: (len-3 iter) [x,y,z]
		posB: (len-3 iter) [x,y,z]
		planeEqn: (ThreeDimPlaneEquation)
			 
	Returns
		 outDist: The out of plane distance between the two points. This is the distance that would remain if we shifted posA along the surface normal such that it was in the same plane (with planeEqn) as posB 
 
	IMPORTANT:
		This is obviously NOT aware of periodic boundaries. Likely you want to pass the nearest images in

	"""
	distFromPlaneA = planeEqn.getSignedDistanceOfPointFromPlane(posA[:3])
	distFromPlaneB = planeEqn.getSignedDistanceOfPointFromPlane(posB[:3])
	interPlaneDist = abs( distFromPlaneA-distFromPlaneB)
	totalDist = vectHelp.getDistTwoVectors(posA[:3], posB[:3])
	outOfPlaneDist = math.sqrt( (totalDist**2) - (interPlaneDist**2) )
	return outOfPlaneDist


def getInterPlaneDistTwoPoints(posA, posB, planeEqn):
	""" Gets the inter-plane distance between two points. This is the distance between planeA and planeB, which are both parralel to planeEqn and contain posA and posB respectively
	
	Args:
		posA: (len-3 iter) [x,y,z] Position of point A
		posB: (len-3 iter) [x,y,z] Position of point B
		planeEqn: (ThreeDimPlaneEquation)

	Returns
		outDist: (float) As description says; the distance between parralel planes containing point A and point B (and both parralel to planeEqn) 
 
	"""
	distFromPlaneA = planeEqn.getSignedDistanceOfPointFromPlane(posA)
	distFromPlaneB = planeEqn.getSignedDistanceOfPointFromPlane(posB)
	distAB = abs( distFromPlaneA - distFromPlaneB )
	return distAB

def getVectorToMoveFromParallelPlanesAToB(planeA, planeB, parallelTol=1e-6):
	""" Gets the vector to move (in the shortest distance) from planeA to planeB.
	
	Args:
		planeA: (ThreeDimPlaneEquation)
		planeB: (ThreeDimPlaneEquation)
		parralelTol: (float) The magnitude of dot products between planeA and planeB normal vectors needs to be this close to 1 (tolerance is here to account for float errors)

	Returns
		 outVector: (len-3 iter) Vector allowing movement from planeA to planeB.
 
	NOTE:
		This works by getting the distance closest to origin for one plane, and getting the distance of that point from the other plane. This only works if the two planes are parralel and would give an essentially arbitrary value otherwise. Thus, I'm not sure what kind of errors you'd get for near-parralel planes (e.g. within numerical error). Suspect not a massive issue since near-parralel planes are unlikely (never?) going to intersect near origin

	Raises:
		 ValueError: If the planes intersect. This wont be perfect due to float errors.
	"""
	#Check for non parallel plane (presumably dot products can be used)
	if not checkPlanesAreParralelSimple(planeA, planeB, parallelTol=parallelTol):
		raise ValueError("planeEqns with coefficients {} and {} ARE NOT PARRALEL".format( planeA.coeffs[:3], planeB.coeffs[:3] ))

	#Get the distance between planes if their parralel
	pointOnB = planeB.getPointClosestToOrigin()
	signedDist = planeA.getSignedDistanceOfPointFromPlane(pointOnB)
	normVector = vectHelp.getUnitVectorFromInpVector(planeA.coeffs[:3])
	tVect = [signedDist*x for x in normVector]
	return tVect


def checkPlanesAreParralelSimple(planeEqnA, planeEqnB, parallelTol=1e-6):
	""" Tries to check whether two planes are parallel by using the dot products of their normal vectors. 
	
	Args:
		planeEqnA: (ThreeDimPlaneEquation) 
		planeEqnB: (ThreeDimPlaneEquation)
		parallelTol: The magnitude of dot products between planeA and planeB normal vectors needs to be this close to 1 (tolerance is here to account for float errors)
			 
	Returns
		isParallel: (Bool) Returns True if planes are parallel and False otherwise
 
	"""
	normVectA, normVectB = [ vectHelp.getUnitVectorFromInpVector(planeEqnA.coeffs[:3]), vectHelp.getUnitVectorFromInpVector(planeEqnB.coeffs[:3]) ]
	absDotProd = abs(vectHelp.getDotProductTwoVectors(normVectA, normVectB))
	if abs(absDotProd-1)>parallelTol:
		return False
	return True



