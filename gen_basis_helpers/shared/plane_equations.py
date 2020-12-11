
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

