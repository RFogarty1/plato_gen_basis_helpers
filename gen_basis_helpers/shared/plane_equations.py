
import itertools as it
import math

import numpy as np

class ThreeDimPlaneEquation():
	""" Class representing a 3-dimension plane equation ax + by +cz = d

	"""
	def __init__(self, a, b, c, d):
		""" Initializer
		
		Args:
			a,b,c,d (floats): Parameters for the plane equation ax + by + cz = d
				 
		"""
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


	def getDistanceOfPointFromPlane(self, inpXyz):
		""" Calculates the distance of a point from the plane
		
		Args:
			inpXyz: (len 3 float iter) [x,y,z] co-ordinates
				 
		Returns
			outDist: (float) The distance between input point and the nearest point on this plane
	 
		"""
		#Step 1 = Find the d value for the parralel plane this lies on; the displacement vector between point and plane is then the vector normal to this plane
		#Step 2  = use the fact the normal vector points the same way as the displacement vecctor to get a length
		dValForThisPoint = self.calcDForInpXyz(inpXyz)
		diffInDVals = abs( dValForThisPoint - self.d )
		lenNormalVectorToThisPlane = math.sqrt( (self.a**2) + (self.b**2) + (self.c**2) )
		outDist = diffInDVals / lenNormalVectorToThisPlane
		return outDist


	def calcDForInpXyz(self, inpXyz):
		""" For a given xyz calculate d. If d=self.d then the point lies on this plane; else it lies on a parralel plane with d being the output to this function
		
		Args:
			inpXyz: (len 3 float iter) co-ordinates for a point
				 
		"""
		assert len(inpXyz) == 3
		return sum( [param*coord for param,coord in it.zip_longest([self.a,self.b,self.c],inpXyz)] )


	@property
	def coeffs(self):
		""" Returns a,b,c,d coeffs (in ax + by +cz = d) as a len-4 iter
		"""
		return [self.a,self.b,self.c,self.d]

