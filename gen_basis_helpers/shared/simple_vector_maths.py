
import itertools as it
import math

import numpy as np

def getUnitVectorFromInpVector(inpVector):
	lenVect = getLenOneVector(inpVector)
	return [x/lenVect for x in inpVector]

def getLenOneVector(vectA):
	return math.sqrt( sum([x**2 for x in vectA]) )

def getDistTwoVectors(vectA,vectB):
	sqrDiff = [ (a-b)**2 for a,b in it.zip_longest(vectA,vectB) ]
	return math.sqrt( sum(sqrDiff) )

def getAngleTwoVectors(vectA,vectB):
	normFactorA = math.sqrt( sum( [x**2 for x in vectA] ) )
	normFactorB = math.sqrt( sum( [x**2 for x in vectB] ) )

	normA = [x/normFactorA for x in vectA]
	normB = [x/normFactorB for x in vectB]

	dotProd = getDotProductTwoVectors(normA,normB)
	return math.degrees( math.acos(dotProd) )

def getDotProductTwoVectors(vectA,vectB):
	return sum( [a*b for a,b in it.zip_longest(vectA,vectB)] )


#https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
def getRotationMatrixLinkingTwoUnitVectors(uVectA, uVectB):
	assert abs(1 - getLenOneVector(uVectA))<1e-3
	assert abs(1 - getLenOneVector(uVectB))<1e-3
	assert len(uVectA)==3
	assert len(uVectB)==3

	crossProd = np.cross(uVectA, uVectB)
	v1,v2,v3 = crossProd
	vx = np.array( [ [0    , -1*v3,  1*v2],
	                 [v3   ,  0   , -1*v1],
	                 [-1*v2, v1   , 0    ] ] )

	vxSquared = vx.dot(vx)
	cosTheta = np.dot(uVectA,uVectB)
	angularFactor = 1 / (1 + cosTheta)
	vxSquaredWithAngular = angularFactor*vxSquared
	rotMatrix = np.identity(3) + vx + vxSquaredWithAngular

	return rotMatrix


