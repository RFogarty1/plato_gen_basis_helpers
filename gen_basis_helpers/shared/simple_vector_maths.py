
import itertools as it
import math


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

