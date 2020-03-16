#!/usr/bin/python3

import itertools as it
import math

def main():
	primA = GauPrim(0.3,1.2,[0.0,0.0,0.0])
#	primB = GauPrim(0.5,1.5,[0.0,0.0,1.2])
	primB = GauPrim.fromZDirOnly(0.5,1.5,1.2)
	primC = combineTwoGauPrims(primA,primB)

	testPos = [1.0,5.4,2.3]

	print("primA(r)*primB(r) = {}".format( primA.evalFunctAtPos(testPos)*primB.evalFunctAtPos(testPos) ) )
	print("primC(r) = {}".format(primC.evalFunctAtPos(testPos)) )


def combineTwoGauPrims(gPrimA,gPrimB):
	gamma = gPrimA.a + gPrimB.a

	#First get new position vector
	scaledVectA = [x*gPrimA.a for x in gPrimA.pos]
	scaledVectB = [x*gPrimB.a for x in gPrimB.pos]
	newPosVect = [(x+y)/gamma for x,y in zip(scaledVectA,scaledVectB)]

	#Now get the prefactor
	difPosVect = [a-b for a,b in zip(gPrimA.pos,gPrimB.pos)]
	sqrDist = sum([x**2 for x in difPosVect])
	prefactor = math.exp( (-1*gPrimA.a*gPrimB.a*sqrDist)/gamma )
	coeffProd = gPrimA.c*gPrimB.c


	return GauPrim(gamma, coeffProd*prefactor,newPosVect)

#Now duplicated in the shared/. The version in shared should be used; this version will be made to simply
#inherit from it later (not removing since its needed for backwards compat.)
class GauPrim():
	def __init__(self,a,c,pos):
		self.a = a
		self.c = c
		self.pos = pos

	@classmethod
	def fromZDirOnly(cls, a:float, c:float, dist:float):
		return cls(a,c,[0,0,dist])

	def evalFunctAtPos(self,pos):
		rVal = getDistanceTwoPos(self.pos,pos)
		return self.c*math.exp(-1*rVal*rVal*self.a)

	def getIntegralAllSpace(self):
		return self.c*(math.sqrt( math.pi/self.a )**3)

def getDistanceTwoPos(posA,posB):
	diff = [x-y for x,y in it.zip_longest(posA,posB)]
	sqrDiff = sum([x**2 for x in diff])
	return math.sqrt(sqrDiff)

if __name__ == '__main__':
	main()
