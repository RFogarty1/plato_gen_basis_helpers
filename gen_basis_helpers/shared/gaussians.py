
import itertools as it
import math

""" Defining the standard gaussian primitive and composite classes """

#Added a base class mainly for the doc-strings
class GauPrimBase():

	def evalFunctAtDists(self, distances):
		""" Evaluate the primitive function at a set of distances from its origin
		
		Args:
			distances: (float iter) A set of distances from the origin of this gaussian
				 
		Returns
			outVals: (float iter) The value of the Gaussian primitive at the input set of distances
	 
		"""
		raise NotImplementedError("")



class GauPrimComposite(GauPrimBase):
	"""Composite of Gaussian primitve functions; meant to represent the linear combination (sum) of the components

	"""

	def __init__(self, objs):
		""" Initializer
		
		Args:
			objs: (iter of GauPrim objects)
				 
		"""
		self.objs = list(objs)

	@property
	def a(self):
		raise NotImplementedError("")

	@a.setter
	def a(self,val):
		raise NotImplementedError("")

	@property
	def c(self):
		raise NotImplementedError("")

	@c.setter
	def c(self,val):
		raise NotImplementedError("")

	@property
	def leaves(self):
		outLeaves = list()
		for x in self.objs:
			outLeaves.extend(x.leaves)
		return outLeaves

	def evalFunctAtDists(self, distances):
		return self._sumOverLeafResults("evalFunctAtDists",distances)

	def getIntegralAllSpace(self):
		return self._sumOverLeafResults("getIntegralAllSpace")

	#For when we need to sum a series of arrays, one per leaf object. try:except could detect is iter or not 
	def _sumOverLeafResults(self, functionName, *args):
		sumVals = getattr(self.objs[0], functionName)(*args)
		for x in self.objs[1:]:
			currVals = getattr(x,functionName)(*args) #May be single-valued or an iterable
			try:
				sumVals = [x+y for x,y in it.zip_longest(sumVals,currVals)]
			except TypeError:
				sumVals += currVals
		return sumVals

	def __eq__(self, other):

		if len(self.leaves) != len(other.leaves):
			return False

		for leafA, leafB in zip(self.leaves,other.leaves):
			if leafA != leafB:
				return False

		return True


class GauPrim(GauPrimBase):
	def __init__(self,a,c,pos):
		self._eqTol = 1e-6
		self.a = a
		self.c = c
		self.pos = pos

	@classmethod
	def fromZDirOnly(cls, a:float, c:float, dist:float):
		return cls(a,c,[0,0,dist])

	@classmethod
	def fromExpAndCoeffOnly(cls, exponent, coeff):
		return cls(exponent, coeff, [0,0,0])

	@property
	def leaves(self):
		""" Returns all leaves (nodes) of composite object. Also implemented on non-composite to allow composite/non-composite to be compatable """
		return [self]

	def evalFunctAtPos(self,pos):
		rVal = _getDistanceTwoPos(self.pos,pos)
		return self.c*math.exp(-1*rVal*rVal*self.a)

	def evalFunctAtDists(self, distances):
		outList = [None for x in distances]
		for idx,x in enumerate(distances):
			outList[idx] = self.c*math.exp(-1*x*x*self.a)
		return outList

	def getIntegralAllSpace(self):
		return self.c*(math.sqrt( math.pi/self.a )**3)


	def __eq__(self, other):

		if len(other.leaves) != 1:
			return False #Must be a composite with multiple primitives; which menas it CANT be equal

		floatAttrs = ["a","c"]
		listAttrs = ["pos"]
		otherLeaf = other.leaves[0] 

		eqTol = min(self._eqTol, otherLeaf._eqTol)

		for currAttr in floatAttrs:
			if abs( getattr(self,currAttr) - getattr(otherLeaf,currAttr) ) > eqTol:
				return False

		for currAttr in listAttrs:
			selfList, otherList = getattr(self,currAttr), getattr(otherLeaf,currAttr)
			if len(selfList) != len(otherList):
				return False
			for x,y in zip(selfList,otherList):
				if abs(x-y) > eqTol:
					return False

		return True

def _getDistanceTwoPos(posA,posB):
	diff = [x-y for x,y in it.zip_longest(posA,posB)]
	sqrDiff = sum([x**2 for x in diff])
	return math.sqrt(sqrDiff)

