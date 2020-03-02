
import itertools as it
from . import base_flow as baseFlow

import plato_pylib.utils.elastic_consts as elasticHelp


class CrystalStrain():
	"""Representation of a crystal strain (with unit strain parameter).

	"""


	def __init__(self, strainVals):
		""" Initializer
		
		Args:
			strainVals: (length 6 iter) Each value represents the coefficient for an individual strain matrix in our chosen basis set (which is a pretty standard one). Setting each to 1 and looking at the strain matrix is probably the simplest way to see the basis used.
				
		"""
		self._eqTol = 1e-5
		self.strainVals = list(strainVals)

	@property
	def strainMatrix(self):
		""" 3x3 np array representing the crystal strain [[xx,xy,xz],[yx,yy,yz],[zx,zy,zz]]
		
		"""
		outMatrix = _getUnitStrainMatrix(1)*self.strainVals[0]
		for idx,x in enumerate(self.strainVals[1:],2):
			outMatrix += x*_getUnitStrainMatrix(idx)
		return outMatrix

	def __eq__(self,other):
		if isinstance(other, CrystalStrain): 
			eqTol = max(self._eqTol, other._eqTol) #We want to use the loosest equality definition.
			diffs = [x-y for x,y in it.zip_longest(self.strainVals, other.strainVals)]
			if all([x<eqTol for x in diffs]):
				return True
			else:
				return False
		return NotImplemented #Delegates the equality to other

	def toStr(self):
		""" Str representation of strain in terms of a linear sum of coefficients and basis functions
		"""
		minStrainCoeff = 1e-4
		nonZeroIndices = list()
		for idx,x in enumerate(self.strainVals):
			if abs(x) > minStrainCoeff:
				nonZeroIndices.append(idx)

		outCoeffs = [x for idx,x in enumerate(self.strainVals) if idx in nonZeroIndices]
		outIndices = [idx for idx,x in enumerate(self.strainVals,1) if idx-1 in nonZeroIndices]
		return "+".join(["{}eps{}".format(coeff,idx) for idx,coeff in it.zip_longest(outIndices,outCoeffs)])



def _getUnitStrainMatrix(matrixNumb:int):
	strainParam = 1.0
	return elasticHelp._STRAIN_MATRIX_DICT[matrixNumb](1.0)



