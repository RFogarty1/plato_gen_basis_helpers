import copy
import itertools as it
import math
import os

import plato_pylib.parseOther.parse_cp2k_basis as parseCP2KBasis

from .. import parse_basis as parseBasHelp
from ...shared import gaussians as gaussians
from ...shared import basis_set_objs as basObjs

class BasisParser(parseBasHelp.OrbitalBasisSetParserTemplate):

	registeredKwargs = set(parseBasHelp.OrbitalBasisSetParserTemplate.registeredKwargs)
	registeredKwargs.add("basisSetFolder")
	registeredKwargs.add("basisSetFile")
	registeredKwargs.add("element")
	registeredKwargs.add("basisNameInFile") #CP2K specific
	registeredKwargs.add("basisAlias")
	registeredKwargs.add("coeffsForUnormalisedGaussians")


	#Note if you pass None it will overwrite
	def _setDefaultInitAttrs(self):
		self.coeffsForUnormalisedGaussians=True

	#We simply need the basisName and the basis functions in an iter to be passed to the constructor. Triple mock should be fine
	def _parseFromSelf(self):
		nativeFmtFullBasis = self._getBasisInNativeFormat()
		allBasisFuncts = getBasisFunctObjsFromCP2KBasisSetNativeFormat(nativeFmtFullBasis, convNormToRawGaussians=self.coeffsForUnormalisedGaussians)
		return basObjs.OrbitalBasisSetStandard(allBasisFuncts, self.basisAlias)

	def _getBasisInNativeFormat(self):
		fullParsedFile = parseCP2KBasis.parseCP2KBasisFile(self._basisPath)
		specificBasis = fullParsedFile.getUniqueBasisSet(self.element, self.basisNameInFile)
		return specificBasis

	@property
	def _basisPath(self):
		return os.path.join(self.basisSetFolder, self.basisSetFile)

	@property
	def basisAlias(self):
		""" The name you want to give to this basis set (may be used, for example, in plotting the radial function) """
		if self._basisAlias is not None:
			return self._basisAlias
		else:
			return self.basisNameInFile

	@basisAlias.setter
	def basisAlias(self,val):
		self._basisAlias = val




def getBasisFunctObjsFromCP2KBasisSetNativeFormat(cp2kBasisSet, convNormToRawGaussians=True):
	""" Extracts basis functions from the CP2K basis set format returned by the cp2k basis parser
	
	Args:
		cp2kBasisSet: (BasisSetCP2K object) 
	 
	Returns
		 basisFuncts: (iter of GauSumOrbitalBasisFunction objects) Each represents a single orbital using a linear combination of gaussian functions
 
	"""
	#Get the functions out
	outBasisFuncts = list()
	for x in cp2kBasisSet.exponentSets:
		outBasisFuncts.extend( _getBasisFunctObjFromExponentSet(x, convNormToRawGaussians=convNormToRawGaussians)  )
	
	#Index the functions with varying l-values
	lCounts = { l:0 for l in set([x.lVal for x in outBasisFuncts]) }

	for x in outBasisFuncts:
		x.nVal = lCounts[x.lVal] + 1
		lCounts[x.lVal] = lCounts[x.lVal] + 1

	return outBasisFuncts


def _getBasisFunctObjFromExponentSet(expSet, convNormToRawGaussians=True):
	""" Converts the exponent-set format usd to store groups of CP2K basis functions into a list of single basis function objects
	
	Args:
		expSet: (plato_pylib ExponentSetCP2K object) Represents a group of basis function which share the same exponents
		convNormToRawGaussians: (Bool) If true the the coefficients can be directly applied to the gaussian functions (the plato way); if false the coefficients need to be multiplied by the NORMALISED gaussian for each exponent (the CP2K way). Note that the normalisations include the angular dependence (i.e. its the full gaussian thats normalised; not just the radial part)
 
	Returns
		 basisFuncts: (iter of GauSumOrbitalBasisFunction objects) Each represents a single orbital using a linear combination of gaussian functions
 
	NOTE: This function should NOT be used directly. It doesnt properly label the n-values of basis functions with varying l-values (sets all to 1) since a hgiher level function is always expected to sort that anyway

	"""
	allExponents = [x for x in expSet.exponents]
	allNVals = 1 #Expect a higher level function to clean this up

	allBasisFuncts = list()
	for coeffs, lVal in it.zip_longest(expSet.coeffs,expSet.lVals):
		if convNormToRawGaussians:
			normCoeffs = [coeff*calcNormConstantForCP2KOnePrimitive(exponent,lVal) for exponent,coeff in it.zip_longest(allExponents,coeffs)]
		else:
			normCoeffs = coeffs
		currPrims = [gaussians.GauPrim.fromExpAndCoeffOnly(a,c) for a,c in it.zip_longest(allExponents,normCoeffs)]
		gauFunct = gaussians.GauPrimComposite(currPrims)
		currBasisFunct = basObjs.GauSumOrbitalBasisFunction(allNVals, lVal, gauFunct)
		allBasisFuncts.append( currBasisFunct )

	return allBasisFuncts


#Taken from plato script to convert plato basis sets into cp2k
def calcNormConstantForCP2KOnePrimitive(exponent:float, angMom:int):
	expZet = 0.25 * ((2*angMom)+3)
	preFac = (2**angMom) * ((2/math.pi)**0.75) 
	normFactor = preFac*(exponent**expZet)
	return normFactor




