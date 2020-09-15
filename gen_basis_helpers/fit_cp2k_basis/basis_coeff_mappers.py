import math
import plato_pylib.plato.parse_gau_files as parseGau

from ..gau_prod_theorem import get_ints_s_expansions as sIntHelp


from . import core

class FitCoeffsToBasisFunctionExponentsAndCoeffsBase():
	""" Has interface __call__(fitCoeffs)->(exponents,coeffs) where exponents and coeffs are of the same length and lead to the full basis function
	"""

	def __call__(self, coeffs):
		raise NotImplementedError("")


#TODO: Ideally swap the coeffs/exponents order in the initializer constructor to match others
class FitCoeffsToBasisFunctionExponentsAndCoeffsMixedOptStandard(FitCoeffsToBasisFunctionExponentsAndCoeffsBase):
	""" Has interface __call__(fitCoeffs)->(exponents,coeffs) where exponents and coeffs are of the same length and lead to the full basis function. NOTE: This implementation requires that exponents MUST be listed first in fitCoeffs. Furthermore coeffs/exponents need to be ordered a certain way (see initializer)
	"""

	def __init__(self, fixedCoeffs, fixedExponents):
		""" Initializer
		
		Args:
			IMPORTANT NOTE: The fixedCoeffs/fixedExponents must both be the first n-values in the basis function (if exponents are a1,a2,a3 you cant fix a1 and a3 with this implementation, or even just a3 on its own)
			fixedCoeffs: (iter of floats) The values of the fixed coefficients. Pass a blank list to vary ALL coefficients
			fixedExponents: (iter of floats) The values of the fixed exponents. Pass a blank list to vary ALL exponents
				 
		"""
		self.fixedCoeffs = list(fixedCoeffs)
		self.fixedExponents = list(fixedExponents)

	def __call__(self, fitCoeffs):
		nTotal = len(fitCoeffs) + len(self.fixedCoeffs) + len(self.fixedExponents)
		assert nTotal%2==0, "Number of coefficients needs to match number of exponents; not possible with {} total variables".format(nTotal)
		nEach = int( int(nTotal)/int(2) )

		nFreeExponents = nEach - len(self.fixedExponents)
		nFreeCoeffs = nEach - len(self.fixedCoeffs)

		outExponents = list(self.fixedExponents) + list(fitCoeffs[:nFreeExponents])
		outCoeffs = list(self.fixedCoeffs) + list(fitCoeffs[nFreeExponents:])

		return outExponents, outCoeffs

class CoeffsToFullBasisFunctionForMixedCoeffExponentOpt(core.CoeffsTransformer):

	def __init__(self, fixedExponents, fixedCoeffs, angMom):
		self.fixedExponents = fixedExponents
		self.fixedCoeffs = fixedCoeffs
		self.angMom = angMom

	def __call__(self, coeffs):
		mapper = FitCoeffsToBasisFunctionExponentsAndCoeffsMixedOptStandard(self.fixedCoeffs, self.fixedExponents)
		exponents, outCoeffs = mapper(coeffs)
		outFunct = parseGau.GauPolyBasis(exponents, [outCoeffs], label=self.angMom)
		return outFunct


class CoeffsToFullBasisSetForMixedCoeffExponentOpt():
	
	def __init__(self, fixedBasisFuncts, fixedExponents, fixedCoeffs, angMom):
		""" Initializer
		
		Args:
			IMPORTANT NOTE: The fixedCoeffs/fixedExponents must both be the first n-values in the basis function (if exponents are a1,a2,a3 you cant fix a1 and a3 with this implementation, or even just a3 on its own)

			fixedBasisFuncts: (iter of plato_pylib GauPolyBasis objects) Each represents 1 fixed basis function. Note that label should be set to the angular momentum value of each orbital
			fixedCoeffs: (iter of floats) The values of the fixed coefficients. Pass a blank list to vary ALL coefficients
			fixedExponents: (iter of floats) The values of the fixed exponents. Pass a blank list to vary ALL exponents
			angMom: (int) The angular momentum of the new basis function
				 
		"""
		self.fixedBasisFuncts = list(fixedBasisFuncts)
		self.fixedExponents = list(fixedExponents)
		self.fixedCoeffs = list(fixedCoeffs)
		self.angMom = angMom
		
	def __call__(self, coeffs):
		currArgs = [self.fixedExponents, self.fixedCoeffs, self.angMom]
		newFunctMapper = CoeffsToFullBasisFunctionForMixedCoeffExponentOpt(*currArgs)
		newFunct = newFunctMapper(coeffs)
		return self.fixedBasisFuncts + [newFunct]


class CoeffsToFullBasisFunctionForFixedExponents(core.CoeffsTransformer):
	"""Class holds the exponents/angular momentum for a basis function and returns the whole basis function when given the curent coefficients for the the exponents

	"""
	def __init__(self, exponents, angMom):
		""" Initializer
		
		Args:
			exponents: (iter of floats) Each represents one exponent
			angMom: (int) Angular momentum for the orbital
				 
		"""
		self.exponents = list(exponents)
		self.angMom = angMom

	def __call__(self, coeffs):
		return parseGau.GauPolyBasis(self.exponents, [coeffs], label=self.angMom)


class CoeffsToFullBasisSetForFixedExponents(core.CoeffsTransformer):
	""" Class holds a full basis set + the angular momentum/exponents of a final basis function being optimised. This allows it to return the FULL basis set when given the current coefficients for the last basis function (the one being optimised)

	"""
	def __init__(self, fixedOrbExpansion, exponents, angMom):
		""" Initializer
		
		Args:
			fixedOrbExpansion: (iter of plato_pylib GauPolyBasis objects) Each represents 1 fixed basis function. Note that label should be set to the angular momentum value of each orbital
			exponents: (iter of floats) Each represents one exponent for the basis function we're optimising (which is NOT present in fixedOrbExpansAngular momentum for the orbitalion)
			angMom: (int) Angular momentum for the orbital we're optimising
		"""
		self.fixedOrbExpansion = list(fixedOrbExpansion)
		self.exponents = list(exponents)
		self.angMom = angMom

	def __call__(self, coeffs):
		newBasisFunct = parseGau.GauPolyBasis(self.exponents, [coeffs], label=self.angMom)
		outObj = self.fixedOrbExpansion + [newBasisFunct]
		return outObj



class CoeffsToNormalisedValuesForMixedCoeffExponentOpt(core.CoeffsTransformer):
	""" Converts basis set coefficients to values leading to a normalised basis function ( <\phi|\phi>=1 )

	"""

	def __init__(self, fixedExponents, fixedCoeffs, angMom):
		""" Initializer
		
		Args:
			IMPORTANT NOTE: The fixedCoeffs/fixedExponents must both be the first n-values in the basis function (if exponents are a1,a2,a3 you cant fix a1 and a3 with this implementation, or even just a3 on its own)
			fixedCoeffs: (iter of floats) The values of the fixed coefficients. Pass a blank list to vary ALL coefficients
			fixedExponents: (iter of floats) The values of the fixed exponents. Pass a blank list to vary ALL exponents
			angMom: (int) The angular momentum of the orbital
				 
		"""
		self.fixedCoeffs = fixedCoeffs
		self.fixedExponents = fixedExponents
		self.angMom = angMom

	def __call__(self, coeffs):
		raise NotImplementedError("")
#		mapToExponentsAndCoeffs = FitCoeffsToBasisFunctionExponentsAndCoeffsMixedOptStandard(self.fixedCoeffs,self.fixedExponents)
#		exponents, gauCoeffs = mapToExponentsAndCoeffs(coeffs)
#		fixedExponentObj = CoeffsToNormalisedValuesFixedExponents(exponents,self.angMom)
#		normalisedCoeffs = fixedExponentObj(gauCoeffs)
#		return fixedExponentObj(gauCoeffs)

class CoeffsToNormalisedValuesFixedExponents(core.CoeffsTransformer):
	""" Converts basis set coefficients to values leading to a normalised basis function ( <\phi|\phi>=1 )

	"""

	def __init__(self, exponents, angMom):
		""" Initializer
		
		Args:
			exponents: (iter of floats) Exponents for the orbital
			angMom: (int) Angular momentum of the orbital

		"""
		self.exponents = list(exponents)
		self.angMom = angMom

	def _getGauPolyBasisFromCoeffs(self, coeffs):
		mapFunct = CoeffsToFullBasisFunctionForFixedExponents(self.exponents, self.angMom)
		return mapFunct(coeffs)

	def _getScaleFactor(self, coeffs):
		gauPolyBasisObj = self._getGauPolyBasisFromCoeffs(coeffs)
		distVal = 0.0

		if (self.angMom == 0):
			overlapVal = sIntHelp.getSelfOverlapMcWedaWeightFromGauPolyBasis(gauPolyBasisObj, distVal)
		else:
			raise ValueError("{} is an unsupported value for self.angMom".format(self.angMom))

		return 1/ math.sqrt(overlapVal)

	def __call__(self, coeffs):
		scaleFactor = self._getScaleFactor(coeffs)
		return [x*scaleFactor for x in coeffs]

