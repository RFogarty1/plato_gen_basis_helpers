
import plato_pylib.plato.parse_gau_files as parseGau

from ..gau_prod_theorem import get_ints_s_expansions as sIntHelp


from . import core


class CoeffsToFullBasisFunctionForFixedExponents(core.CoeffsTransformer):
	"""Class holds the exponents/angular momentum for a basis function and returns the whole basis function when given the curent coefficients for the the exponents

	"""
	def __init__(self, exponents, angMom):
		""" Initializer
		
		Args:
			exponents: (iter of floats) Each represents one exponent
			angMom: (int) Angular momentum for the orbital
				 
		"""
		self.exponents = exponents
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
		self.exponents = exponents
		self.angMom = angMom

	def __call__(self, coeffs):
		newBasisFunct = parseGau.GauPolyBasis(self.exponents, [coeffs], label=self.angMom)
		outObj = self.fixedOrbExpansion + [newBasisFunct]
		return outObj


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

		return 1/overlapVal

	def __call__(self, coeffs):
		scaleFactor = self._getScaleFactor(coeffs)
		return [x*scaleFactor for x in coeffs]

