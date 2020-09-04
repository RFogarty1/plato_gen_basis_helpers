
import plato_pylib.parseOther.parse_cp2k_basis as parseCP2KBasis

from . import core

class BasisCoeffUpdaterCP2KStandard(core.CoeffObserver):
	""" Class is responsible for converting changes in coefficients (for the basis functions being optimised) into changes in the relevant basis file in CP2K


	"""
	def __init__(self, filePath, eleName, basisNames, coeffToBasisMapper ):
		self.filePath = filePath
		self.eleName = eleName
		self.basisNames = basisNames
		self.coeffToBasisMapper = coeffToBasisMapper
		self.coeffs = None

	def updateCoeffs(self, coeffs):
		self.coeffs = coeffs
		self._updateUsingCurrentCoeffs()

	def _updateUsingCurrentCoeffs(self):
		outFileObj = self._getCP2KBasisFullFileObj()
		parseCP2KBasis.writeBasisFileFromParsedBasisFileObj(self.filePath, outFileObj)

	def _getBasisSetCP2KFormat(self):
		basisStdFmt = self.coeffToBasisMapper(self.coeffs)
		angMomVals = [x.label for x in basisStdFmt]
		getCP2KBasisFunct = parseCP2KBasis.getCP2KBasisFromPlatoOrbitalGauPolyBasisExpansion #Just to make lines shorter
		cp2kBasis = getCP2KBasisFunct(basisStdFmt, angMomVals, self.eleName, basisNames=self.basisNames)
		return cp2kBasis

	#TODO: Ideally it would be nice to try parsing the initial self.filePath and incorporating the new basis into it; instead of OVERWRITING the whole file
	def _getCP2KBasisFullFileObj(self):
		ourBasis = self._getBasisSetCP2KFormat()
		outObj = parseCP2KBasis.ParsedBasisFileCP2K(self.filePath,[ourBasis])
		return outObj
