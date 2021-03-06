

""" Objects originally meant to make it easy for me to plot the radial distribution of gaussian atomic orbital basis functions """

import math

class OrbitalBasisSetBase():
	"""Class represnting an orbital-based basis set for a single atom (e.g. for when using basis functions similar to atomic orbitals)
	"""

	@property
	def basisName(self):
		""" Str that acts as a label for this basis set """
		raise NotImplementedError("")

	@property
	def orbitalBasisFunctions(self):
		""" iter of OrbitalBasisFunctionBase objects. """
		raise NotImplementedError("")


	def getBasisFunctionsWithLVal(self):
		""" Get all the basis functions for a given angular momentum.
		
		Args:
			lVal: (int) Orbital angular momentum, s=0, l=1, d=2 etc.
				 
		Returns
			outBasisFuncts: (iter of OrbitalBasisFunctionBase objects) All these objects have l=lVal angular momentum AND are ordered by their n-value
 
		"""
		raise NotImplementedError("")


class OrbitalBasisFunctionBase():
	"""Class representing a single atomic-orbital (or similar) basis function
	"""

	@property
	def nVal(self):
		""" Integer used to differentiate basis functions with the same angular momentum. No deeper meaning
		
		"""
		raise NotImplementedError("")

	@property
	def lVal(self):
		""" (int) Angular momentum for this orbital (s=0,p=1,d=2)
		"""
		raise NotImplementedError("")

	def getRadialValsAtDists(self, xVals):
		""" Return the radial values of the orbital basis function at various distances from origin
		
		Args:
			xVals: (float iter) List of distances from origin to get radial values at
				 
		Returns
			yVals: (float iter) Radial value of this basis function at distances xVals from origin
	 
		"""
		raise NotImplementedError("")








class OrbitalBasisSetStandard(OrbitalBasisSetBase):

	def __init__(self, orbitalBasisFunctions, basisName):
		""" Initializer
		
		Args:
			orbitalBasisFunctions: (iter of OrbitalBasisFunctionBase) Each represents a single basis function and is capable of calculating its radial values
			basisName: (str) Alias used to identify this basis set
				 
		"""
		self._basisName = basisName
		self._orbitalBasisFunctions = orbitalBasisFunctions


	@property
	def orbitalBasisFunctions(self):
		return self._orbitalBasisFunctions

	@property
	def basisName(self):
		return self._basisName

	def getBasisFunctionsWithLVal(self, lVal):
		outList = list()
		for x in self.orbitalBasisFunctions:
			if x.lVal==lVal:
				outList.append(x)
		return outList



class GauSumOrbitalBasisFunction(OrbitalBasisFunctionBase):
	
	def __init__(self, nVal, lVal, gauFunct):
		""" Initializer
		
		Args:
			nVal: (int) Used to differentiate basis functions with the same angular momentum. No deeper meaning.
			lVal: (int) Angular momentum for this orbital (s=0,p=1,d=2)
			gauFunct: (GauPrimBase object) This has a evalFunctAtDists(self, distances) method, which gives the total radial component that far from origin. This will generally be a composite object
				 
		"""
		self._nVal = nVal
		self._lVal = lVal
		self.gauFunct = gauFunct

	def getRadialValsAtDists(self, xVals, multByRl=False):
		""" Get the radial part of the basis function at distances given in xVals
		
		Args:
			xVals: (float iter) List of distances from origin (i.e. from the atom centre)
			multByRl: (Bool, default=False) Whether to multiply the result by r^{l}. This should give the true radial distribution of the orbital; we actually represent the "radial" part in a way where some angular momenta terms are incorporated for convenience. Setting to True gives the same as the values in the plato *.bas files
	 
		Returns
			radialVals: (float iter) Radial values of the wavefunction
		"""
		outYVals = self.gauFunct.evalFunctAtDists(xVals)
		if multByRl:
			lVal = self.lVal
#			normFactor = math.sqrt( (4*math.pi) / ((2*lVal)+1) )
			normFactor = 1.0
			for idx,val in enumerate(outYVals):
				outYVals[idx] *= (xVals[idx]**lVal)*normFactor

		return outYVals


	@property
	def nVal(self):
		return self._nVal

	@nVal.setter
	def nVal(self,val):
		self._nVal = val

	@property
	def lVal(self):
		return self._lVal


