

import itertools as it
import types

import numpy as np

from ..shared import surfaces as surfHelp
from . import base_flow as baseFlow


class StackingFaultWorkflow(baseFlow.BaseLabelledWorkflow):
	""" Workflow for calculating stacking fault energies """


	def __init__(self, perfectStructCalcObj, calcObjs, dispVals, fitterObj, energyType=None):
		""" Initializer
		
		Args:
			perfectCalcObj: (CalcMethod object) Represents the perfect structure (with no stacking fault). Used to define the zero energy mark and isnt ALWAYS the same as displacement=0 case
			calcObjs: (iter of CalcMethod objects)
			dispVals: (iter of floats) These need to be in the correct order for calcObjs (the calcObj representing displacement dispVals[x] needs to be at calcObjs[x])
			fitterObj: (StackingFaultFitterBase object) Carries out the fit of displacement factor vs energies; needed to get unstable stacking fault at minimum
				 
		"""
		self.perfectStructCalcObj = perfectStructCalcObj
		self.calcObjs = list(calcObjs)
		self.dispVals = list(dispVals)
		self.energyType = "electronicTotalE" if energyType is None else energyType
		self.fitterObj = fitterObj
		self._output = types.SimpleNamespace( **{k:None for k in self.namespaceAttrs} )

	@property	
	def namespaceAttrs(self):
		return ["displacements","totalEnergiesRaw", "energiesPerSurfaceAreaRaw", "totalEForPerfectStruct",
		        "ePerAreaForPerfectStruct","stackFaultEnergiesRaw", "fitResult"]

	@property
	def output(self):
		return self._output

	def _getTotalEnergies(self):
		totalEnergies = list()
		for dVal, calcObj in it.zip_longest(self.dispVals,self.calcObjs):
			parsedFile = calcObj.parsedFile
			currEnergy = getattr(parsedFile.energies,self.energyType)
			totalEnergies.append( currEnergy )
		return totalEnergies

	def _getEnergiesPerSurfaceArea(self):
		totalEnergies = self._getTotalEnergies()
		surfaceAreas = [_getABSurfaceAreaFromParsedFile(x.parsedFile) for x in self.calcObjs]
		ePerSurfaceArea = [ e/surfArea for e,surfArea in it.zip_longest(totalEnergies,surfaceAreas) ]
		return ePerSurfaceArea

	def _getTotalEnergyPerfectStruct(self):
		return getattr(self.perfectStructCalcObj.parsedFile.energies, self.energyType)

	def _getEPerSurfaceAreaForPerfectStruct(self):
		return self._getTotalEnergyPerfectStruct() / _getABSurfaceAreaFromParsedFile(self.perfectStructCalcObj.parsedFile)

	def _getStackFaultEnergies(self):
		refEnergyPerSurfArea = self._getEPerSurfaceAreaForPerfectStruct()
		otherEnergyPerSurfArea = self._getEnergiesPerSurfaceArea()
		return [x-refEnergyPerSurfArea for x in otherEnergyPerSurfArea]

	def run(self):
		self._output.totalEForPerfectStruct = self._getTotalEnergyPerfectStruct()
		self._output.ePerAreaForPerfectStruct = self._getEPerSurfaceAreaForPerfectStruct()
		self._output.totalEnergiesRaw = self._getTotalEnergies()
		self._output.energiesPerSurfaceAreaRaw = self._getEnergiesPerSurfaceArea()
		self._output.displacements = self.dispVals
		self._output.stackFaultEnergiesRaw = self._getStackFaultEnergies()
		self._output.fitResult = self.fitterObj.fitStackingFaultEnergies( self.dispVals, self._output.stackFaultEnergiesRaw )

def _getABSurfaceAreaFromParsedFile(parsedFile):
	inpCell = parsedFile.unitCell
	nLayers, lenAbsVac = 1,0
	surfObj = surfHelp.GenericSurface(inpCell, nLayers, lenAbsoluteVacuum=lenAbsVac)
	return surfObj.surfaceArea



class StackingFaultFitResultsBase():

	@property
	def goodnessOfFit(self):
		raise NotImplementedError("")

	@property
	def unstableFaultEnergy(self):
		raise NotImplementedError("")

	@property
	def intrinsicFaultEnergy(self):
		raise NotImplementedError("")

	@property
	def fitFaultVals(self):
		raise NotImplementedError("")

	@property
	def fitDispVals(self):
		raise NotImplementedError("")

	@property
	def fitFaultValsAtInputDisps(self):
		raise NotImplementedError("")

	@property
	def fitParams(self):
		raise NotImplementedError("")

class StackingFaultFitterBase():
	""" Class used to carry out fits to stacking fault data
 
	"""

	def fitStackingFaultEnergies(self, dispVals, structFaultVals):
		""" Fits the input stacking fault energies.
		
		Args:
			dispVals: (iter of floats) Set of values, generally between 0 and 1, indicating the degree of displacement for each structure
			structFaultVals: (iter of floats) Each entry is energy per surface area minus value for the perfect structure
				 
		Returns
			fitResult: StackingFaultFitResultsBase object, contains results of the fit which should include intrinsic and unstable stacking fault energies
	 
		"""
		raise NotImplementedError("")


class StackingFaultFitterPolyStandard(StackingFaultFitterBase):

	def __init__(self, order):
		""" Initialiser
		
		Args:
			order: (int) the highest power of to use in the polynomial fit. 0 means a constant value, 1 means a linear fit, 2 means a quadratic etc. 
				 
		"""
		self.order = order
		self.stepValGetStackEnergy = 0.001
		self.stepValOutputDists = 0.001


	def _getFitCoeffs(self, dispVals, structFaultVals):
		outCoeffs = np.polyfit(dispVals, structFaultVals, self.order)
		return outCoeffs

	def _getPolyFunctionFromCoeffs(self, coeffs):
		revCoeffs = [x for x in reversed(coeffs)] #Need to get low to high order (i.e x^{0} coeff first)
		outFunct = np.polynomial.polynomial.Polynomial( revCoeffs )
		return outFunct

	def _getUnstableFaultEnergy(self, fitFunct):
		xRange = np.arange(0, 1, self.stepValGetStackEnergy)
		fitVals = [fitFunct(x) for x in xRange]
		return max(fitVals)

	def _getOutputDistVals(self, dispVals):
		minVal = min(dispVals) if min(dispVals) < 0 else 0
		maxVal = max(dispVals) if max(dispVals) > 1 else 1
		fitDispVals = np.arange(minVal,maxVal,self.stepValOutputDists)
		return fitDispVals

	def _calcGoodnessOfFitVals(self, stackFaultVals, fitValsAtStackFault):
		absRelErrors = [ abs( (act-fit)/act ) for act,fit in it.zip_longest(stackFaultVals, fitValsAtStackFault) ]
		return sum(absRelErrors) / len(absRelErrors)

	def fitStackingFaultEnergies(self, dispVals, structFaultVals):
		#1) Get coefficients for the fit
		outCoeffs = self._getFitCoeffs(dispVals, structFaultVals)
		polyFunct = self._getPolyFunctionFromCoeffs(outCoeffs)

		#Get everything needed for the fitResult object
		fitResDict = dict()
		fitResDict["fitParams"] = outCoeffs
		fitResDict["fitFaultValsAtInputDisps"] = [polyFunct(x) for x in dispVals]
		fitResDict["intrinsicFaultEnergy"] = polyFunct(1)
		fitResDict["unstableFaultEnergy"] = self._getUnstableFaultEnergy(polyFunct)
		fitResDict["fitDispVals"] = self._getOutputDistVals(dispVals)
		fitResDict["fitFaultVals"] = [polyFunct(x) for x in fitResDict["fitDispVals"]]
		fitResDict["goodnessOfFit"] = self._calcGoodnessOfFitVals(structFaultVals, fitResDict["fitFaultValsAtInputDisps"])

		return types.SimpleNamespace(**fitResDict)


