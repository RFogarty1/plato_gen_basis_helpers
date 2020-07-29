
import copy
import itertools as it
import types

import numpy as np

from ..shared import surfaces as surfHelp
from . import base_flow as baseFlow



class StackingFaultWorkflowTwoStructs(baseFlow.BaseLabelledWorkflow):
	"""Workflow for calculating a stacking fault energy when given the final and initial structures
	"""

	def __init__(self, perfectStructObj, stackFaultStructObj):
		""" Initializer
		
		Args:
			perfectStructObj: (CalcMethod object) Object used for calculating energy of the system without a stacking fault
			stackFaultStructObj: (CalcMethod object) Object used for calculating energy of the system with a stacking fault

		"""
		self.perfectStructObj = perfectStructObj
		self.stackFaultStructObj = stackFaultStructObj
		self.eType = "electronicTotalE"
		self._output = types.SimpleNamespace( **{k:None for k in self.namespaceAttrs} )
		self._writeInpFiles()

	@property	
	def namespaceAttrs(self):
		return ["stackFaultEnergy"]

	@property
	def output(self):
		return [self._output]

	@property
	def preRunShellComms(self):
		allComms = [self.perfectStructObj.runComm, self.stackFaultStructObj.runComm]
		return allComms


	def _writeInpFiles(self):
		self.perfectStructObj.writeFile()
		self.stackFaultStructObj.writeFile()

	def run(self):
		ePerAreaPerfect = self._getEnergyPerAreaForPerfectStruct()
		ePerAreaFaulted = self._getEnergyPerAreaForStackFault()
		self._output.stackFaultEnergy = ePerAreaFaulted-ePerAreaPerfect

	def _getEnergyPerAreaForStackFault(self):
		energyFaulted = getattr(self.stackFaultStructObj.parsedFile.energies, self.eType)
		surfArea = _getABSurfaceAreaFromParsedFile(self.stackFaultStructObj.parsedFile)
		return energyFaulted/surfArea

	def _getEnergyPerAreaForPerfectStruct(self):
		energyPerfect = getattr(self.perfectStructObj.parsedFile.energies, self.eType)
		surfAreaPerfect = _getABSurfaceAreaFromParsedFile(self.perfectStructObj.parsedFile)
		return energyPerfect/surfAreaPerfect




#Probably useless junk
class StackingFaultWorkflow(baseFlow.BaseLabelledWorkflow):
	""" Workflow for calculating stacking fault energies including unstable and stable. """


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
		self._writeInpFiles()

	def _writeInpFiles(self):
		self.perfectStructCalcObj.writeFile()
		[x.writeFile() for x in self.calcObjs]

	@property
	def preRunShellComms(self):
		allComms = [self.perfectStructCalcObj.runComm]
		for x in self.calcObjs:
			allComms.append( x.runComm )
		return allComms

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
	inpCell = copy.deepcopy(parsedFile.unitCell)
	inpCell.fractCoords = [[0,0,0,"x"]] #Need some kind of fractional co-ordinates to create the surface object
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

	@property
	def fitFunct(self):
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

	def __init__(self, order, searchRangeIntrinsic=None, searchRangeUnstable=None):
		""" Initialiser
		
		Args:
			order: (int) the highest power of to use in the polynomial fit. 0 means a constant value, 1 means a linear fit, 2 means a quadratic etc.
			searchRangeIntrinsic: (len 2 float iter) [min,max] displacements to look for a minimum in (defaults to whole range supplied) 
			searchRangeUnstable: (len 2 float iter) [min,max] displacmenets to look for a maximum in (defaults to whole range supplied)	 
		"""
		self.order = order
		self.stepValGetStackEnergy = 0.001
		self.stepValOutputDists = 0.001
		self.searchRangeIntrinsic = searchRangeIntrinsic
		self.searchRangeUnstable = searchRangeUnstable


	def _getFitCoeffs(self, dispVals, structFaultVals):
		outCoeffs = np.polyfit(dispVals, structFaultVals, self.order)
		return outCoeffs

	def _getPolyFunctionFromCoeffs(self, coeffs):
		revCoeffs = [x for x in reversed(coeffs)] #Need to get low to high order (i.e x^{0} coeff first)
		outFunct = np.polynomial.polynomial.Polynomial( revCoeffs )
		return outFunct

#	def _getUnstableFaultEnergy(self, fitFunct):
#		xRange = np.arange(0, 1, self.stepValGetStackEnergy)
#		fitVals = [fitFunct(x) for x in xRange]
#		return max(fitVals)

	def _getUnstableFaultEnergy(self, actDispVals, actStackVals, fitDispVals, fitStackVals):
		if self.searchRangeUnstable is None:
			searchRange = [ min( [min(actDispVals),min(fitDispVals)] ), max( [max(actDispVals), max(fitDispVals)] ) ] 
		else:
			searchRange = self.searchRangeUnstable

		actDispFiltered, actStackFiltered = self._getFilteredDispAndStackValsWithinRange(actDispVals, actStackVals, searchRange)
		fitDispFiltered, fitStackFiltered = self._getFilteredDispAndStackValsWithinRange(fitDispVals, fitStackVals, searchRange)

		maxFromActVals = max(actStackFiltered)
		maxFromDispVals = max(fitStackFiltered)
		return max( [maxFromActVals, maxFromDispVals] )


	def _getIntrinsicFaultEnergy(self, actDispVals, actStackVals, fitDispVals, fitStackVals):
		if self.searchRangeIntrinsic is None:
			searchRange = [ min( [min(actDispVals),min(fitDispVals)] ), max( [max(actDispVals), max(fitDispVals)] ) ] 
		else:
			searchRange = self.searchRangeIntrinsic

		actDispFiltered, actStackFiltered = self._getFilteredDispAndStackValsWithinRange(actDispVals, actStackVals, searchRange)
		fitDispFiltered, fitStackFiltered = self._getFilteredDispAndStackValsWithinRange(fitDispVals, fitStackVals, searchRange)

		minFromActVals = min(actStackFiltered)
		minFromDispVals = min(fitStackFiltered)
		return min( [minFromActVals, minFromDispVals] )


	def _getFilteredDispAndStackValsWithinRange( self, dispVals, stackVals, searchRange ):
		outDisp, outStack = list(), list()
		minDisp, maxDisp = min(searchRange), max(searchRange)
		for dVal, sVal in it.zip_longest(dispVals, stackVals):
			if (dVal >= minDisp) and (dVal <= maxDisp):
				outDisp.append( dVal )
				outStack.append( sVal ) 
		return outDisp, outStack

	def _getOutputDistVals(self, dispVals):
		minVal = min(dispVals)
		maxVal = max(dispVals)
		fitDispVals = np.arange(minVal,maxVal,self.stepValOutputDists)
		return fitDispVals

	def _calcGoodnessOfFitVals(self, stackFaultVals, fitValsAtStackFault):
		absRelErrors = list()
		for act,fit in it.zip_longest(stackFaultVals, fitValsAtStackFault):
			if abs(act) < 1e-10:
				pass
			else:
				currError = abs(act-fit)/act
				absRelErrors.append(currError)

		try:
			outVal = sum(absRelErrors) / len(absRelErrors)
		except ZeroDivisionError:
			outVal = None

		return outVal

	def fitStackingFaultEnergies(self, dispVals, structFaultVals):
		#1) Get coefficients for the fit
		outCoeffs = self._getFitCoeffs(dispVals, structFaultVals)
		polyFunct = self._getPolyFunctionFromCoeffs(outCoeffs)

		#Get everything needed for the fitResult object
		fitResDict = dict()
		fitResDict["fitParams"] = outCoeffs
		fitResDict["fitFaultValsAtInputDisps"] = [polyFunct(x) for x in dispVals]

		fitDispVals = self._getOutputDistVals(dispVals)
		fitFaultVals = [polyFunct(x) for x in fitDispVals]

		fitResDict["fitDispVals"], fitResDict["fitFaultVals"] = fitDispVals, fitFaultVals

		fitResDict["intrinsicFaultEnergy"] = self._getIntrinsicFaultEnergy( dispVals, structFaultVals, fitDispVals, fitFaultVals )
		fitResDict["unstableFaultEnergy"] = self._getUnstableFaultEnergy( dispVals, structFaultVals, fitDispVals, fitFaultVals )


		fitResDict["goodnessOfFit"] = self._calcGoodnessOfFitVals(structFaultVals, fitResDict["fitFaultValsAtInputDisps"])
		fitResDict["fitFunct"] = polyFunct

		return types.SimpleNamespace(**fitResDict)


