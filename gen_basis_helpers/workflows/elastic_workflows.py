
import itertools as it
import types
from . import base_flow as baseFlow
from ..shared import unit_convs as unitConvs

import plato_pylib.utils.elastic_consts as elasticHelp






class StressStrainWorkflow(baseFlow.BaseLabelledWorkflow):
	""" Workflow representing calculation of a single stress-strain curve
	"""


	def __init__(self, calcObjs, strainCoeffs, strain, label=None, eType="electronicTotalE"):
		""" Initializer
		
		Args:
			calcObjs: (iter of CalcMethod objects) Each of thes objects represents a single point energy calculation
			strainCoeffs: (iter of float) Each represents the strain applied to the eqm. geom. Length should be the same as calcObjs, and the ordering of values should be linked
			strain: (CrystalStrain object) This contains information on the unit-strain (i.e. the strain applied if strainCoeff=1.0)
			eType: (Str) Type of energy we want; Corresponds to attributes on plato_pylib Energies object. Default is total electronic energy
			
	
		"""
		self.calcObjs = list(calcObjs)
		self.strainCoeffs = list(strainCoeffs)
		self.strain = strain
		self.eType = eType
		self._output = [ types.SimpleNamespace( **{k:None for k in self.namespaceAttrs[0]} ) ]

		self._writeInpFiles()

	def _writeInpFiles(self):
		for x in self.calcObjs:
			x.writeFile()

	@property
	def preRunShellComms(self):
		outComms = list()
		for x in self.calcObjs:
			outComms.append( x.runComm ) #Each command is a single string
		return outComms

	@property
	def namespaceAttrs(self):
		return [["fitFunct","secondDeriv","actVals","strainVsEnergy"]]

	@property
	def output(self):
		""" Iter of objects (length 1 for this leaf-class) representing the output of the StressStrainWorkflow

		Attrs:
			actVals: (nx2 iter) Calculated strain-stress values for strain values corresponding to self.strainCoeffs. Units should be GPa
			strainVsEnergy: (nx2 iter) x values are self.strainCoeffs, y values are the total energies from the calculations. Units should be eV
			fitFunct: (function, accepts iter of x-values) When given iter of x-values(strains) it returns the stresses expected from a quadratic fit to the strain vs stress curve nearby strain=0
			secondDeriv: The second-derivate of the stress-strain curve; calculated by fitting an x**2 parabola near the zero-strain point

		"""
		return self._output

	def run(self):
		self._output[0].strainVsEnergy = self._getStrainVsEnergyFromCalcs()
		self._output[0].actVals = self._getStrainVsStressFromCalcs()
		fitObj = self._getFitObjFromStressStrain( self.output[0].actVals )
		self._output[0].fitFunct = fitObj.getFittedValuesForXVals
		self._output[0].secondDeriv = fitObj.secondDeriv

	def _getStrainVsEnergyFromCalcs(self):
		allEnergies = list()
		for x in self.calcObjs:
			allEnergies.append( getattr(x.parsedFile.energies,self.eType) )
		minEnergy = min(allEnergies)
		allEnergies = [x-minEnergy for x in allEnergies]
		return [ [x,y] for x,y in it.zip_longest(self.strainCoeffs, allEnergies) ] #Need to conv to stress + get units

	def _getStrainVsStressFromCalcs(self):
		allStress = list()
		for x in self.calcObjs:
			parsedFile = x.parsedFile
			energy = getattr(parsedFile.energies,self.eType)
			volume = parsedFile.unitCell.volume
			allStress.append( energy/volume )

		allStress = [self._applyInpUnitsToGPaConversionFactor(x) for x in allStress]
		strainVsStress = [ [x,y] for x,y in it.zip_longest(self.strainCoeffs, allStress) ]
		return strainVsStress

	def _applyInpUnitsToGPaConversionFactor(self, stressVal):
		evToJoule = unitConvs.EV_TO_JOULE
		bohrToMetre = unitConvs.BOHR_TO_METRE
		pascalToGpa = unitConvs.PASCAL_TO_GPA
		elasticConvFactor = (evToJoule / (bohrToMetre**3))*pascalToGpa
		return stressVal*elasticConvFactor

	def _getFitObjFromStressStrain(self, stressStrain):
		return elasticHelp.polyFitAndGetSecondDeriv(stressStrain)

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



