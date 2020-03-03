
import collections
import itertools as it
import types

import numpy as np

from . import base_flow as baseFlow
from ..shared import unit_convs as unitConvs
from ..shared import creator_resetable_kwargs as baseCreator

import plato_pylib.utils.elastic_consts as elasticHelp


class HcpElasticWorkflowCreator(baseCreator.CreatorWithResetableKwargsTemplate):
	""" Factory for creating HcpElasticConstantWorkflow objects
	"""
	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs)
	registeredKwargs.add("baseGeom")
	registeredKwargs.add("strainValues")
	registeredKwargs.add("structType") #Used to get the strains we need
	registeredKwargs.add("creator")

	def _createFromSelf(self):
		raise NotImplementedError("")



class HcpElasticConstantsWorkflow(baseFlow.BaseLabelledWorkflow):
	"""Workflow representing calculations used to get elastic constants for a hcp crystal

	"""

	def __init__(self, stressStrainFlows, label=None):
		""" Initializer
		
		Args:
			stressStrainFlows: (iter of 5 StressStrainWorkflow objects) Each of these objects needs to represent one required stress strain curves. In terms of the basis functions \eps_{1} to \eps_{6} we need:
			                  a) \eps_{3}
			                  b) \eps_{1} + \eps_{2}
			                  c) \eps_{1} + \eps_{2} + \eps_{3}
			                  d) 2( \eps_{4} + \eps_{5} )
			                  e) 2( \eps_{6} )

		Raises:
			ValueError: If not all strains are correct (including if the number of stress/strains is wrong) 
		"""
		self.stressStrainFlows = self._getStrainFlowInOrderFromInputList( stressStrainFlows )

		self._output = [ types.SimpleNamespace( **{k:None for k in self.namespaceAttrs[0]} ) ]

	def _getStrainFlowInOrderFromInputList(self, inpList):
		if len(inpList) != 5:
			raise ValueError("Need 5 strains to calculate Hcp elastic constants but {} given".format(len(inpList)))

		expStrainObjs = getRequiredStrainObjsForStructType("hcp") #This actually defines the ordering aswell. Unit-tests (at time of writing) will catch if the order changes though (i.e. will fail if the order of strains in this list changes)

		outList = list()
		for x in inpList:
			currObj = x.strain
			if any([currObj==expObj for expObj in expStrainObjs]):
				outList.append(x)
			else:
				raise ValueError("{} is an invalid strain for HcpElasticConstantsWorkflow".format(currObj.toStr()))

		#Make sure all strains are different
		allStrainObjs = [x.strain for x in inpList]
		for idx,strain in enumerate(allStrainObjs):
			otherStrains = [x for idxB,x in enumerate(allStrainObjs) if idxB!=idx]
			if any([x==strain for x in otherStrains]):
				raise ValueError("Duplicate strains found")

		#Make sure the strains are always in the correct order
		orderedOutList = list()
		for expObj in expStrainObjs:
			for currObj in outList:
				if expObj==currObj.strain:
					orderedOutList.append(currObj)

		return orderedOutList

	@property
	def namespaceAttrs(self):
		return [["elasticConsts","stressStrainData"]]

	@property
	def preRunShellComms(self):
		outList = list()
		for x in self.stressStrainFlows:
			outList.extend( x.preRunShellComms )
		return outList

	@property
	def output(self):
		""" Iter of objects (length 1 for this leaf-class) representing the output of the HcpElasticConstantsWorkflow
		
		Attrs:
			elasticConsts: (OrderedDict) Keys are elastic constant labels (11,12,33,44 corresponding to c_{11}, c_{12} etc.) and values are calculated elastic constants
			stressStrainData: (iter) Each entry is the .output of one of the stress-strain workflows used. The initializer of this class gives the strains used in the same order as output here
				
		"""
		return self._output


	def run(self):
		self.output[0].elasticConsts = self._getElasticDict()
		self.output[0].stressStrainData = [x.output[0] for x in self.stressStrainFlows]

	def _getElasticDict(self):
		secondDerivs = [x.output[0].secondDeriv for x in self.stressStrainFlows]
		keyOrder = ["11","12","13","33","44"]
		coeffMatrix = np.array([ [0, 0,0,1,0],
		                         [2, 2,0,0,0],
		                         [2, 2,4,1,0],
		                         [0, 0,0,0,8],
		                         [2,-2,0,0,0] ])
		elasticVals = np.linalg.inv(coeffMatrix) @ np.array(secondDerivs)

		return collections.OrderedDict( [[key,elastic] for key,elastic in it.zip_longest(keyOrder,elasticVals)] )



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
			diffs = [abs(x-y) for x,y in it.zip_longest(self.strainVals, other.strainVals)]
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


def getRequiredStrainObjsForStructType(structKey):
	keyToFunctDict = {"hcp": _getHcpStrainObjects}
	return keyToFunctDict[structKey]()



def _getHcpStrainObjects():
	expStrains = [ [0,0,1,0,0,0], 
	               [1,1,0,0,0,0],
	               [1,1,1,0,0,0],
	               [0,0,0,2,2,0],
	               [0,0,0,0,0,2] ]

	expStrainObjs = [CrystalStrain(x) for x in expStrains]
	return expStrainObjs


def _getUnitStrainMatrix(matrixNumb:int):
	strainParam = 1.0
	return elasticHelp._STRAIN_MATRIX_DICT[matrixNumb](1.0)



