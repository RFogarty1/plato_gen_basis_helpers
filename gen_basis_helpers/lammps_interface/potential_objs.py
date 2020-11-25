
import collections


class LammpsPotential():
	""" Should contain all info needed to define a potential for LAMMPS

	"""
	def __init__(self, globalPotOpts, potParams, mapPotParams):
		""" Initializer
		
		Args:
			globalPotOpts: (GlobalPotOptions object) Contains options including the type of potential to use
			potParams: (PotentialParamsLJAndHarmonicBase) Contains specific parameters (e.g. force constants) 
			mapPotParams: (f(potParams)->commandDict) Generally callable class. Used to translate potParams into a command dictionary.

		"""
		self.globalPotOpts = globalPotOpts
		self.potParams = potParams
		self.mapPotParams = mapPotParams


	def _getCommandsFromGlobalPotOpts(self):
		commAttrCombos = collections.OrderedDict( [ ["pair_style","pairStyle"], ["kspace_style","kSpaceStyle"], ["bond_style","bondStyle"],
		                                            ["angle_style","angleStyle"], ["dihedral_style","dihedralStyle"], ["improper_style","improperStyle"] ] )
		outDict = collections.OrderedDict()
		for key,attr in commAttrCombos.items():
			val = getattr(self.globalPotOpts,attr)
			if val is not None:
				outDict[key] = val
		return outDict

	@property
	def commandDict(self):
		outDict = self._getCommandsFromGlobalPotOpts()
		extraOptsDict = self.mapPotParams(self.potParams)
		outDict.update(extraOptsDict)
		return outDict

#TODO: Make a function to just create a normal one of these for TIP3P
class GlobalPotOptions():

	def __init__(self, pairStyle=None, kSpaceStyle=None, bondStyle=None, angleStyle=None,
	             dihedralStyle=None, improperStyle=None):
		""" Initializer
		
		Args:
			pairStyle: (str) Passed to pair_style command
			kSpaceStyle: (str) Passed to kspace_style command
			bondStyle: (str) Passed to bond_style command (e.g. harmonic,none)
			angleStyle: (str) Passed to angle_style command (e.g. harmonic,none)
			dihedralStyle: (str) Passed to dihedral_style command (e.g. harmonic,none)
			improperStyle: (str) Passed to improper_style command (e.g. harmonic,none)
		"""
		self.pairStyle = pairStyle
		self.kSpaceStyle = kSpaceStyle
		self.bondStyle = bondStyle
		self.angleStyle = angleStyle
		self.dihedralStyle = dihedralStyle
		self.improperStyle = improperStyle



class PotentialParamsLJAndHarmonicBase():

	@property
	def bondPots(self):
		""" Iter of InternalMoleculeHarmonicPotential objects defining potentials between bonded atoms"""
		return list()

	@property
	def anglePots(self):
		""" Iter of InternalMoleculeHarmonicPotential objects defining potentials from angular distortions between bonded atoms"""
		return list()

	@property
	def lennardJonesPots(self):
		return list()


#TODO
class MapLJAndHarmonicParamsToCommandDict():
	""" Callable class for getting a LAMMPS command dictionary from a PotentialParamsLJAndHarmonicBase. See self.mapParamsToCommandDict for the interface

	"""

	def __init__(self, eleToParamsMap, numbFmt="{:.4f}"):
		""" Iniitalizer
		
		Args:
			eleToParamsMap: (TypeMaps object) This contains dictionaries mapping elements to type indices, as well as (optionally) bonds and angles to type indices. These optional ones will only be looked in when the potential object contains the relevant potentaisl (e.g. if not harmonic angle potentials are defined, the mapper wont look in eleToParamsMap.angleToTypeIdx 
			numbFmt: (Str) Format string for writing the params to file
	 
		"""
		self.eleToParamsMap = eleToParamsMap
		self.numbFmt = numbFmt

	def mapParamsToCommandDict(self, paramObj):
		""" Gets LAMMPS commands from a potential params object
		
		Args:
			paramObj: (PotentialParamsLJAndHarmonicBase object) Contains elemental parameters info 
				 
		Returns
			outDict: (OrderedDict) Each key is a command for a LAMMPS script file, while values are the strings associated with the command
	 
		"""
		outDict = collections.OrderedDict()

		#1st get the pair coeffs
		pairCoeffStr = ""
		for idx,ljObj in enumerate(paramObj.lennardJonesPots):
			currTypeIndices = [self.eleToParamsMap.eleToTypeIdx[x] for x in ljObj.elements]
			if idx==0:
				fmtStr = "{} {} " + self.numbFmt + " " + self.numbFmt 
			else:
				fmtStr = "\npair_coeff " + "{} {} " + self.numbFmt + " " + self.numbFmt 
			pairCoeffStr += fmtStr.format(*currTypeIndices, ljObj.epsilon, ljObj.sigma)

		if pairCoeffStr != "":
			outDict["pair_coeff"] = pairCoeffStr

		#Get the bond coeffs
		bondCoeffStr = self._getInternalHarmonicCoordStr(paramObj.bondPots, self.eleToParamsMap.bondToTypeIdx,"bond_coeff")
		if bondCoeffStr != "":
			outDict["bond_coeff"] = bondCoeffStr

		#Get the angle coeffs
		angleCoeffStr = self._getInternalHarmonicCoordStr(paramObj.anglePots, self.eleToParamsMap.angleToTypeIdx,"angle_coeff")
		if angleCoeffStr != "":
			outDict["angle_coeff"] = angleCoeffStr

		#TODO: I should really try to get dihedrals/improper to make this general; but probably will never need that so....

		return outDict


	def _getInternalHarmonicCoordStr(self, potObjs, mapToTypeIdxDict, keyVal):
		outStr = ""
		for idx,potObj in enumerate(potObjs):
			currTypeIdx = mapToTypeIdxDict[tuple(potObj.elements)]
			if idx==0:
				fmtStr = "{} " + self.numbFmt + " " + self.numbFmt
			else:
				fmtStr = "\n"+ keyVal + " " +  "{} " + self.numbFmt + " " + self.numbFmt
			outStr += fmtStr.format(currTypeIdx, potObj.forceConstant, potObj.eqmVal)

		return outStr

	def __call__(self, paramObj):
		return self.mapParamsToCommandDict(paramObj)


class PriceTIP3PPotential(PotentialParamsLJAndHarmonicBase):
	
	def __init__(self, ljOO=None, ljOH=None, ljHH=None, eqmLengthOH=0.9572,
	             forceConstOH=450, eqmAngleHOH=104.52, forceConstHOH=55):
		""" Initializer
		
		Args:
			ljOO: lennard jones params for OO [epsilon,sigma]. Default is [0.102,3.188]
			ljOH: lennard jones params for OH [epsilon,sigma]. Default is [0,0]
			ljHH: lennard jones params for HH [epsilon, sigma]. Default is [0,0]
			eqmLengthOH: Length of eqm OH bondlength 
			forceConstOH: Force constant for OH bond
			eqmAngleHOH: Angle (Degrees) of the equilibrium H-O-H bond
			forceConstHOH: Force constant for the H-O-H bond

		"""
		self.ljOO = [0.102,3.188] if ljOO is None else ljOO
		self.ljOH = [0,0] if ljOH is None else ljOH
		self.ljHH = [0,0] if ljHH is None else ljHH
		self.eqmLengthOH = eqmLengthOH
		self.forceConstOH = forceConstOH
		self.eqmAngleHOH = eqmAngleHOH
		self.forceConstHOH = forceConstHOH

	#Empty list by default would make it easier to extend
	@property
	def bondPots(self):
		outPot = InternalMoleculeHarmonicPotential(["O","H"], self.eqmLengthOH, self.forceConstOH)
		return [outPot]

	@property
	def anglePots(self):
		outPot = InternalMoleculeHarmonicPotential(["H","O","H"], self.eqmAngleHOH, self.forceConstHOH)
		return [outPot]

	@property
	def lennardJonesPots(self):
		potOO = LennardJonesPotential(["O","O"], *self.ljOO)
		potOH = LennardJonesPotential(["O","H"], *self.ljOH)
		potHH = LennardJonesPotential(["H","H"], *self.ljHH)
		return [potOO, potOH, potHH]


class Jorg1983TIP3PPotential(PriceTIP3PPotential):

	def __init__(self, ljOO=None, ljOH=None, ljHH=None, eqmLengthOH=0.9572,
	             forceConstOH=450, eqmAngleHOH=104.52, forceConstHOH=55):
		""" Initializer
		
		Args:
			ljOO: lennard jones params for OO [epsilon,sigma]. Default is [0.1521,3.1507]
			ljOH: lennard jones params for OH [epsilon,sigma]. Default is [0,0]
			ljHH: lennard jones params for HH [epsilon, sigma]. Default is [0,0]
			eqmLengthOH: Length of eqm OH bondlength 
			forceConstOH: Force constant for OH bond
			eqmAngleHOH: Angle (Degrees) of the equilibrium H-O-H bond
			forceConstHOH: Force constant for the H-O-H bond

		"""
		self.ljOO = [0.1521,3.1507] if ljOO is None else ljOO
		self.ljOH = [0,0] if ljOH is None else ljOH
		self.ljHH = [0,0] if ljHH is None else ljHH
		self.eqmLengthOH = eqmLengthOH
		self.forceConstOH = forceConstOH
		self.eqmAngleHOH = eqmAngleHOH
		self.forceConstHOH = forceConstHOH


class InternalMoleculeHarmonicPotential():

	def __init__(self, elements, eqmVal, forceConstant):
		""" Initializer
		
		Args:
			elements: (len-x iter) Containing the elements in the bond or angle
			eqmVal: (float) Length of the equilibrium bond or angle of eqm angle etc
			forceConstant: (float) Force constant for the bond/angle (higher=harder to stretch/bend)
 
		"""
		self._eqTol = 1e-5
		self.elements = list(elements)
		self.eqmVal = eqmVal
		self.forceConstant = forceConstant

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)
		numericalAttrs = ["forceConstant", "eqmVal"]
		if self.elements != other.elements:
			return False

		for attr in numericalAttrs:
			if abs( getattr(self,attr) - getattr(other,attr) ) > eqTol:
				return False 

		return True


class LennardJonesPotential():

	def __init__(self, elements, epsilon, sigma):
		""" Initializer
		
		Args:
			elements: (len-2 iter) Contains the elements this interaction is between
			epsilon: Epsilon in 4\epsilon[ (\sigma/r)^12 - (\sigma/r)^6 ]
			sigma: Sigma in 4\epsilon[ (\sigma/r)^12 - (\sigma/r)^6 ]
 
		"""
		self._eqTol = 1e-5
		self.elements = list(elements)
		self.epsilon = epsilon
		self.sigma = sigma 

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)
		numericalAttrs = ["epsilon","sigma"]
		if self.elements != other.elements:
			return False

		for attr in numericalAttrs:
			if abs( getattr(self,attr) - getattr(other,attr) ) > eqTol:
				return False 

		return True


def getElementToChargesDictPriceTIP3P():
	return {"O":-0.830, "H":0.415}

def getElementToChargesDictJorgensen1983():
	return {"O":-0.834, "H":0.417}


