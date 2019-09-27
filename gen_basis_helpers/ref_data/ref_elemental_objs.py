
import itertools as it
from plato_pylib.plato.plato_paths import PlatoPathDescriptor, PlatoModelFolders


class AngMomMissing(KeyError):
	pass

class RefElementalDataBase():
	""" Class for holding reference data for a pure-elemental system """


	@property
	def modelFiles(self):
		""" plato pylib PlatoModelFolders object, handles paths to model datafolders """
		raise NotImplementedError


	def getExptGeom(self, key):
		""" Get UnitCell for experimental geometry
		
		Args:
			key: structure key; hcp/bcc/fcc are standard options. Likely set to default/single possible value in subclasses 
				
		Returns
			outStruct: plato_pylib UnitCell object containing the optimised geometry
		
		"""
		raise NotImplementedError


	def getPlaneWaveGeom(self, key):
		""" Get UnitCell for plane-wave optimised geometry
		
		Args:
			key: structure key; hcp/bcc/fcc are standard options. More may be added in subclasses
				
		Returns
			outStruct: plato_pylib UnitCell object containing the optimised geometry
		
		"""
		raise NotImplementedError

	def getStructsForEos(self, key):
		""" Get UnitCells for carrying out equation of states
		
		Args:
			key: structure key; hcp/bcc/fcc are standard options. More may be added in subclasses
				
		Returns
			outStructs: list of plato_pylib UnitCell object containing the geometries to use for calculating and energy vs volume curve
		
		"""
		raise NotImplementedError


	def getEosFitDict(self,key, eosModel="murnaghan"):
		""" Get dictionary containing eos fit data. 
		
		Args:
			key: structure key; hcp/bcc/fcc are standard options. More may be added in subclasses
			eosModel(Optional) = key to pass to ASE denoting eos model to use
		Returns
			outDict: dict containing info such as v0,b0,e0

		"""
		raise NotImplementedError


	def getSelfInterstitialPlaneWaveStruct(self, structType, interstitialType, relaxType, cellSize):
		""" Get structure of a self interstitial
		
		Args:
			structType: str, e.g. hcp/bcc/fcc
			interstitialType: str description of the type of interstitial (e.g. octahedral/tetrahedral)
			relaxType: str, unrelaxed/relaxed_constant_p/relaxed_constant_v
			cellSize: str, with cell dimensions (xdim_ydim_zdim)

		Returns
			outStruct: plato_pylib UnitCell structure
		
		Raises:
			Errors
		"""
		raise NotImplementedError


	def getVacancyPlaneWaveStruct(self, structType, relaxType, cellSize):
		""" Get structure for single vacancy
		
		Args:
			structType: str, e.g. hcp/bcc/fcc
			relaxType: str, unrelaxed/relaxed_constant_p/relaxed_constant_v
			cellSize: str, with cell dimensions (xdim_ydim_zdim)

		Returns
			outStruct: plato_pylib UnitCell structure
		
		Raises:
			Errors
		"""
		raise NotImplementedError

	def getPlaneWaveDosData(self, structKey):
		""" Get energy vs density of states numpy array for structure determined by structKey 

		Args:
			structKey: key that lets the code fetch the correct density of states

		Returns:
			dos: np array, 1st row is energies while second is density of states

		Raises:
			KeyError: If structKey not found
		"""
		raise NotImplementedError

class ShellAngMomMapper():
	""" Maps shell indices to and from angular momentum values for elements """
	def __init__(self, eleAngMomDict:dict):
		""" Create ShellAngMomMapper obj. 
		
		Args:
			eleAngMomDict: dict, keys are elements (case insensitive); values are itertables which list the angular momenta of shells in order
		Notes:
			This structure starts with numbering at zero, so the first shell is index 0

		"""
		self.eleAngMomDict = {k.lower():list(v) for k,v in eleAngMomDict.items()} 


	def getShellToAngMomDict(self, eleStr):
		""" Return a dictionary that maps shell indices to angular momenta values
		
		Args:
			eleStr: str(case-insensitive), should be symbol chemical element of interest
				
		Returns
			outDict: dict; keys are shell indices, values are angular momenta
		
		Raises:
			Errors
		"""
		currAngMoms = self.eleAngMomDict[eleStr.lower()] 
		currDict = {k:v for k,v in it.zip_longest( range(len(currAngMoms)), currAngMoms )}
		return currDict

	def getShellIdxToAngMom(self, eleStr, shellIdx):
		""" Return angular momentum value of given shell index of given element
		
		Args:
			eleStr: str(case-insensitive), should be symbol chemical element of interest
			shellIdx: Index of the shell of interest (the first shell has index of 0)
	
		Returns
			angMom: int, The angular momentum of the shell. 0 for s-type, 1 for p-type etc.
		
		Raises:
			Errors
		"""

		return self.eleAngMomDict[eleStr.lower()][shellIdx]

	def getShellIndicesForAngMom(self, eleStr, angMomVal):
		""" Return list of shell indices which have a given angular momentum value for given element
		
		Args:
			eleStr: str(case-insensitive), should be symbol chemical element of interest
			angMomVal: Angular momentum value of interest (0=s-type, 1=p-type etc)
				
		Returns
			shellIndices: list, all shell indices with given ang mom. Empty list if none found
		
		Raises:
			Errors
		"""

		angMomVals = self.eleAngMomDict[eleStr.lower()]
		outIndices = list()
		for idx,val in enumerate(angMomVals):
			if val == angMomVal:
				outIndices.append(idx)

		return outIndices


	def getSingleShellIndexFromAngMom(self, eleStr, angMomVal, nthVal=0):
		""" Return a single shell index for a given angular momentum. Takes the nth value found
		
		Args:
			eleStr: str(case-insensitive), should be symbol chemical element of interest
			angMomVal: Angular momentum value of interest (0=s-type, 1=p-type etc)
			nthVal: int, which occurance to take. 0= take first, 1=take second etc.
				
		Returns
			shellIdx: nth Shell index which has given angular momentum
		
		Raises:
			KeyError: If eleStr isnt found
			AngMomMissing (KeyError): If the angular momentum isnt found (including if its found in the list, but the nthVal of it isnt)
		"""

		allIndices = self.getShellIndicesForAngMom(eleStr, angMomVal)
		try:
			outIndex = allIndices[nthVal]
		except IndexError:
			raise AngMomMissing("angMomVal = {} not found at nthVal {}".format(angMomVal,nthVal))

		return outIndex
