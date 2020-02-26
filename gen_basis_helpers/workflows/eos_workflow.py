
import types

import plato_pylib.utils.fit_eos as fitEos
from . import base_flow as baseFlow



class EosWorkflow(baseFlow.BaseLabelledWorkflow):

	def __init__(self, calcObjs, fitFunction, label, eType="electronicTotalE"):
		""" Initializer for EosWorkflow; used to calculate energy-volume curves and associated data
		
		Args:
			calcObjs: (iter of CalcMethod objects) 
			fitFunction: f(volumes, energies)->outDict. volumes and energies are both iterables. For the outDict keys required you'll have to look inside this class for now (sorry). Standard function to use can be found in this module though
			label: (StandardLabel object) Used to differentiate between workflows
			eType: (Str) Type of energy we want; Corresponds to attributes on plato_pylib Energies object. Default is total electronic energy
	
		Raises:
			Errors
		"""

		self._calcObjs = calcObjs #invisible to the composite workflows
		self._fitFunction = fitFunction #Invisible to the composite workflows
		self._label = label
		self._eType = eType

	@property
	def preRunShellComms(self):
		outComms = list()
		for x in self._calcObjs:
			outComms.append( x.runComm ) #Each command is a single string
		return outComms
	
	@property
	def label(self):
		return [self._label]

	@property
	def namespaceAttrs(self):
		return [["data"]]

	def run(self):
		volumes, energies = list(), list()
		for x in self._calcObjs:
			currVol = x.parsedFile.unitCell.volume
			currEnergy = getattr(x.parsedFile.energies, self._eType)
			ePerAtom = currEnergy / x.parsedFile.numbAtoms
		self.output = [types.SimpleNamespace( data=self._fitFunction(volumes,energies) )] #List to allow composite pattern use



class StandardEosFitFunction():
	""" Callable class representing an eos fit. The interface of the call function is:
	
	Args:
		volumes: (float iter) Volume per atom values for structures used to calculate energy-volume curves. Units = bohr^3
		energies: (float iter) Energy per atom values for structures used to calculated energy-volume curves. Units = eV
	Returns
		fitDict: A dictionary containing various parameters relating to the EoS fit. GPa used for bulk-modulus, bohr^3 per atom for bolume
	
	"""

	def __init__(self, eosStr="murnaghan", maxFev=10000):
		""" Initializer for setting up a standard eos-fit function
		
		Args:
			eoStr: (str) Representation of which equation of state to use. This is (at time of writing) passed directly to the Atomic Simulation Environment implementation
			maxFev: (int) Maximum number of function evaluations to use when fitting the EoS.
				
		"""
		self.eosStr = eosStr
		self.maxFev = maxFev


	def __call__(self, volumes, energies):
		return fitEos.getBulkModFromVolsAndEnergiesBohrAndEvUnits(volumes, energies, eosModel=self.eosStr, maxFev=self.maxFev)

