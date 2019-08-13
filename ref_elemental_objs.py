
from plato_pylib.plato.plato_paths import PlatoPathDescriptor, PlatoModelFolders


class RefElementalDataBase():
	""" Class for holding reference data for a pure-elemental system """


	@property
	def modelFiles(self):
		""" plato pylib PlatoModelFolders object, handles paths to model datafolders """
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
		return NotImplementedError


