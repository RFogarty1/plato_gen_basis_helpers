
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

