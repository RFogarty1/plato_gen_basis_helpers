
import types

from . import base_flow as baseFlow


class SurfaceEnergyWorkflow(baseFlow.BaseLabelledWorkflow):
	""" Workflow for calculating a surface energy """

	def __init__(self, surfCalcObj, bulkCalcObj, surfaceAreaFromUnitCell, label=None):
		""" Initializer
		
		Args:
			surfCalcObj: (CalcMethod object) Object used for calculating energy of the surface slab
			bulkCalcObj: (CalcMethod object) Object used for calculating energy of the bulk system
			surfaceAreaFromUnitCell: (Function f(plato_pylib unitCell object)) Takes a unitCell object and returns a surface area. Recommended route is to use the SurfaceObject function with vacuum=0.0 and nLayers=1, but may need to be careful about making assumptions on how the lattice parameters/angles map to lattice vectors 
				
		"""
		self.eType = "electronicTotalE"
		self.bulkCalcObj = bulkCalcObj
		self.surfCalcObj = surfCalcObj
		self.surfaceAreaFromUnitCell = surfaceAreaFromUnitCell
		self.output = [types.SimpleNamespace()]

	def run(self):
		parsedFile = self.surfCalcObj.parsedFile
		numbSurfaceAtoms = parsedFile.numbAtoms
		uCell = parsedFile.unitCell
		surfaceArea = self.surfaceAreaFromUnitCell(uCell)
		self.output[0].surfaceEnergy = (numbSurfaceAtoms/(2*surfaceArea)) * (self._energyPerAtomSurface - self._energyPerAtomBulk)

	@property
	def namespaceAttrs(self):
		return [ ["surfaceEnergy"] ]	
	
	@property
	def _energyPerAtomBulk(self):
		return self._getEnergiesPerAtomFromCalcObj(self.bulkCalcObj)

	@property
	def _energyPerAtomSurface(self):
		return self._getEnergiesPerAtomFromCalcObj(self.surfCalcObj)

	def _getEnergiesPerAtomFromCalcObj(self, calcObj):
		parsedFile = calcObj.parsedFile
		totalEnergy = getattr(parsedFile.energies, self.eType)
		nAtoms = parsedFile.numbAtoms
		return totalEnergy/nAtoms


class SurfaceAreaFromUnitCellFunct():

	#Create surface from surfObjClass with with vacuum and [1,1,1] dimensions
	#Then its trivial to get the surface area from it	
	def __init__(self,surfObjClass):
		raise NotImplementedError("")


	def __call__(self, uCellObj):
		raise NotImplementedError("")
