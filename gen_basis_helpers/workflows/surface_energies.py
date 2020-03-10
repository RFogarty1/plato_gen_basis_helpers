
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

		self._writeInpFiles()

	def _writeInpFiles(self):
		for x in [self.bulkCalcObj,self.surfCalcObj]:
			x.writeFile()

	def run(self):
		parsedFile = self.surfCalcObj.parsedFile
		numbSurfaceAtoms = parsedFile.numbAtoms
		uCell = parsedFile.unitCell
		surfaceArea = self.surfaceAreaFromUnitCell(uCell)
		self.output[0].surfaceEnergy = (numbSurfaceAtoms/(2*surfaceArea)) * (self._energyPerAtomSurface - self._energyPerAtomBulk)
		self.output[0].surfEPerAtom = self._energyPerAtomSurface
		self.output[0].bulkEPerAtom = self._energyPerAtomBulk

	@property
	def preRunShellComms(self):
		return [self.bulkCalcObj.runComm,self.surfCalcObj.runComm]


	@property
	def namespaceAttrs(self):
		return [ ["surfaceEnergy","surfEPerAtom","bulkEPerAtom"] ]	
	
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
		""" Initializer
		
		Args:
			surfObjClass: (BaseSurface class, NOT INSTANCE) This has the ability to map a unitCell into a surface object. This will likely need to include any rotations required such that the surface is along the z-dir (though at time of writing only the hcp0001 class is implemented, and this doesnt require any rotating)
		
		Note the callable function takes a plato_pylib unitCell object as its sole argument.

		I'm also assuming the surface object initialiser has the signature (self,bulkUCell,nLaters,lenVac)
		
		"""
		self.surfClass = surfObjClass


	def __call__(self, uCellObj):
		nLayers, lenVac = 1,0
		surfaceObj = self.surfClass(uCellObj, nLayers, lenVac)
		return surfaceObj.surfaceArea

