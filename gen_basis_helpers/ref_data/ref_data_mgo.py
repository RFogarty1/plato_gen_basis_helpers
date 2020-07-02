

import os

from ..castep import castep_creator as castepCreator
from ..shared import config_vars as configVars
from . import ref_elemental_objs as refEleObjs
from . import helpers_ref_data as helpers
import plato_pylib.shared.ucell_class as uCell


BASE_FOLDER = os.path.join( configVars.CASTEP_DB_PATH, "mg_o" )


def createMgOReferenceDataObj():
	return MgOReferenceDataObj()



class MgOReferenceDataObj(refEleObjs.RefElementalDataBase):

	def __init__(self):
		pass

	def getExptGeom(self,key="rocksalt"):
		return getExptRockSaltStruct()


	def getStructsForEos(self, structKey):
		return _getUCellsForBulkModCalcsStandard(structKey)

	def getEosFitDict(self,key,eos="murnaghan"):
		return getPlaneWaveEosFitDict(key,eos=eos)

	def getPlaneWaveSurfaceEnergy(self, structKey):
		return getPlaneWaveSurfEnergy(structKey)


def getExptRockSaltStruct():
	""" Data taken from 10.2183/pjab.55.43 (Sasaki 1979)
	"""
	return _getExptRockSaltStruct()

def _getExptRockSaltStruct():
	a = 4.217 #Angstrom, the only lattice parameter needed
	outUCell = _createMgORocksaltStructFromLattParam(a)
	outUCell.convAngToBohr()
	return outUCell

# Structures to use for bulk mod calculations (to be consistent)
def _getUCellsForBulkModCalcsStandard(structType:str):
	structTypeToFunct = {"rocksalt": _getRockSaltModUCellsStandard}
	return structTypeToFunct[structType.lower()]()


def _getRockSaltModUCellsStandard():
	#approximately 5% around the expt value; so volumes more like 15% around it
	aVals = [4.00, 4.04, 4.08, 4.12, 4.16, 4.20, 4.22,
	         4.24, 4.28, 4.32, 4.36, 4.40, 4.44]

	uCells = [_createMgORocksaltStructFromLattParam(a) for a in aVals]
	for x in uCells:
		x.convAngToBohr()
	return uCells


def _createMgORocksaltStructFromLattParam(a):
	""" Gets a plato_pylib unit cell structure for MgO when given just the lattice parameter
	
	Args:
		a: Lattice parameter, whatever units you want (output cell will have the same units)
			 
	Returns
		 outCell: (plato_pylib UnitCell) Contains the PRIMITIVE cell (2 atoms)
 
	"""

	cellVecs = [ [0.0  , 0.5*a, 0.5*a],
	             [0.5*a, 0.0  , 0.5*a],
	             [0.5*a, 0.5*a, 0.0  ] ]

	fractPos = [ [0.0, 0.0, 0.0,"Mg"], [0.5,0.5,0.5,"O"] ]
	return uCell.UnitCell.fromLattVects(cellVecs, fractCoords=fractPos)

#Plane wave geometry
def getPlaneWaveGeomParsedFileObject(structType="rocksalt"):
	structTypeToFunct = {"rocksalt":_getPlaneWaveRockSaltParsedFileObj}
	return structTypeToFunct[structType.lower()]()

def _getPlaneWaveRockSaltParsedFileObj():
	refPath = os.path.join(BASE_FOLDER,"opt_geom","rocksalt","geom_opt.castep")
	return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)

def getPlaneWaveGeom(structType="rocksalt"):
	return getPlaneWaveGeomParsedFileObject(structType).unitCell


#EoS data
def getPlaneWaveEosFitDict(structType:str, eos="murnaghan"):
	structTypeToFunct = {"rocksalt": _getPlaneWaveEosDictRocksalt}
	return structTypeToFunct[structType](eos)

def _getPlaneWaveEosDictRocksalt(eos):
	outFolder = os.path.join(BASE_FOLDER,"eos","rocksalt")
	return helpers.getEosFitDictFromEosCastepFolder(outFolder)

def getPlaneWaveSurfEnergy(structKey):
	outDict = {"rocksalt001": _getSurfaceEnergyRocksalt001}
	return outDict[structKey.lower()]()

def _getSurfaceEnergyRocksalt001():
	refBaseFolder = os.path.join(BASE_FOLDER, "surface_energies", "rocksalt_001")
	bulkFile = os.path.join(refBaseFolder, "bulk_calc", "conv_val_1100pt000.castep")
	surfFile = os.path.join(refBaseFolder, "surface_calc", "conv_val_1100pt000.castep")
	return helpers.getCastepRefRocksalt001SurfacEnergyFromSurfAndBulkFilePaths(surfFile,bulkFile)


def getPlaneWaveSurfaceParsedFileObject(surfType, nLayers=None, relaxType="unrelaxed"):
	""" Returns a ParsedFile object for castep calculation on a surface. Note the unrelaxed geometries are built from from the castep optimised bulk cell
	
	Args:
		surfType (str): The type of surface; examples are hcp0001 and hcp10m10 
		nLayers (int, optional): The number of surface layers used in the calculation. Default will vary for different surfTypes; but will represent converged structures
		relaxType (str, optional): String denoting the type of relaxation applied. Default="unrelaxed"; other options are "constant_volume" for now

	Returns
		parsedFile (ParsedFile object): This contains the geometry and total energy of the requested structure
	
	"""
	structTypeDefaultNLayers = {"rocksalt_001":4, "rocksalt_110":4}
	if nLayers is None:
		nLayers = structTypeDefaultNLayers[surfType]

	structTypeToFunct = { ("rocksalt_001", 4, "constant_volume"): _getRocksalt001PlaneWaveRelaxedParsedFile_4layers,
	                      ("rocksalt_110", 4, "constant_volume"): _getRocksalt110PlaneWaveRelaxedParsedFile_4layers }

	return structTypeToFunct[(surfType.lower(), nLayers, relaxType.lower())]()

def _getRocksalt001PlaneWaveRelaxedParsedFile_4layers():
	refPath = os.path.join(BASE_FOLDER,"surface_energies","rocksalt_001","relaxed_surface_calc","geom_opt.castep")
	return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)

def _getRocksalt110PlaneWaveRelaxedParsedFile_4layers():
	refPath = os.path.join(BASE_FOLDER,"surface_energies","rocksalt_110","relaxed_surface_calc","geom_opt.castep")
	return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)



