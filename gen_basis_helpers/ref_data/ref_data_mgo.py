

import os
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
		return _getExptRockSaltStruct()


	def getStructsForEos(self, structKey):
		return _getUCellsForBulkModCalcsStandard(structKey)

	def getEosFitDict(self,key,eos="murnaghan"):
		return getPlaneWaveEosFitDict(key,eos=eos)

	def getPlaneWaveSurfaceEnergy(self, structKey):
		return getPlaneWaveSurfEnergy(structKey)


def _getExptRockSaltStruct():
	""" Data taken from 10.2183/pjab.55.43 (Sasaki 1979)
	"""
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



