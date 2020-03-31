


import os

import plato_pylib.parseOther.parse_qe_files as parseQE

from ..shared import config_vars as configVars
from ..shared import register_key_decorator as regKeyDeco

BASE_FOLDER = os.path.join( configVars.CASTEP_DB_PATH, "adsorbates" )


#Reference dictionaries
_SINGLE_ADSORBATE_GEOM_DICT   = dict()
_SINGLE_ADSORBATE_ENERGY_DICT = dict()

#Registration decorators
registerGeomDeco = regKeyDeco.RegisterKeyValDecorator(_SINGLE_ADSORBATE_GEOM_DICT, forceKeysToCase="lower")
registerTotalEnergyDeco = regKeyDeco.RegisterKeyValDecorator(_SINGLE_ADSORBATE_ENERGY_DICT, forceKeysToCase="lower") 

#These come from quantum espresso
def getPlaneWaveGasPhaseAdsorbateGeom(adsorbate):
	""" Get plane-wave optimised geometry for a gas phase molecule/atom
	
	Args:
		adsorbate: (str) Key for the required molecule/atom (e.g. h2o). Case insensitive
			 
	Returns
		 outCell: (plato_pylib UnitCell object) Optimised plane-wave geometry
 
	"""
	return _SINGLE_ADSORBATE_GEOM_DICT[adsorbate.lower()]()



def getPlaneWaveGasPhaseAdsorbateEnergy(adsorbate):
	""" Get plane-wave total energy for a gas phase molecule/atom
	
	Args:
		adsorbate: (str) Key for the required molecule/atom (e.g. h2o). Case insensitive

	Returns
		outEnergt: (float) The total energy for that molecule/atom
 
	"""
	return _SINGLE_ADSORBATE_ENERGY_DICT[adsorbate.lower()]()


@registerGeomDeco("h2")
def _getH2Geom():
	filePath = os.path.join(BASE_FOLDER,"h2","hydrogen.out")
	return parseQE.parseQuantumEspressoOutfile(filePath)["unitCell"] 

@registerTotalEnergyDeco("h2")
def _getH2Energy():
	filePath = os.path.join(BASE_FOLDER,"h2","hydrogen.out")
	return parseQE.parseQuantumEspressoOutfile(filePath)["energies"].electronicTotalE 


@registerGeomDeco("hydroxyl")
def _():
	filePath = os.path.join(BASE_FOLDER,"hydroxyl", "hydroxyl.out")
	return parseQE.parseQuantumEspressoOutfile(filePath)["unitCell"] 

@registerTotalEnergyDeco("hydroxyl")
def _():
	filePath = os.path.join(BASE_FOLDER,"hydroxyl", "hydroxyl.out")
	return parseQE.parseQuantumEspressoOutfile(filePath)["energies"].electronicTotalE 


@registerGeomDeco("h2o")
def _():
	filePath = os.path.join(BASE_FOLDER,"h2o", "water.out")
	return parseQE.parseQuantumEspressoOutfile(filePath)["unitCell"] 

@registerTotalEnergyDeco("h2o")
def _():
	filePath = os.path.join(BASE_FOLDER,"h2o", "water.out")
	return parseQE.parseQuantumEspressoOutfile(filePath)["energies"].electronicTotalE 



