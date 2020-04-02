

import functools
import os

import plato_pylib.parseOther.parse_qe_files as parseQE

from ..shared import config_vars as configVars
from ..shared import register_key_decorator as regKeyDeco

BASE_FOLDER = os.path.join( configVars.CASTEP_DB_PATH, "mg", "plus_adsorbates" )


#Create reference dictionaries
_SURFACE_PLUS_SINGLE_ADSORBATE_GEOM_DICT = dict()
_SURFACE_PLUS_SINGLE_ADSORBATE_ENERGY_DICT = dict()

#Create decorators to register functions to the dictionaties (functions mean lazy-evaluation, which is usually what we want)
registerGeomDeco  = regKeyDeco.RegisterKeyValDecorator(_SURFACE_PLUS_SINGLE_ADSORBATE_GEOM_DICT)
registerTotalEnergyDeco = regKeyDeco.RegisterKeyValDecorator(_SURFACE_PLUS_SINGLE_ADSORBATE_ENERGY_DICT)

class RegistrationKey():
	
	def __init__(self, *,surfaceKey, adsorbateKey, siteKey, cellDims, surfRelax=True, adsorbRelax=True):
		self.attrNames = ["surfaceKey", "adsorbateKey", "siteKey", "cellDims", "surfRelax", "adsorbRelax"]
		self.surfaceKey = surfaceKey.lower()
		self.adsorbateKey = adsorbateKey.lower() if adsorbateKey is not None else adsorbateKey
		self.siteKey = siteKey.lower() if siteKey is not None else siteKey
		self.cellDims = tuple(cellDims) #Needs to be hashable (lists arent hashable)
		self.surfRelax = surfRelax
		self.adsorbRelax = adsorbRelax

	def __eq__(self,other):
		for currAttr in self.attrNames:
			if getattr(self,currAttr) != getattr(other,currAttr):
				return False
		return True

	def __hash__(self):
		return hash( tuple([getattr(self,x) for x in self.attrNames]) )

@functools.singledispatch
def getPlaneWaveSingleAdsorbateGeom(surface, adsorbate, site, cellDims, surfRelax=True, adsorbRelax=True):
	""" Returns plane-wave DFT optimised geometry for Mg + an asorbate molecule
	
	Args:
		surface: (str) The type of surface we're using; e.g. hcp0001
		adsorbate: (str) Denotes the type of molecule that adsorbs (e.g. "H" or "OH")
		site: (str) The adsorption site
		cellDims: (len 3 int iter) Cell dimensions in each direction (e.g. [2,2,4])
		surfRelax: (Bool) Whether the surface was relaxed
		adsorbRelax: (Bool) Whether the adsorbed molecule was relaxed
	
	NOTE: Options are combined into keys; available keys (from which options can possibly be deduced) can be found in _SURFACE_PLUS_SINGLE_ADSORBATE_GEOM_DICT.keys()			 
	Returns
		 outEnergy: (float) Total electronic energy of the full system
	"""
	regKey = RegistrationKey( surfaceKey=surface, adsorbateKey=adsorbate, siteKey=site, cellDims=cellDims, surfRelax=surfRelax, adsorbRelax=adsorbRelax ) 
	return _SURFACE_PLUS_SINGLE_ADSORBATE_GEOM_DICT[regKey]()

@getPlaneWaveSingleAdsorbateGeom.register(RegistrationKey)
def _(regKey):
    return _SURFACE_PLUS_SINGLE_ADSORBATE_GEOM_DICT[regKey]()


@functools.singledispatch
def getPlaneWaveAdsorptionSingleAdsorbateTotalEnergy(surface, adsorbate, site, cellDims, surfRelax=True, adsorbRelax=True):
	""" Returns plane-wave DFT energy of Mg + an adsorbate molecule. This is the TOTAL ENERGY of one structure, NOT an adsorption energy
	
	Args:
		surface: (str) The type of surface we're using; e.g. hcp0001
		adsorbate: (str) Denotes the type of molecule that adsorbs (e.g. "H" or "OH")
		site: (str) The adsorption site
		cellDims: (len 3 int iter) Cell dimensions in each direction (e.g. [2,2,4])
		surfRelax: (Bool) Whether the surface was relaxed
		adsorbRelax: (Bool) Whether the adsorbed molecule was relaxed
	
	NOTE: Options are combined into keys; available keys (from which options can possibly be deduced) can be found in _SURFACE_PLUS_SINGLE_ADSORBATE_ENERGY_DICT.keys()
		 
	Returns
		 outGeom: (plato_pylib UnitCell object) Total electronic energy of the full system
	"""
	regKey = RegistrationKey( surfaceKey=surface, adsorbateKey=adsorbate, siteKey=site, cellDims=cellDims, surfRelax=surfRelax, adsorbRelax=adsorbRelax ) 
	return _SURFACE_PLUS_SINGLE_ADSORBATE_ENERGY_DICT[regKey]()

@getPlaneWaveAdsorptionSingleAdsorbateTotalEnergy.register(RegistrationKey)
def _(regKey):
	return _SURFACE_PLUS_SINGLE_ADSORBATE_ENERGY_DICT[regKey]()

# Mg with hydrogen adsorbed
@registerGeomDeco(RegistrationKey(surfaceKey="hcp0001", adsorbateKey="H".lower(), siteKey="octahedral", cellDims=[2,2,4]))
def _getGeomHcp0001HOctahedral_224():
	outPath = os.path.join(BASE_FOLDER, "hydrogen", "hcp0001", "cell_2_2_4","Mg_H_ocathedral.out")
	outUCell = parseQE.parseQuantumEspressoOutfile(outPath)["unitCell"]
	return outUCell

@registerTotalEnergyDeco(RegistrationKey(surfaceKey="hcp0001", adsorbateKey="H".lower(), siteKey="octahedral", cellDims=[2,2,4]))
def _getEnergyHcp0001HOctahedral_224():
	outPath = os.path.join(BASE_FOLDER, "hydrogen", "hcp0001", "cell_2_2_4","Mg_H_ocathedral.out")
	outEnergy = parseQE.parseQuantumEspressoOutfile(outPath)["energies"].electronicTotalE
	return outEnergy


@registerGeomDeco(RegistrationKey(surfaceKey="hcp0001", adsorbateKey="H".lower(), siteKey="tetrahedral", cellDims=[2,2,4]))
def _getGeomHcp0001HOctahedral_224():
	outPath = os.path.join(BASE_FOLDER, "hydrogen", "hcp0001", "cell_2_2_4","Mg_H_tetrahedral.out")
	outUCell = parseQE.parseQuantumEspressoOutfile(outPath)["unitCell"]
	return outUCell

@registerTotalEnergyDeco(RegistrationKey(surfaceKey="hcp0001", adsorbateKey="H".lower(), siteKey="tetrahedral", cellDims=[2,2,4]))
def _getEnergyHcp0001HOctahedral_224():
	outPath = os.path.join(BASE_FOLDER, "hydrogen", "hcp0001", "cell_2_2_4","Mg_H_tetrahedral.out")
	outEnergy = parseQE.parseQuantumEspressoOutfile(outPath)["energies"].electronicTotalE
	return outEnergy



# Mg with OH adsorbed (0.25 and 1.0 coverage)
@registerGeomDeco(RegistrationKey(surfaceKey="hcp0001", adsorbateKey="OH".lower(), siteKey="fcc-hollow", cellDims=[2,2,4]))
def _():
	outPath = os.path.join(BASE_FOLDER, "hydroxyl", "hcp0001", "cell_2_2_4", "fcc_hollow", "Mg_hydroxyl_1.out")
	return parseQE.parseQuantumEspressoOutfile(outPath)["unitCell"]

@registerTotalEnergyDeco(RegistrationKey(surfaceKey="hcp0001", adsorbateKey="OH".lower(), siteKey="fcc-hollow", cellDims=[2,2,4]))
def _():
	outPath = os.path.join(BASE_FOLDER, "hydroxyl", "hcp0001", "cell_2_2_4", "fcc_hollow", "Mg_hydroxyl_1.out")
	return parseQE.parseQuantumEspressoOutfile(outPath)["energies"].electronicTotalE


@registerGeomDeco(RegistrationKey(surfaceKey="hcp0001-OH-0pt75-fcc-hollow", adsorbateKey="OH".lower(), siteKey="fcc-hollow", cellDims=[2,2,4]))
def _():
	outPath = os.path.join(BASE_FOLDER, "hydroxyl", "hcp0001", "cell_2_2_4", "fcc_hollow", "Mg_hydroxyl_4.out")
	return parseQE.parseQuantumEspressoOutfile(outPath)["unitCell"]


@registerTotalEnergyDeco(RegistrationKey(surfaceKey="hcp0001-OH-0pt75-fcc-hollow", adsorbateKey="OH".lower(), siteKey="fcc-hollow", cellDims=[2,2,4]))
def _():
	outPath = os.path.join(BASE_FOLDER, "hydroxyl", "hcp0001", "cell_2_2_4", "fcc_hollow", "Mg_hydroxyl_4.out")
	return parseQE.parseQuantumEspressoOutfile(outPath)["energies"].electronicTotalE



#Mg reference surfaces (not neccesarily blank)
@registerGeomDeco(RegistrationKey(surfaceKey="hcp0001", adsorbateKey=None, siteKey=None, cellDims=[2,2,4]))
def _():
	outPath = os.path.join(BASE_FOLDER, "plain_surfaces", "hcp0001", "cell_2_2_4", "Mg_pureslab.out")
	return parseQE.parseQuantumEspressoOutfile(outPath)["unitCell"]

@registerTotalEnergyDeco(RegistrationKey(surfaceKey="hcp0001", adsorbateKey=None, siteKey=None, cellDims=[2,2,4]))
def _():
	outPath = os.path.join(BASE_FOLDER, "plain_surfaces", "hcp0001", "cell_2_2_4", "Mg_pureslab.out")
	return parseQE.parseQuantumEspressoOutfile(outPath)["energies"].electronicTotalE


@registerGeomDeco(RegistrationKey(surfaceKey="hcp0001-OH-0pt75-fcc-hollow", adsorbateKey=None, siteKey=None, cellDims=[2,2,4]))
def _():
	outPath = os.path.join(BASE_FOLDER, "plain_surfaces", "hcp0001", "cell_2_2_4", "with_adsorbed", "hydroxyl_fcc_hollow","Mg_hydroxyl_3.out")
	return parseQE.parseQuantumEspressoOutfile(outPath)["unitCell"]

@registerTotalEnergyDeco(RegistrationKey(surfaceKey="hcp0001-OH-0pt75-fcc-hollow", adsorbateKey=None, siteKey=None, cellDims=[2,2,4]))
def _():
	outPath = os.path.join(BASE_FOLDER, "plain_surfaces", "hcp0001", "cell_2_2_4", "with_adsorbed", "hydroxyl_fcc_hollow","Mg_hydroxyl_3.out")
	return parseQE.parseQuantumEspressoOutfile(outPath)["energies"].electronicTotalE




