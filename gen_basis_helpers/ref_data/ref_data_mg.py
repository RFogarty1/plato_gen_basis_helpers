#!/usr/bin/python3

#Purpose of this code it to provide access to reference structures + similar to use for calculations

""" Provides access to reference structures and data for pure Mg """

import itertools as it
import os
import math
import sys
import types

from ..shared import config_vars as configVars
from ..castep import castep_creator as castepCreator
from ..job_utils import dos_helpers as dosHelp
from . import ref_elemental_objs as refEleObjs
from . import helpers_ref_data as helpers
import plato_fit_integrals.initialise.create_surf_energies_workflows as surfFlow
import gen_basis_helpers.shared.surfaces as surf
import plato_pylib.parseOther.parse_castep_files as parseCastep
import plato_pylib.plato.mod_plato_inp_files as modInp
import plato_pylib.plato.plato_paths as platoPaths
import plato_pylib.plato.parse_tbint_files as parseTbint
import plato_pylib.shared.ucell_class as UCell
import plato_pylib.utils.defects as defects
import plato_pylib.utils.supercell as supCell
import numpy as np


#tb1Model = os.path.join("Mg_bases_spd_att6","rc_7pt3","tb1_mcweda")
tb1Model = os.path.join("Mg_models","two_body_2019")
dft2Model = str(tb1Model)
dftModel = str(tb1Model)

BASE_FOLDER = os.path.join( configVars.CASTEP_DB_PATH, "mg" )


def createMgReferenceDataObj():
	tb1ModAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(tb1Model)
	dft2ModelAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(dft2Model)
	dftModelAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(dftModel,dtype="dft")
	modelHolder = platoPaths.PlatoModelFolders(tb1Path=tb1ModAbs, dft2Path=dft2ModelAbs, dftPath=dftModelAbs)
	return MgReferenceDataObj(modelHolder)


def createMgAngMomShellIndices():
	dft2ModelAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(dft2Model)
	dft2AdtPath = os.path.join(dft2ModelAbs, "Mg.adt")
	shellToAngMomDict = parseTbint.parseAdtFile(dft2AdtPath)["shellToAngMom"]
	
	angMomIndices = list()
	for key in range(len(shellToAngMomDict.keys())):
		angMomIndices.append( shellToAngMomDict[key] )

	return angMomIndices


#Overall object holding ref Data
class MgReferenceDataObj(refEleObjs.RefElementalDataBase):

	def __init__(self, modelHolder):
		self._modelHolder = modelHolder

	@property
	def modelFiles(self):
		return self._modelHolder

	def getExptGeom(self,key):
		keyToFunct = {"hcp": getExptStructAsUCell}
		return keyToFunct[key]()

	def getPlaneWaveGeom(self,key):
		return getPlaneWaveGeom(key)

	def getStructsForEos(self,key):
		return getUCellsForBulkModCalcs(key)
	
	def getEosFitDict(self,key,eos="murnaghan"):
		return getPlaneWaveEosFitDict(key,eos=eos)

	def getSelfInterstitialPlaneWaveStruct(self, structType, interstitialType, relaxType, cellSize):
		return getInterstitialPlaneWaveStruct(structType, interstitialType, relaxType, cellSize)

	def getSelfInterstitialPlaneWaveFormationEnergy(self, structType, interstitialType, relaxType, cellSize):
		return getInterstitialPlaneWaveFormationEnergy(structType, interstitialType, relaxType, cellSize)

	def getVacancyPlaneWaveStruct(self, structType, relaxType, cellSize):
		return getVacancyPlaneWaveStruct(structType, relaxType, cellSize)

	def getVacancyPlaneWaveEnergy(self, structType, relaxType, cellSize):
		return getVacancyPlaneWaveEnergy(structType, relaxType, cellSize)

	def getPlaneWaveDosData(self, structKey):
		return getDosPlaneWave(structKey)

	def getPlaneWaveDosGeom(self, structKey):
		return getDosPlaneWaveGeom(structKey)

	def getPlaneWaveAtomTotalEnergy(self, charge=0):
		if charge==0:
			return _getMgAtomPlaneWaveTotalEnergy()
		elif charge==2:
			return _getMg2PlusIonPlaneWaveTotalEnergy()
		else:
			raise ValueError("charge = {} is an invalid value".format(charge))

		inpFolder = os.path.join(BASE_FOLDER,"atom_calc")
		inpFiles = helpers.getCastepOutPathsForFolder(inpFolder)
		assert len(inpFiles)==1
		outEnergy = parseCastep.parseCastepOutfile(inpFiles[0])["energies"].electronicTotalE
		return outEnergy

	def getPlaneWaveDissocSepVsEnergy(self, relativeE=True, inBohr=True):
		if relativeE:
			zeroEnergy = 2*self.getPlaneWaveAtomTotalEnergy()
		else:
			zeroEnergy = 0.0
		sepsVsEnergies = getPlaneWaveDissocSepVsTotalE(inBohr=inBohr)
		for x in sepsVsEnergies:
			x[1] -= zeroEnergy
		return sepsVsEnergies

	def getPlaneWaveSurfaceEnergy(self, structKey):
		return getPlaneWaveSurfEnergy(structKey)


#Plane-wave reference atom energies
def _getMgAtomPlaneWaveTotalEnergy():
	inpFolder = os.path.join(BASE_FOLDER,"atom_calc")
	inpFiles = helpers.getCastepOutPathsForFolder(inpFolder)
	assert len(inpFiles)==1
	outEnergy = parseCastep.parseCastepOutfile(inpFiles[0])["energies"].electronicTotalE
	return outEnergy

def _getMg2PlusIonPlaneWaveTotalEnergy():
	inpFolder = os.path.join(BASE_FOLDER,"atom_calc","ion_calcs")
	inpFiles = helpers.getCastepOutPathsForFolder(inpFolder)
	assert len(inpFiles)==1
	outEnergy = parseCastep.parseCastepOutfile(inpFiles[0])["energies"].electronicTotalE
	return outEnergy


# Experimental Structure
def getExptStructAsUCell():
	''' From Walker 1959: DOI=10.1016/0001-6160(59)90090-2 Units of returned ucell are in bohr'''
	return _getMgExptHcpAsUCell()

def _getMgExptHcpAsUCell():
	''' From Walker 1959: DOI=10.1016/0001-6160(59)90090-2 Units of returned ucell are in bohr'''
	
	cOverA = 1.624
	aVal = 3.20922 #In Angstrom
	lattParams = [aVal,aVal,aVal*cOverA]

	outCell = helpers.getPerfectHcpMinimalUCell("Mg")
	outCell.setLattParams(lattParams)
	outCell.convAngToBohr()

	return outCell



def getPlaneWaveStableStackingFaultSlabParsedFileObject(surfType, faultType, nLayers=None,relaxType="constant_volume"):
	""" Get a parsed file object for a stable stacking fault geometry
		
	Args:
		surfType (str): hcp0001 only option at the moment
		faultType (str): "no_fault", "I2", "I1", "T2"
		nLayers (int): Number of layers used in the slab calculation
		relaxType (str): How the cell was relaxed; constant_volume is likely always what you want

	"""
	structTypeDefaultNLayers = {"hcp0001":10}
	if nLayers is None:
		nLayers = structTypeDefaultNLayers[surfType]

	if faultType=="no_fault":
		return getPlaneWaveSurfaceParsedFileObject(surfType,nLayers=nLayers, relaxType=relaxType)

	structTypeToFunct = { ("hcp0001", "i2", 10, "constant_volume"): _getHcp0001PlaneWaveI2FaultParsedFile_10layers_constantVol,
	                      ("hcp0001", "i1", 10, "constant_volume"): _getHcp0001PlaneWaveI1FaultParsedFile_10layers_constantVol,
	                      ("hcp0001", "t2", 10, "constant_volume"): _getHcp0001PlaneWaveT2FaultParsedFile_10layers_constantVol }

	return structTypeToFunct[ (surfType, faultType.lower(), nLayers, relaxType) ] ()

def _getHcp0001PlaneWaveI2FaultParsedFile_10layers_constantVol():
	refPath = os.path.join(BASE_FOLDER, "stacking_faults", "hcp_i2", "stable_geom", "nlayers_10_absvac_10_ang", "geom_opt.castep")
	return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)

def _getHcp0001PlaneWaveI1FaultParsedFile_10layers_constantVol():
	refPath = os.path.join(BASE_FOLDER, "stacking_faults", "hcp_i1", "stable_geom", "nlayers_10_absvac_10_ang", "geom_opt.castep")
	return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)

def _getHcp0001PlaneWaveT2FaultParsedFile_10layers_constantVol():
	refPath = os.path.join(BASE_FOLDER, "stacking_faults", "hcp_t2_fault", "stable_geom", "nlayers_10_absvac_10_ang", "geom_opt.castep")
	return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)

def getPlaneWaveSurfaceParsedFileObject(surfType, nLayers=None, relaxType="unrelaxed"):
	""" Returns a ParsedFile object for castep calculation on a surface. Note the unrelaxed geometries are built from from the castep optimised bulk cell
	
	Args:
		surfType (str): The type of surface; examples are hcp0001 and hcp10m10 
		nLayers (int, optional): The number of surface layers used in the calculation. Default will vary for different surfTypes; but will represent converged structures
		relaxType (str, optional): String denoting the type of relaxation applied. Default="unrelaxed"; other options are "constant_volume" for now

	Returns
		parsedFile (ParsedFile object): This contains the geometry and total energy of the requested structure
	
	"""
	structTypeDefaultNLayers = {"hcp0001":10, "hcp10m10":16, "hcp10m10_long_termination":16}
	if nLayers is None:
		nLayers = structTypeDefaultNLayers[surfType]

	structTypeToFunct = { ("hcp0001" , 10,"unrelaxed"): _getMgHcp0001PlaneWaveUnrelaxedParsedFile_10layers,
	                      ("hcp0001" , 10,"constant_volume"  ): _getMgHcp0001PlaneWaveRelaxedParsedFile_10layers,
	                      ("hcp10m10", 16,"unrelaxed"): _getMgHcp10m10PlaneWaveUnrelaxedParsedFile_16layers,
	                      ("hcp10m10", 16, "constant_volume" ): _getMgHcp10m10PlaneWaveRelaxedParsedFile_16layers,
	                      ("hcp10m10_long_termination", 16, "constant_volume"): _getMgHcp10m10_longTerminationPlaneWaveRelaxedParsedFile_16layers
	                     }

	return structTypeToFunct[(surfType,nLayers,relaxType)]()

def _getMgHcp0001PlaneWaveUnrelaxedParsedFile_10layers():
    refPath = os.path.join(BASE_FOLDER,"surface_energies", "hcp0001", "castep_geom", "unrelaxed", "surface_n10_vac_18pt90.castep")
    return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)

def _getMgHcp0001PlaneWaveRelaxedParsedFile_10layers():
    refPath = os.path.join(BASE_FOLDER,"surface_energies", "hcp0001", "castep_geom", "relaxed", "nlayers_10_absvac_10_ang", "geom_opt.castep")
    return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)

def _getMgHcp10m10PlaneWaveUnrelaxedParsedFile_16layers():
	refPath = os.path.join(BASE_FOLDER,"surface_energies", "hcp10m10", "castep_geom", "unrelaxed", "surface_n16_vac_18pt90.castep")
	return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)

def _getMgHcp10m10PlaneWaveRelaxedParsedFile_16layers():
	refPath = os.path.join(BASE_FOLDER,"surface_energies", "hcp10m10", "castep_geom", "relaxed", "nlayers_16_absvac_10_ang", "geom_opt.castep")
	return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)

def _getMgHcp10m10_longTerminationPlaneWaveRelaxedParsedFile_16layers():
	refPath = os.path.join(BASE_FOLDER,"surface_energies", "hcp10m10_long_termination", "castep_geom", "relaxed", "nlayers_16_absvac_10_ang", "geom_opt.castep")
	return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)

#INTERFACE FUNCTION
def getPlaneWaveGeom(structType:str):
	structTypeToFunct = {"hcp": _getMgPlaneWaveHcpGeomAsUCell,
	                     "bcc": _getMgPlaneWaveBCCGeomAsUCell,
	                     "fcc": _getMgPlaneWaveFCCGeomAsUCell,
	                     "atom":_getMgAtomGeomAsUCell}
	return structTypeToFunct[structType.lower()]()

def getPlaneWaveGeomParsedFileObject(structType:str):
	structTypeToFunct = {"hcp": _getMgPlaneWaveHcpGeomParsedFile}
	return structTypeToFunct[structType.lower()]()


def _getMgPlaneWaveHcpGeomParsedFile():
	refPath = os.path.join(BASE_FOLDER,"opt_geoms","hcp","Mg_hcp_SPE_otf_10el_usp_PP_6pt06.castep")
	return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)


#TODO: These functions should be removed and a wrapper around getPlaneWaveGeomParsedFileObject should be used to fetch geometries
# Optimised Plane-wave structures (PBE)
def _getMgPlaneWaveHcpGeomAsUCell():
	refPath = os.path.join(BASE_FOLDER,"opt_geoms","hcp","Mg_hcp_SPE_otf_10el_usp_PP_6pt06.geom")
	initUCell = parseCastep.parseCastepGeomFile(refPath)[-1]["unitCell"]
	return initUCell

def _getMgPlaneWaveBCCGeomAsUCell():
	refPath = os.path.join(BASE_FOLDER,"opt_geoms","bcc", "Mg_bcc_opt_otf_10el_usp_PP.geom")
	initUCell = parseCastep.parseCastepGeomFile(refPath)[-1]["unitCell"]
	return initUCell

def _getMgPlaneWaveFCCGeomAsUCell():
	refPath = os.path.join(BASE_FOLDER,"opt_geoms", "fcc", "Mg_fcc_opt_otf_10el_usp_PP.geom")
	initUCell = parseCastep.parseCastepGeomFile(refPath)[-1]["unitCell"]
	return initUCell

#Single atom geom used to get plane-wave atomic energy
def _getMgAtomGeomAsUCell():
	inpFolder = os.path.join(BASE_FOLDER,"atom_calc")
	inpFiles = helpers.getCastepOutPathsForFolder(inpFolder)
	assert len(inpFiles)==1
	outCell = parseCastep.parseCastepOutfile(inpFiles[0])["unitCell"]
	outCell.convAngToBohr()
	return outCell



# Structures to use for bulk mod calculations (to be consistent)
def getUCellsForBulkModCalcs(structType:str):
	structTypeToFunct = {"hcp": _getHCPBulkModUCells,
	                     "bcc": _getBCCBulkModUCells,
	                     "fcc": _getFCCBulkModUCells}
	return structTypeToFunct[structType.lower()]()

def _getHCPBulkModUCells():
	refFolder = os.path.join(BASE_FOLDER,"eos","hcp")
	return helpers.getUCellsFromCastepBulkModFolder(refFolder)

def _getBCCBulkModUCells():
	refFolder = os.path.join(BASE_FOLDER,"eos","bcc")
	return helpers.getUCellsFromCastepBulkModFolder(refFolder)

def _getFCCBulkModUCells():
	refFolder = os.path.join(BASE_FOLDER,"eos","fcc")
	return helpers.getUCellsFromCastepBulkModFolder(refFolder)


#Random structures:
def getHcpComprStructAsUCell():
	refPath = os.path.join(BASE_FOLDER,"eos", "hcp", "Mg_hcp_SPE_otf_10el_usp_PP_5pt81.castep")
	initUCell = parseCastep.parseCastepOutfile(refPath)["unitCell"]
	initUCell.convAngToBohr()
	return initUCell


def getBccComprStructAsUCell():
	refPath = os.path.join(BASE_FOLDER,"eos", "bcc", "Mg_bcc_opt_otf_10el_usp_PP_10_5pt592.castep")
	initUCell = parseCastep.parseCastepOutfile(refPath)["unitCell"]
	initUCell.convAngToBohr()
	return initUCell


#Plane wave bulk mod fits using ASE (Actually does the fit while function is called)
def getPlaneWaveEosFitDict(structType:str, eos="murnaghan"):
	structTypeToFunct = {"hcp": _getPlaneWaveEosDictHCP,
	                     "fcc": _getPlaneWaveEosDictFCC,
	                     "bcc": _getPlaneWaveEosDictBCC}
	return structTypeToFunct[structType](eos)

def _getPlaneWaveEosDictHCP(eos):
	outFolder = os.path.join(BASE_FOLDER,"eos","hcp")
	return helpers.getEosFitDictFromEosCastepFolder(outFolder)

def _getPlaneWaveEosDictFCC(eos):
	outFolder = os.path.join(BASE_FOLDER,"eos","fcc")
	return helpers.getEosFitDictFromEosCastepFolder(outFolder)

def _getPlaneWaveEosDictBCC(eos):
	outFolder = os.path.join(BASE_FOLDER,"eos","bcc")
	return helpers.getEosFitDictFromEosCastepFolder(outFolder)



#Density of states for planeWave Calcs
def getDosPlaneWave(structType:str):
	structTypeToFunct = {"hcpExpt".lower(): _getHCPDosPlaneWaveExptGeom,
	                     "fcc": _getFCCDosPlaneWave,
	                     "bcc": _getBCCDosPlaneWave,
	                     "hcpCompr".lower(): _getHcpDosPlaneWaveComprGeomA,
	                     "bccCompr".lower():_getBccDosPlaneWaveComprGeomA}

	return structTypeToFunct[structType.lower()]()

def _getHCPDosPlaneWaveExptGeom():
	outFile = os.path.join(BASE_FOLDER, "dos", "hcp_eqm", "Mg_hcp.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])

def _getFCCDosPlaneWave():
	outFile = os.path.join(BASE_FOLDER, "dos", "fcc_eqm", "Mg_fcc.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])

def _getBCCDosPlaneWave():
	outFile = os.path.join(BASE_FOLDER, "dos", "bcc_eqm", "Mg_bcc.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])


def _getHcpDosPlaneWaveComprGeomA():
	outFile = os.path.join(BASE_FOLDER, "dos", "hcp_compr", "Mg_hcp.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])

def _getBccDosPlaneWaveComprGeomA():
	outFile = os.path.join(BASE_FOLDER, "dos", "bcc_compr", "Mg_bcc.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])



#Density og states as geom
def getDosPlaneWaveGeom(structType:str):

	structTypeToFunct = {"hcpExpt".lower(): getExptStructAsUCell,
	                     "hcpCompr".lower(): _getHcpComprGeomForDos,
	                     "bccCompr".lower(): _getBccComprGeomForDos,
	                     "bcc".lower(): _getBccEqmGeomFosDos}

	return structTypeToFunct[structType.lower()]()


def _getHcpComprGeomForDos():
	outFile = os.path.join(BASE_FOLDER, "dos", "hcp_compr", "Mg_hcp_SPE_otf_10el_usp_PP_5pt81.castep")
	return helpers.getUCellInBohrFromCastepOutFile(outFile)

def _getBccComprGeomForDos():
	outFile = os.path.join(BASE_FOLDER, "dos", "bcc_compr", "Mg_bcc_opt_otf_10el_usp_PP_10_5pt592.castep")
	return helpers.getUCellInBohrFromCastepOutFile(outFile)

def _getBccEqmGeomFosDos():
	outFile = os.path.join(BASE_FOLDER, "dos", "bcc_eqm", "Mg_bcc_opt_otf_10el_usp_PP_10_5pt842.castep")
	return helpers.getUCellInBohrFromCastepOutFile(outFile)


#Dissociation energy curves
def getPlaneWaveDissocSepVsTotalE(inBohr=True):
	outFolder = os.path.join(BASE_FOLDER,"dissoc_curve")
	outFiles = helpers.getCastepOutPathsForFolder(outFolder)
	outEnergies = [parseCastep.parseCastepOutfile(x)["energies"].electronicTotalE for x in outFiles]
	outSeps = [helpers.getDimerSepFromCastepOutFile(x,inBohr=True) for x in outFiles]
	outList = list()
	for sep,e in it.zip_longest(outSeps,outEnergies):
		outList.append( [sep,e] ) 

	return outList


def getInterstitialPlaneWaveParsedFile(structType, interstitialType, relaxType, cellSize):
	paramsToOutFunct = {("hcp","no_inter", "plane_wave_geom", "3_3_2"): _getHcpPlaneWaveNoInter332,
	                    ("hcp","octahedral","relaxed_constant_pressure", "3_3_2"): _getHcpPlaneWaveParsedFile_interOctaRelaxedConstPressure332,
	                    ("hcp","split","relaxed_constant_pressure","3_3_2"): _getHcpPlaneWaveParsedFile_interSplitRelaxedConstPressure332,
	                    ("hcp","octahedral","unrelaxed","3_3_2"): _getHcpPlaneWaveParsedFile_interOctaUnrelaxed332,
	                    ("hcp","tetrahedral_old","unrelaxed","3_3_2"): _getHcpPlaneWaveParsedFile_interTetraOldUnrelaxed332
}
	return paramsToOutFunct[(structType,interstitialType,relaxType,cellSize)]()

def _getHcpPlaneWaveNoInter332():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "no_inter", "no_inter.castep")
	parsedFile = castepCreator.getParsedFileObjFromCastepOutputFile(refFile)
	return parsedFile

def _getHcpPlaneWaveParsedFile_interOctaUnrelaxed332():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "unrelaxed", "octa", "octa_unrelaxed.castep")
	parsedFile = castepCreator.getParsedFileObjFromCastepOutputFile(refFile)
	return parsedFile

def _getHcpPlaneWaveParsedFile_interTetraOldUnrelaxed332():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "unrelaxed", "hcp_tetra_inter.castep")
	parsedFile = castepCreator.getParsedFileObjFromCastepOutputFile(refFile)
	return parsedFile

#Structures to use for intersitials
def getInterstitialPlaneWaveStruct(structType:"str, e.g. hcp", interstitialType:"str, octahedral or tetrahedral",
                                   relaxType:"str, unrelaxed or relaxed", cellSize:"Str with dims, e.g 3_3_2"):

	return getInterstitialPlaneWaveParsedFile(structType, interstitialType, relaxType, cellSize).unitCell


def _getHcpPlaneWaveParsedFile_interOctaRelaxedConstPressure332():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "relaxed", "constant_p", "octa","mg_octa_self_inter.castep")
	parsedFile = castepCreator.getParsedFileObjFromCastepOutputFile(refFile)
	return parsedFile

def _getHcpPlaneWaveParsedFile_interSplitRelaxedConstPressure332():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "relaxed", "constant_p", "tetra", "hcp_tetra_constant_p_inter_spe.castep") #file called tetra, but actually split
	parsedFile = castepCreator.getParsedFileObjFromCastepOutputFile(refFile)
	return parsedFile


def getInterstitialPlaneWaveFormationEnergy(structType, interstitialType, relaxType, cellSize):
	parsedInter = getInterstitialPlaneWaveParsedFile(structType, interstitialType, relaxType, cellSize)
	parsedNoInter = getInterstitialPlaneWaveParsedFile(structType, "no_inter", "plane_wave_geom", cellSize)
	return helpers.getDefectEnergyFromCastepNoDefectAndDefectParsedFileObjs(parsedNoInter, parsedInter)

def getVacancyPlaneWaveStruct(structType, relaxType, cellSize):
	baseUCell = getPlaneWaveGeom(structType)
	cellDims = [int(x) for x in cellSize.split("_")]

	if relaxType == "unrelaxed":
		vacCell = supCell.superCellFromUCell(baseUCell, cellDims)
		defects.makeVacancyUnitCell(vacCell)
		return vacCell
	elif relaxType == "novac":
		return supCell.superCellFromUCell(baseUCell, cellDims)
	elif relaxType == "relaxed_constant_pressure":
		return _getVacancyHcpPlaneWaveStruct_relaxedConstantPressure()
	else:
		raise NotImplementedError("Only unrelaxed vancancies are currently implemented")


def _getVacancyHcpPlaneWaveStruct_relaxedConstantPressure():
	vacFilePath = os.path.join(BASE_FOLDER, "vacancy", "relaxed", "constant_p", "vacancy_opt_500ev_8_8_6.castep")
	return helpers.getUCellInBohrFromCastepOutFile(vacFilePath)

def getVacancyPlaneWaveEnergy(structType, relaxType, cellSize):
	paramsToEnergyDict = {("hcp", "unrelaxed", "3_3_2"): _getVacancyEnergyHcpUnrelaxed332,
	                      ("hcp", "relaxed_constant_pressure","3_3_2"): _getVacancyEnergyHcpRelaxedConstantPressure332}
	return paramsToEnergyDict[(structType, relaxType, cellSize)]()

def _getVacancyEnergyHcpUnrelaxed332():
	baseFolder = os.path.join(BASE_FOLDER, "vacancy", "unrelaxed", "3_3_2")
	vacFilePath = os.path.join(baseFolder, "vacancy_500ev_8_8_6.castep" )
	bulkFilePath = os.path.join(baseFolder, "novac_500ev_8_8_6.castep")
	outEnergy = helpers.getDefectEnergyFromCastepNoDefectAndDefectFiles(bulkFilePath, vacFilePath)
	return outEnergy


def _getVacancyEnergyHcpRelaxedConstantPressure332():
	vacFilePath = os.path.join(BASE_FOLDER, "vacancy", "relaxed", "constant_p", "vacancy_opt_500ev_8_8_6.castep")
	bulkFilePath = os.path.join(BASE_FOLDER, "interstitial", "no_inter", "no_inter.castep")
	return helpers.getDefectEnergyFromCastepNoDefectAndDefectFiles(bulkFilePath, vacFilePath)

def getPlaneWaveSurfEnergy(structKey):
	outDict = {"hcp0001": _getSurfaceEnergyHcp0001}

	return outDict[structKey.lower()]()


def _getSurfaceEnergyHcp0001():
	refBaseFolder = os.path.join(BASE_FOLDER,"surface_energies", "hcp0001")
	bulkModFile = os.path.join(refBaseFolder, "Mg_hcp_SPE_otf_10el_usp_PP_6pt06.castep")
	surfFile = os.path.join(refBaseFolder, "Mg24_k1.castep")

	return helpers.getCastepRefHcp0001SurfaceEnergyFromSurfAndBulkFilePaths(surfFile,bulkModFile)


#struct is hcp0001
def getPWUnrelaxedStackFaultDispsVsStackFaultValsForSlab(structType, faultType, direction, nLayers=None):
	""" Gets displacements (units of lattice param) against unrelaxed stacking fault energies (eV a_0^{-2}) for a set of plane-wave geometries
	
	Args:
		structType (str): Surface type, e.g hcp0001
		faultType (str): Label for the type of fault, e.g. i2,i1,t2
		direction (str): Displacement direction, e.g. 10m10
		nLayers (opt,int): Number of layers in the cell. Sensible default chosen if not passed
			 
	Returns
		data: SimpleNamespace containing displacement values and stacking fault values for the chosen structure
 
	"""
	structTypeDefaultNLayers = {"hcp0001":10}
	if nLayers is None:
		nLayers = structTypeDefaultNLayers[(structType)]

	structTypeToFunct = { ("hcp0001", "i2", "10m10", 10): _getHcp0001i2_10m10_dispValsVsStackFaultEnergies_10Layers,
	                      ("hcp0001", "i1", "10m10", 10): _getHcp0001i1_10m10_dispValsVsStackFaultEnergies_10Layers,
	                      ("hcp0001", "t2", "10m10", 10): _getHcp0001t2_10m10_dispValsVsStackFaultEnergies_10Layers }

	return structTypeToFunct[(structType, faultType, direction, nLayers)]() 


def _getHcp0001i2_10m10_dispValsVsStackFaultEnergies_10Layers():
	dataFolder = os.path.join(BASE_FOLDER, "stacking_faults", "hcp_i2", "unrelaxed_10m10")
	dispVals, stackFaultVals = helpers.getDispsAndStackFaultEnergiesFromCastepFilesInFolder(dataFolder)
	outStruct = types.SimpleNamespace( dispVals=dispVals, stackFaultVals=stackFaultVals )
	return outStruct

def _getHcp0001i1_10m10_dispValsVsStackFaultEnergies_10Layers():
	dataFolder = os.path.join(BASE_FOLDER, "stacking_faults", "hcp_i1", "unrelaxed_10m10")
	dispVals, stackFaultVals = helpers.getDispsAndStackFaultEnergiesFromCastepFilesInFolder(dataFolder)
	outStruct = types.SimpleNamespace( dispVals=dispVals, stackFaultVals=stackFaultVals )
	return outStruct

def _getHcp0001t2_10m10_dispValsVsStackFaultEnergies_10Layers():
	dataFolder = os.path.join(BASE_FOLDER, "stacking_faults", "hcp_t2_fault", "unrelaxed_10m10")
	dispVals, stackFaultVals = helpers.getDispsAndStackFaultEnergiesFromCastepFilesInFolder(dataFolder)
	outStruct = types.SimpleNamespace( dispVals=dispVals, stackFaultVals=stackFaultVals )
	return outStruct

