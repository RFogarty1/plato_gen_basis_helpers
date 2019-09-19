#!/usr/bin/python3

#Purpose of this code it to provide access to reference structures + similar to use for calculations

""" Provides access to reference structures and data for pure Mg """

import os
import math
import sys

from ..job_utils import dos_helpers as dosHelp
from . import ref_elemental_objs as refEleObjs
from . import helpers_ref_data as helpers
import plato_pylib.parseOther.parse_castep_files as parseCastep
import plato_pylib.plato.mod_plato_inp_files as modInp
import plato_pylib.plato.plato_paths as platoPaths
import plato_pylib.plato.parse_tbint_files as parseTbint
import plato_pylib.shared.ucell_class as UCell
import plato_pylib.utils.defects as defects
import plato_pylib.utils.supercell as supCell
import numpy as np


tb1Model = os.path.join("Mg_bases_spd_att6","rc_7pt3","tb1_mcweda")
dft2Model = str(tb1Model)
dftModel = str(tb1Model)

BASE_FOLDER = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/testing_models/castep_database/mg"

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

	def getVacancyPlaneWaveStruct(self, structType, relaxType, cellSize):
		return getVacancyPlaneWaveStruct(structType, relaxType, cellSize)



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


#INTERFACE FUNCTION
def getPlaneWaveGeom(structType:str):
	structTypeToFunct = {"hcp":_getMgPlaneWaveHcpGeomAsUCell,
	                     "bcc": _getMgPlaneWaveBCCGeomAsUCell,
	                     "fcc": _getMgPlaneWaveFCCGeomAsUCell}
	return structTypeToFunct[structType.lower()]()


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
	refPath = os.path.join(BASE_FOLDER,"eos", "Mg_bcc_opt_otf_10el_usp_PP_10_5pt592.castep")
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



#Structures to use for intersitials

def getInterstitialPlaneWaveStruct(structType:"str, e.g. hcp", interstitialType:"str, octahedral or tetrahedral",
                                   relaxType:"str, unrelaxed or relaxed", cellSize:"Str with dims, e.g 3_3_2"):

	paramsToStructDict = {("hcp","tetrahedral","unrelaxed","3_3_2"):_getHcpPlaneWaveStruct_interTetraUnrelaxed332(),
	                      ("hcp","octahedral" ,"unrelaxed","3_3_2"):_getHcpPlaneWaveStruct_interOctaUnrelaxed332(),
	                      ("hcp","octahedral", "relaxed_constant_pressure","3_3_2"):_getHcpPlaneWaveStruct_interOctaRelaxedConstPressure332()}

	return paramsToStructDict[(structType,interstitialType,relaxType,cellSize)]



def _getHcpPlaneWaveStruct_interTetraUnrelaxed332():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "unrelaxed", "hcp_tetra_inter.castep")
	pasedUCell = parseCastep.parseCastepOutFile(refRef)["unitCell"]
	parsedUCell.convAngToBohr()
	return parsedUCell


def _getHcpPlaneWaveStruct_interOctaUnrelaxed332():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "unrelaxed", "hcp_octa_inter.castep")
	parsedUCell = parseCastep.parseCastepOutFile(refRef)["unitCell"]
	parsedUCell.convAngToBohr()
	return parsedUCell


#NOTE: This geometry is wrong, acidentally obtained with the wrong k-point mesh being used in the final optimisation, hence this will give a bug until i replace the file
def _getHcpPlaneWaveStruct_interOctaRelaxedConstPressure332():
	refFile = os.path.join(BASE_FOLDER,"interstitial","relaxed","constant_p","hcp_octa_inter.geom")
	outUCell = parseCastep.parseCastepGeomFile(refFile)[-1]["unitCell"]
	return outUCell



def getVacancyPlaneWaveStruct(structType, relaxType, cellSize):
	baseUCell = getPlaneWaveGeom(structType)
	cellDims = [int(x) for x in cellSize.split("_")]

	if relaxType == "unrelaxed":
		vacCell = supCell.superCellFromUCell(baseUCell, cellDims)
		defects.makeVacancyUnitCell(vacCell)
		return vacCell
	else:
		raise NotImplementedError("Only unrelaxed vancancies are currently implemented")


