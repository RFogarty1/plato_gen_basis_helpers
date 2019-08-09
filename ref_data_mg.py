#!/usr/bin/python3

#Purpose of this code it to provide access to reference structures + similar to use for calculations

""" Provides access to reference structures and data for pure Mg """

import os
import math
import sys

import dos_helpers as dosHelp
import ref_elemental_objs as refEleObjs
import plato_pylib.shared.ucell_class as UCell
import plato_pylib.parseOther.parse_castep_files as parseCastep
import plato_pylib.plato.mod_plato_inp_files as modInp
import plato_pylib.plato.plato_paths as platoPaths
import numpy as np


sys.path.append("/media/ssd1/rf614/Work/usr_scripts/coding/Plato_Analysis_Lib_Functions")
import fit_bulk_mod as fitBMod

tb1Model = os.path.join("Mg_bases_spd_att6","rc_7pt3","tb1_mcweda")
dft2Model = str(tb1Model)
dftModel = str(tb1Model)



def createMgReferenceDataObj():
	basePath = "/media/ssd1/rf614/Work/Plato"
	tb1ModAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(tb1Model)
	dft2ModelAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(dft2Model)
	dftModelAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(dftModel,dtype="dft")
	modelHolder = platoPaths.PlatoModelFolders(tb1Path=tb1ModAbs, dft2Path=dft2ModelAbs, dftPath=dftModelAbs)
	return MgReferenceDataObj(modelHolder)

#Overall object holding ref Data
class MgReferenceDataObj(refEleObjs.RefElementalDataBase):

	def __init__(self, modelHolder):
		self._modelHolder = modelHolder

	@property
	def modelFiles(self):
		return self._modelHolder

	def getPlaneWaveGeom(self,key):
		return getPlaneWaveGeom(key)
	

# Experimental Structure
def getExptStructAsUCell():
	''' From Walker 1959: DOI=10.1016/0001-6160(59)90090-2 Units of returned ucell are in bohr'''
	return getMgExptHcpAsUCell()

def _getMgExptHcpAsUCell():
	''' From Walker 1959: DOI=10.1016/0001-6160(59)90090-2 Units of returned ucell are in bohr'''
	
	cOverA = 1.624
	aVal = 3.20922 #In Angstrom
	lattParams = [aVal,aVal,aVal*cOverA]

	outCell = _getPerfectHcpMinimalUCell("Mg")
	outCell.setLattParams(lattParams)
	outCell.convAngToBohr()

	return outCell



def _getPerfectHcpMinimalUCell(element="Mg"):
	lattVects = [ [ 1.00000000, 0.00000000, 0.00000000],
	              [-0.50000000, 0.86602540, 0.00000000],
	              [ 0.00000000, 0.00000000, 1.00000000] ]
	
	fractCoords = [ [0.0,0.0,0.0,element], [ 1/3, 2/3, 0.5, element] ]
	idealCoverA = math.sqrt( 8/3 )

	perfectHcpCell = UCell.UnitCell.fromLattVects(lattVects, fractCoords=fractCoords)
	perfectHcpCell.setLattParams([1.0,1.0,idealCoverA])
	return perfectHcpCell


#INTERFACE FUNCTION
def getPlaneWaveGeom(structType:str):
	structTypeToFunct = {"hcp":_getMgPlaneWaveHcpGeomAsUCell,
	                     "bcc": _getMgPlaneWaveBCCGeomAsUCell,
	                     "fcc": _getMgPlaneWaveFCCGeomAsUCell}
	return structTypeToFunct[structType.lower()]()


# Optimised Plane-wave structures (PBE)
def _getMgPlaneWaveHcpGeomAsUCell():
	refFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/opt_geoms"
	refPath = os.path.join(refFolder,"hcp","Mg_hcp_SPE_otf_10el_usp_PP_6pt06.geom")
	initUCell = parseCastep.parseCastepGeomFile(refPath)[-1]["unitCell"]
	return initUCell

def _getMgPlaneWaveBCCGeomAsUCell():
	refFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/opt_geoms"
	refPath = os.path.join(refFolder, "bcc", "Mg_bcc_opt_otf_10el_usp_PP.geom")
	initUCell = parseCastep.parseCastepGeomFile(refPath)[-1]["unitCell"]
	return initUCell

def _getMgPlaneWaveFCCGeomAsUCell():
	refFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/opt_geoms"
	refPath = os.path.join(refFolder, "fcc", "Mg_fcc_opt_otf_10el_usp_PP.geom")
	initUCell = parseCastep.parseCastepGeomFile(refPath)[-1]["unitCell"]
	return initUCell


# Structures to use for bulk mod calculations (to be consistent)
def getUCellsForBulkModCalcs(structType:str):
	structTypeToFunct = {"hcp": _getHCPBulkModUCells,
	                     "bcc": _getBCCBulkModUCells,
	                     "fcc": _getFCCBulkModUCells}
	return structTypeToFunct[structType.lower()]()

def _getHCPBulkModUCells():
	baseRefFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/bulk_modulii/pure_mg"
	refFolder = os.path.join(baseRefFolder,"hcp")
	return _getUCellsFromBulkModFolder(refFolder)

def _getBCCBulkModUCells():
	baseRefFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/bulk_modulii/pure_mg"
	refFolder = os.path.join(baseRefFolder,"bcc")
	return _getUCellsFromBulkModFolder(refFolder)

def _getFCCBulkModUCells():
	baseRefFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/bulk_modulii/pure_mg"
	refFolder = os.path.join(baseRefFolder,"fcc")
	return _getUCellsFromBulkModFolder(refFolder)

def _getUCellsFromBulkModFolder(refFolder):
	casOutFiles = [os.path.join(refFolder,x) for x in os.listdir(refFolder) if x.endswith('.castep')]
	parsedUCells = [parseCastep.parseCastepOutfile(x)["unitCell"] for x in casOutFiles]
	[x.convAngToBohr() for x in parsedUCells]
	return parsedUCells

#Random structures:
def getHcpComprStructAsUCell():
	refFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/bulk_modulii/pure_mg"
	refPath = os.path.join(refFolder, "hcp", "Mg_hcp_SPE_otf_10el_usp_PP_5pt81.castep")
	initUCell = parseCastep.parseCastepOutfile(refPath)["unitCell"]
	initUCell.convAngToBohr()
	return initUCell


def getBccComprStructAsUCell():
	refFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/bulk_modulii/pure_mg"
	refPath = os.path.join(refFolder, "bcc", "Mg_bcc_opt_otf_10el_usp_PP_10_5pt592.castep")
	initUCell = parseCastep.parseCastepOutfile(refPath)["unitCell"]
	initUCell.convAngToBohr()
	return initUCell


#Plane wave bulk mod fits using ASE (Actually does the fit while function is called)
def getPlaneWaveEosFitDict(structType:str):
	structTypeToFunct = {"hcp": _getPlaneWaveEosDictHCP,
	                     "fcc": _getPlaneWaveEosDictFCC,
	                     "bcc": _getPlaneWaveEosDictBCC}
	return structTypeToFunct[structType]()

def _getPlaneWaveEosDictHCP():
	refFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/bulk_modulii/pure_mg"
	outFolder = os.path.join(refFolder,"hcp")
	fileList = [os.path.join(outFolder,x) for x in os.listdir(outFolder) if x.endswith('.castep')]
	return _getEosDictFromFilePaths(fileList)

def _getPlaneWaveEosDictFCC():
	refFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/bulk_modulii/pure_mg"
	outFolder = os.path.join(refFolder,"fcc")
	fileList = [os.path.join(outFolder,x) for x in os.listdir(outFolder) if x.endswith('.castep')]
	return _getEosDictFromFilePaths(fileList)

def _getPlaneWaveEosDictBCC():
	refFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/bulk_modulii/pure_mg"
	outFolder = os.path.join(refFolder,"bcc")
	fileList = [os.path.join(outFolder,x) for x in os.listdir(outFolder) if x.endswith('.castep')]
	return _getEosDictFromFilePaths(fileList)


def _getEosDictFromFilePaths(filePaths):
	outDict = fitBMod.getBulkModFromOutFilesAseWrapper(filePaths, eos="birchmurnaghan")
	return outDict

#Density of states for planeWave Calcs
def getDosPlaneWave(structType:str):
	structTypeToFunct = {"hcpExpt".lower(): _getHCPDosPlaneWaveExptGeom,
	                     "fcc": _getFCCDosPlaneWave,
	                     "bcc": _getBCCDosPlaneWave,
	                     "hcpCompr".lower(): _getHcpDosPlaneWaveComprGeomA,
	                     "bccCompr".lower():_getBccDosPlaneWaveComprGeomA}

	return structTypeToFunct[structType.lower()]()

def _getHCPDosPlaneWaveExptGeom():
	baseRefFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/some_dos/hcp_eqm"
	outFile = os.path.join(baseRefFolder, "Mg_hcp.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])

def _getFCCDosPlaneWave():
	baseRefFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/some_dos/fcc_eqm"
	outFile = os.path.join(baseRefFolder, "Mg_fcc.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])

def _getBCCDosPlaneWave():
	baseRefFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/some_dos/bcc_eqm"
	outFile = os.path.join(baseRefFolder, "Mg_bcc.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])


def _getHcpDosPlaneWaveComprGeomA():
	baseRefFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/some_dos/hcp_compr"
	outFile = os.path.join(baseRefFolder, "Mg_hcp.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])

def _getBccDosPlaneWaveComprGeomA():
	baseRefFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/some_dos/bcc_compr"
	outFile = os.path.join(baseRefFolder, "Mg_bcc.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])



#Structures to use for intersitials

def getInterstitialPlaneWaveStruct(structType:"str, e.g. hcp", interstitialType:"str, octahedral or tetrahedral",
                                   relaxType:"str, unrelaxed or relaxed", cellSize:"Str with dims, e.g 3_3_2"):

	paramsToStructDict = {("hcp","tetrahedral","unrelaxed","3_3_2"):_getHcpPlaneWaveStruct_interTetraUnrelaxed332(),
	                      ("hcp","octahedral" ,"unrelaxed","3_3_2"):_getHcpPlaneWaveStruct_interOctaUnrelaxed332(),
	                      ("hcp","octahedral", "relaxed_constant_pressure","3_3_2"):_getHcpPlaneWaveStruct_interOctaRelaxedConstPressure332()}

	return paramsToStructDict[(structType,interstitialType,relaxType,cellSize)]



def _getHcpPlaneWaveStruct_interTetraUnrelaxed332():
	refFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/interstitial/ref_files"
	refFile = os.path.join(refFolder,"mg_tetrahedral_interstitial_3_3_2.cell")
	return _getUCellFromCrystalMakerCastepOutFile(refFile)


def _getHcpPlaneWaveStruct_interOctaUnrelaxed332():
	refFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/interstitial/ref_files"
	refFile = os.path.join(refFolder,"mg_octahedral_interstitial_3_3_2.cell")
	return _getUCellFromCrystalMakerCastepOutFile(refFile)

def _getHcpPlaneWaveStruct_interOctaRelaxedConstPressure332():
	refFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/build_database/fermi_dirac_smearing/interstitial/relaxed/const_p"
	refFile = os.path.join(refFolder,"hcp_octa_inter.geom")
	outUCell = parseCastep.parseCastepGeomFile(refFile)[-1]["unitCell"]
	return outUCell

def _getUCellFromCrystalMakerCastepOutFile(refFile):
	tokenizedFile = parseCastep.tokenizeCastepCellFileAndRemoveBlockFromKeys(refFile)

	#Step 1 = get lattice vectors via the cell parameters. The fract co-ords should work fine for this (checked for octa supercell of hcp)
	cellParamStrList = tokenizedFile["lattice_abc"].strip().split()
	lattParams = [float(x) for x in cellParamStrList[1:4]]
	lattAngles = [float(x) for x in cellParamStrList[4:]]
	outUCell = UCell.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	outUCell.convAngToBohr()

	#Step 2 = get the fractional co-ords. Note I'm pretty sure this relies on me picking the "right" lattice vectors.
	outUCell.fractCoords = parseCastep._getFractCoordsFromTokenizedCellFile(tokenizedFile)

	return outUCell




