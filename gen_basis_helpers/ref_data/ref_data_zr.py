
""" Provides access to reference structures and data for pure Zr """

import os

import numpy as np 

from ..job_utils import dos_helpers as dosHelp
from . import helpers_ref_data as helpers
from . import ref_elemental_objs as refEleObjs
import plato_pylib.plato.mod_plato_inp_files as modInp
import plato_pylib.plato.parse_tbint_files as parseTbint
import plato_pylib.plato.plato_paths as platoPaths
import plato_pylib.shared.ucell_class as UCell

import plato_pylib.parseOther.parse_castep_files as parseCastep

tb1Model = os.path.join("Zr_models","two_body_2019") 
dft2Model = str(tb1Model)
dftModel = str(tb1Model) #Note havnet even got this stub version yet


BASE_FOLDER = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/testing_models/castep_database/zr"



def createZrReferenceDataObj():
	basePath = "/media/ssd1/rf614/Work/Plato"
	tb1ModAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(tb1Model)
	dft2ModelAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(dft2Model)
	dftModelAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(dftModel,dtype="dft")
	modelHolder = platoPaths.PlatoModelFolders(tb1Path=tb1ModAbs, dft2Path=dft2ModelAbs, dftPath=dftModelAbs)
	return ZrReferenceDataObj(modelHolder)

def createZrAngMomShellIndices():
	dft2ModelAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(dft2Model)
	dft2AdtPath = os.path.join(dft2ModelAbs, "Zr.adt")
	shellToAngMomDict = parseTbint.parseAdtFile(dft2AdtPath)["shellToAngMom"]
	
	angMomIndices = list()
	for key in range(len(shellToAngMomDict.keys())):
		angMomIndices.append( shellToAngMomDict[key] )

	return angMomIndices



#Overall object holding ref Data
class ZrReferenceDataObj(refEleObjs.RefElementalDataBase):

	def __init__(self, modelHolder):
		self._modelHolder = modelHolder

	@property
	def modelFiles(self):
		return self._modelHolder

	def getPlaneWaveGeom(self,key):
		return getPlaneWaveGeom(key)

	def getExptGeom(self,key):
		keyToFunct = {"hcp": getExptStructAsUCell}
		return keyToFunct[key]()


	def getStructsForEos(self,key):
		return getUCellsForBulkModCalcs(key)

	def getEosFitDict(self,key,eos="murnaghan"):
		return getPlaneWaveEosFitDict(key,eos=eos)


	def getSelfInterstitialPlaneWaveStruct(self, structType, interstitialType, relaxType, cellSize):
		return getInterstitialPlaneWaveStruct(structType, interstitialType, relaxType, cellSize)


	def getVacancyPlaneWaveStruct(self, structType, relaxType, cellSize):
		return getVacancyPlaneWaveStruct(structType, relaxType, cellSize)

	def getPlaneWaveDosData(self, structKey):
		return getDosPlaneWaveData(structKey)

	def getPlaneWaveDosGeom(self, structKey):
		return getDosPlaneWaveGeom(structKey)


def getExptStructAsUCell():
	''' a=3.233 A and c=5.182 A, making c/a = 1.063  - originally from W.B. Pearson, A Handbook of Lattice Spacings and Structures ofMetals (Pergamon Press, Oxford, 1967). '''

	aVal = 3.233
	cVal = 5.182
	lattParams = [aVal, aVal, cVal] 

	outCell = helpers.getPerfectHcpMinimalUCell("Zr")
	outCell.setLattParams(lattParams)
	outCell.convAngToBohr()

	return outCell


def getPlaneWaveGeom(structType:str):
	refFolder = os.path.join(BASE_FOLDER, "opt_geoms", structType)
	structTypeToFileName = {"hcp": "Zr_hcp_opt.castep",
	                        "bcc": "Zr_bcc_opt.castep",
	                        "fcc": "Zr_fcc_opt.castep"}
	refPath = os.path.join(refFolder, structTypeToFileName[structType])
	uCell = parseCastep.parseCastepOutfile(refPath)["unitCell"]
	uCell.convAngToBohr()
	return uCell


def getUCellsForBulkModCalcs(structType:str):
	refFolder = os.path.join(BASE_FOLDER, "eos", structType)
	return helpers.getUCellsFromCastepBulkModFolder(refFolder)


#Plane wave bulk mod fits using ASE (Actually does the fit while function is called)
def getPlaneWaveEosFitDict(structType:str, eos="murnaghan"):
	outFolder = os.path.join(BASE_FOLDER, "eos", structType)
	return helpers.getEosFitDictFromEosCastepFolder(outFolder)



def getInterstitialPlaneWaveStruct(structType:"str, e.g. hcp", interstitialType:"str, octahedral or tetrahedral",
                                   relaxType:"str, unrelaxed or relaxed", cellSize:"Str with dims, e.g 3_3_2"):
	paramsToStructDict = {("hcp","octahedral", "relaxed_constant_pressure","3_3_2"): _getHcpPlaneWaveStruct_interOctaRelaxedConstPressure332(),
	                      ("hcp","octahedral" ,"unrelaxed","3_3_2"): _getHcpPlaneWaveStruct_interOctaUnrelaxed(),
	                      ("hcp","tetrahedral","unrelaxed","3_3_2"):_getHcpPlaneWaveStruct_interTetraUnrelaxed332()}
	return paramsToStructDict[(structType,interstitialType,relaxType,cellSize)]


def _getHcpPlaneWaveStruct_interOctaRelaxedConstPressure332():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "relaxed", "constant_p", "Zr_hcp_strain6_0-Ointerstitial.castep")
	parsedUCell = parseCastep.parseCastepOutfile(refFile)["unitCell"]
	parsedUCell.convAngToBohr()
	return parsedUCell


def _getHcpPlaneWaveStruct_interOctaUnrelaxed():
	''' Data taken from first step of the relaxation optimisation, filename Zr_hcp_strain6_0-Ointerstitial.castep'''
	lattParams = [9.696, 9.696, 10.294] #Angstroms
	lattAngles = [90.0,90.0,120.0]
	fractCoords = [
	               [0.111110, 0.222223, 0.125000, "Zr"],
	               [0.111110, 0.222223, 0.625000, "Zr"],
	               [0.111110, 0.555557, 0.125000, "Zr"],
	               [0.111110, 0.555557, 0.625000, "Zr"],
	               [0.111110, 0.888890, 0.125000, "Zr"],
	               [0.111110, 0.888890, 0.625000, "Zr"],
	               [0.444443, 0.222223, 0.125000, "Zr"],
	               [0.444443, 0.222223, 0.625000, "Zr"],
	               [0.444443, 0.555557, 0.125000, "Zr"],
	               [0.444443, 0.555557, 0.625000, "Zr"],
	               [0.444443, 0.888890, 0.125000, "Zr"],
	               [0.444443, 0.888890, 0.625000, "Zr"],
	               [0.777777, 0.222223, 0.125000, "Zr"],
	               [0.777777, 0.222223, 0.625000, "Zr"],
	               [0.777777, 0.555557, 0.125000, "Zr"],
	               [0.777777, 0.555557, 0.625000, "Zr"],
	               [0.777777, 0.888890, 0.125000, "Zr"],
	               [0.777777, 0.888890, 0.625000, "Zr"],
	               [0.222223, 0.111110, 0.375000, "Zr"],
	               [0.222223, 0.111110, 0.875000, "Zr"],
	               [0.222223, 0.444443, 0.375000, "Zr"],
	               [0.222223, 0.444443, 0.875000, "Zr"],
	               [0.222223, 0.777777, 0.375000, "Zr"],
	               [0.222223, 0.777777, 0.875000, "Zr"],
	               [0.555557, 0.111110, 0.375000, "Zr"],
	               [0.555557, 0.111110, 0.875000, "Zr"],
	               [0.555557, 0.444443, 0.375000, "Zr"],
	               [0.555557, 0.444443, 0.875000, "Zr"],
	               [0.555557, 0.777777, 0.375000, "Zr"],
	               [0.555557, 0.777777, 0.875000, "Zr"],
	               [0.888890, 0.111110, 0.375000, "Zr"],
	               [0.888890, 0.111110, 0.875000, "Zr"],
	               [0.888890, 0.444443, 0.375000, "Zr"],
	               [0.888890, 0.444443, 0.875000, "Zr"],
	               [0.888890, 0.777777, 0.375000, "Zr"],
	               [0.888890, 0.777777, 0.875000, "Zr"],
	               [0.333333, 0.666667, 0.250000, "Zr"]
	              ]
	outUCell = UCell.UnitCell(lattParams=lattParams , lattAngles=lattAngles, fractCoords=fractCoords)
	outUCell.convAngToBohr()
	return outUCell


def _getHcpPlaneWaveStruct_interTetraUnrelaxed332():
	lattParams = [9.696, 9.696, 10.294] #Angstroms
	lattAngles = [90.0,90.0,120.0]
	fractCoords = [
	               [0.111110, 0.222223, 0.125000, "Zr"],
	               [0.111110, 0.222223, 0.625000, "Zr"],
	               [0.111110, 0.555557, 0.125000, "Zr"],
	               [0.111110, 0.555557, 0.625000, "Zr"],
	               [0.111110, 0.888890, 0.125000, "Zr"],
	               [0.111110, 0.888890, 0.625000, "Zr"],
	               [0.444443, 0.222223, 0.125000, "Zr"],
	               [0.444443, 0.222223, 0.625000, "Zr"],
	               [0.444443, 0.555557, 0.125000, "Zr"],
	               [0.444443, 0.555557, 0.625000, "Zr"],
	               [0.444443, 0.888890, 0.125000, "Zr"],
	               [0.444443, 0.888890, 0.625000, "Zr"],
	               [0.777777, 0.222223, 0.125000, "Zr"],
	               [0.777777, 0.222223, 0.625000, "Zr"],
	               [0.777777, 0.555557, 0.125000, "Zr"],
	               [0.777777, 0.555557, 0.625000, "Zr"],
	               [0.777777, 0.888890, 0.125000, "Zr"],
	               [0.777777, 0.888890, 0.625000, "Zr"],
	               [0.222223, 0.111110, 0.375000, "Zr"],
	               [0.222223, 0.111110, 0.875000, "Zr"],
	               [0.222223, 0.444443, 0.375000, "Zr"],
	               [0.222223, 0.444443, 0.875000, "Zr"],
	               [0.222223, 0.777777, 0.375000, "Zr"],
	               [0.222223, 0.777777, 0.875000, "Zr"],
	               [0.555557, 0.111110, 0.375000, "Zr"],
	               [0.555557, 0.111110, 0.875000, "Zr"],
	               [0.555557, 0.444443, 0.375000, "Zr"],
	               [0.555557, 0.444443, 0.875000, "Zr"],
	               [0.555557, 0.777777, 0.375000, "Zr"],
	               [0.555557, 0.777777, 0.875000, "Zr"],
	               [0.888890, 0.111110, 0.375000, "Zr"],
	               [0.888890, 0.111110, 0.875000, "Zr"],
	               [0.888890, 0.444443, 0.375000, "Zr"],
	               [0.888890, 0.444443, 0.875000, "Zr"],
	               [0.888890, 0.777777, 0.375000, "Zr"],
	               [0.888890, 0.777777, 0.875000, "Zr"],
	               [0.222222, 0.777778, 0.187500, "Zr"]
	              ]
	outUCell = UCell.UnitCell(lattParams=lattParams, lattAngles=lattAngles, fractCoords=fractCoords)
	outUCell.convAngToBohr()
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



def getDosPlaneWaveData(structType:str):

	structTypeToFunct = {"hcpExpt".lower(): _getHCPDosPlaneWaveExptGeom,
	                     "bcc": _getBCCDosPlaneWave,
	                     "hcpCompr".lower(): _getHcpDosPlaneWaveComprGeomA,
	                     "bccCompr".lower():_getBccDosPlaneWaveComprGeomA}

	return structTypeToFunct[structType.lower()]()



def _getHCPDosPlaneWaveExptGeom():
	outFile = os.path.join(BASE_FOLDER, "dos", "hcp_expt", "Zr_hcp.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])


def _getHcpDosPlaneWaveComprGeomA():
	outFile = os.path.join(BASE_FOLDER, "dos", "hcp_compr", "Zr_hcp.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])


def _getBCCDosPlaneWave():
	outFile = os.path.join(BASE_FOLDER, "dos", "bcc_eqm", "Zr_bcc.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])


def _getBccDosPlaneWaveComprGeomA():
	outFile = os.path.join(BASE_FOLDER, "dos", "bcc_compr", "Zr_bcc_compr.fixed.dat")
	return np.array(dosHelp.parseOptaDosDatFile(outFile)["dosdata"])


def getDosPlaneWaveGeom(structType:str):
	if structType == "bcc":
		return getPlaneWaveGeom("bcc")

	structTypeToFunct = {"hcpExpt".lower(): getExptStructAsUCell,
	                     "hcpCompr".lower(): _getHcpComprGeomForDos,
	                     "bccCompr".lower(): _getBccComprGeomForDos}

	return structTypeToFunct[structType.lower()]()

def _getHcpComprGeomForDos():
	outFile = os.path.join(BASE_FOLDER, "dos", "hcp_compr", "Zr_hcp_compr.castep")
	return helpers.getUCellInBohrFromCastepOutFile(outFile)


def _getBccComprGeomForDos():
	outFile = os.path.join(BASE_FOLDER, "dos", "bcc_compr", "Zr_bcc_compr_spe.castep")
	return helpers.getUCellInBohrFromCastepOutFile(outFile)

