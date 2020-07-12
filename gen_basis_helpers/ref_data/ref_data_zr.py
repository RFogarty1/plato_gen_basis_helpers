
""" Provides access to reference structures and data for pure Zr """

import itertools as it
import os
import types

import numpy as np 

from ..castep import castep_creator as castepCreator
from ..shared import config_vars as configVars
from ..job_utils import dos_helpers as dosHelp
from ..shared import unit_convs as unitConvs
from . import helpers_ref_data as helpers
from . import ref_elemental_objs as refEleObjs
import plato_pylib.plato.mod_plato_inp_files as modInp
import plato_pylib.plato.parse_tbint_files as parseTbint
import plato_pylib.plato.plato_paths as platoPaths
import plato_pylib.shared.ucell_class as UCell
import plato_pylib.utils.defects as defects
import plato_pylib.utils.supercell as supCell

import plato_pylib.parseOther.parse_castep_files as parseCastep

tb1Model = os.path.join("Zr_models","two_body_2019") 
dft2Model = str(tb1Model)
dftModel = str(tb1Model) #Note havnet even got this stub version yet


BASE_FOLDER =  os.path.join( configVars.CASTEP_DB_PATH, "zr" )



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

	def getStructsForEosWithPlaneWaveEnergies(self,key, perAtom=True, eUnits="eV"):
		return getUCellsAndEnergiesForBulkModCalcs(key,perAtom=perAtom, eUnits=eUnits)

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
		return getDosPlaneWaveData(structKey)

	def getPlaneWaveDosGeom(self, structKey):
		return getDosPlaneWaveGeom(structKey)

	def getPlaneWaveAtomTotalEnergy(self):
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



def getUCellsAndEnergiesForBulkModCalcs(structType,perAtom=True, eUnits="eV"):
	refFolder = os.path.join(BASE_FOLDER, "eos", structType)
	parsedFiles = [parseCastep.parseCastepOutfile(x) for x in helpers.getCastepOutPathsForFolder(refFolder)]
	allUCells = [x["unitCell"] for x in parsedFiles]
	[x.convAngToBohr() for x in allUCells]
	allEnergies = [x["energies"].electronicTotalE for x in parsedFiles]
	if perAtom:
		allEnergies = [ energy/parsed["numbAtoms"] for energy,parsed in it.zip_longest(allEnergies, parsedFiles) ]
	outObjs = list()
	for cell,energy in it.zip_longest(allUCells,allEnergies):
		currObj = types.SimpleNamespace(uCell=cell, energy=energy)
		outObjs.append(currObj)

	#Do any energy conversions needed
	if eUnits.lower()=="ev":
		pass
	elif eUnits.lower()=="ryd":
		print("Converting to Rydberg units")
		convUnits = unitConvs.EV_TO_RYD
		for x in outObjs:
			x.energy *= convUnits
	else:
		raise AttributeError("{} is an invalid value for eUnits".format(eUnits))

	return outObjs


#Plane wave bulk mod fits using ASE (Actually does the fit while function is called)
def getPlaneWaveEosFitDict(structType:str, eos="murnaghan"):
	outFolder = os.path.join(BASE_FOLDER, "eos", structType)
	return helpers.getEosFitDictFromEosCastepFolder(outFolder)



def getInterstitialPlaneWaveParsedFile(structType, interstitialType, relaxType, cellSize):
	paramsToOutFunct = {("hcp","no_inter", "plane_wave_geom", "3_3_2"): _getHcpPlaneWaveNoInter332,
	                    ("hcp","basal_tetra","relaxed_constant_pressure","3_3_2"): _getParsedFileHcpPlaneWaveBasalTetra_constantPressure,
	                    ("hcp","octahedral","relaxed_constant_pressure","3_3_2"): _getParsedFileHcpPlaneWaveOctahedral_constantPressure,
	                    ("hcp","split","relaxed_constant_pressure","3_3_2"): _getParsedFileHcpPlaneWaveSplitInter_constantPressure,
	                    ("hcp","basal_octa","relaxed_constant_pressure","3_3_2"):_getParsedFileHcpPlaneWaveBasalOctaInter_constantPressure,
	                    ("hcp","crowdion","relaxed_constant_pressure","3_3_2"):_getParsedFileHcpPlaneWAveCrowdionInter_constantPressure}
	return paramsToOutFunct[(structType,interstitialType,relaxType,cellSize)]()


def _getHcpPlaneWaveNoInter332():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "no_inter", "Zr_hcp_no_inter.castep")
	parsedFile = castepCreator.getParsedFileObjFromCastepOutputFile(refFile)
	return parsedFile

def _getParsedFileHcpPlaneWaveBasalTetra_constantPressure():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "relaxed", "constant_p", "basal_tetra", "geom_opt.castep")
	parsedFile = castepCreator.getParsedFileObjFromCastepOutputFile(refFile)
	return parsedFile

def _getParsedFileHcpPlaneWaveOctahedral_constantPressure():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "relaxed", "constant_p", "Zr_hcp_strain6_0-Ointerstitial.castep")
	parsedFile = castepCreator.getParsedFileObjFromCastepOutputFile(refFile)
	return parsedFile

def _getParsedFileHcpPlaneWaveSplitInter_constantPressure():
	#Originally was tetrahedral; but relaxed to split structure
	refFile = os.path.join(BASE_FOLDER,"interstitial", "relaxed", "constant_p", "Zr_hcp_strain6_0-Tinterstitial.castep")
	parsedFile = castepCreator.getParsedFileObjFromCastepOutputFile(refFile)
	return parsedFile

def _getParsedFileHcpPlaneWaveBasalOctaInter_constantPressure():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "relaxed", "constant_p", "basal_octa", "geom_opt.castep")
	parsedFile = castepCreator.getParsedFileObjFromCastepOutputFile(refFile)
	return parsedFile

def _getParsedFileHcpPlaneWAveCrowdionInter_constantPressure():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "relaxed", "constant_p", "crowdion", "geom_opt.castep")
	parsedFile = castepCreator.getParsedFileObjFromCastepOutputFile(refFile)
	return parsedFile

def getInterstitialPlaneWaveStruct(structType:"str, e.g. hcp", interstitialType:"str, octahedral or tetrahedral",
                                   relaxType:"str, unrelaxed or relaxed", cellSize:"Str with dims, e.g 3_3_2"):

	try:
		outVal = getInterstitialPlaneWaveParsedFile(structType, interstitialType, relaxType, cellSize).unitCell
		return outVal
	except KeyError:
		pass

	paramsToStructDict = {("hcp","octahedral" , "relaxed_constant_pressure","3_3_2"): _getHcpPlaneWaveStruct_interOctaRelaxedConstPressure332,
	                      ("hcp","tetrahedral", "relaxed_constant_pressure","3_3_2"): _getHcpPlaneWaveStruct_interTetraRelaxedConstantPressure332,
	                      ("hcp","octahedral" , "unrelaxed","3_3_2"): _getHcpPlaneWaveStruct_interOctaUnrelaxed,
	                      ("hcp","tetrahedral", "unrelaxed","3_3_2"):_getHcpPlaneWaveStruct_interTetraUnrelaxed332}
	return paramsToStructDict[(structType,interstitialType,relaxType,cellSize)]()


def _getHcpPlaneWaveStruct_interOctaRelaxedConstPressure332():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "relaxed", "constant_p", "Zr_hcp_strain6_0-Ointerstitial.castep")
	parsedUCell = parseCastep.parseCastepOutfile(refFile)["unitCell"]
	parsedUCell.convAngToBohr()
	return parsedUCell

def _getHcpPlaneWaveStruct_interTetraRelaxedConstantPressure332():
	refFile = os.path.join(BASE_FOLDER, "interstitial", "relaxed", "constant_p", "Zr_hcp_strain6_0-Tinterstitial.castep")
	parsedUCell = parseCastep.parseCastepOutfile(refFile)["unitCell"]
	parsedUCell.convAngToBohr()
	return parsedUCell


def _getHcpPlaneWaveStruct_interOctaUnrelaxed():
	''' Data taken from first step of the relaxation optimisation, filename Zr_hcp_strain6_0-Ointerstitial.castep'''
	lattParams = [9.696, 9.696, 10.294] #Angstroms
	lattAngles = [90.0,90.0,120.0]
	fractCoords = [
	               [0.111110, 0.222223, 0.125000],
	               [0.111110, 0.222223, 0.625000],
	               [0.111110, 0.555557, 0.125000],
	               [0.111110, 0.555557, 0.625000],
	               [0.111110, 0.888890, 0.125000],
	               [0.111110, 0.888890, 0.625000],
	               [0.444443, 0.222223, 0.125000],
	               [0.444443, 0.222223, 0.625000],
	               [0.444443, 0.555557, 0.125000],
	               [0.444443, 0.555557, 0.625000],
	               [0.444443, 0.888890, 0.125000],
	               [0.444443, 0.888890, 0.625000],
	               [0.777777, 0.222223, 0.125000],
	               [0.777777, 0.222223, 0.625000],
	               [0.777777, 0.555557, 0.125000],
	               [0.777777, 0.555557, 0.625000],
	               [0.777777, 0.888890, 0.125000],
	               [0.777777, 0.888890, 0.625000],
	               [0.222223, 0.111110, 0.375000],
	               [0.222223, 0.111110, 0.875000],
	               [0.222223, 0.444443, 0.375000],
	               [0.222223, 0.444443, 0.875000],
	               [0.222223, 0.777777, 0.375000],
	               [0.222223, 0.777777, 0.875000],
	               [0.555557, 0.111110, 0.375000],
	               [0.555557, 0.111110, 0.875000],
	               [0.555557, 0.444443, 0.375000],
	               [0.555557, 0.444443, 0.875000],
	               [0.555557, 0.777777, 0.375000],
	               [0.555557, 0.777777, 0.875000],
	               [0.888890, 0.111110, 0.375000],
	               [0.888890, 0.111110, 0.875000],
	               [0.888890, 0.444443, 0.375000],
	               [0.888890, 0.444443, 0.875000],
	               [0.888890, 0.777777, 0.375000],
	               [0.888890, 0.777777, 0.875000],
	               [0.333333, 0.666667, 0.250000]
	              ]
	elementList = ["Zr" for x in fractCoords]
	outUCell = UCell.UnitCell(lattParams=lattParams, lattAngles=lattAngles, fractCoords=fractCoords, elementList=elementList)
	outUCell.convAngToBohr()
	return outUCell


def _getHcpPlaneWaveStruct_interTetraUnrelaxed332():
	lattParams = [9.696, 9.696, 10.294] #Angstroms
	lattAngles = [90.0,90.0,120.0]
	fractCoords = [
	               [0.111110, 0.222223, 0.125000],
	               [0.111110, 0.222223, 0.625000],
	               [0.111110, 0.555557, 0.125000],
	               [0.111110, 0.555557, 0.625000],
	               [0.111110, 0.888890, 0.125000],
	               [0.111110, 0.888890, 0.625000],
	               [0.444443, 0.222223, 0.125000],
	               [0.444443, 0.222223, 0.625000],
	               [0.444443, 0.555557, 0.125000],
	               [0.444443, 0.555557, 0.625000],
	               [0.444443, 0.888890, 0.125000],
	               [0.444443, 0.888890, 0.625000],
	               [0.777777, 0.222223, 0.125000],
	               [0.777777, 0.222223, 0.625000],
	               [0.777777, 0.555557, 0.125000],
	               [0.777777, 0.555557, 0.625000],
	               [0.777777, 0.888890, 0.125000],
	               [0.777777, 0.888890, 0.625000],
	               [0.222223, 0.111110, 0.375000],
	               [0.222223, 0.111110, 0.875000],
	               [0.222223, 0.444443, 0.375000],
	               [0.222223, 0.444443, 0.875000],
	               [0.222223, 0.777777, 0.375000],
	               [0.222223, 0.777777, 0.875000],
	               [0.555557, 0.111110, 0.375000],
	               [0.555557, 0.111110, 0.875000],
	               [0.555557, 0.444443, 0.375000],
	               [0.555557, 0.444443, 0.875000],
	               [0.555557, 0.777777, 0.375000],
	               [0.555557, 0.777777, 0.875000],
	               [0.888890, 0.111110, 0.375000],
	               [0.888890, 0.111110, 0.875000],
	               [0.888890, 0.444443, 0.375000],
	               [0.888890, 0.444443, 0.875000],
	               [0.888890, 0.777777, 0.375000],
	               [0.888890, 0.777777, 0.875000],
	               [0.222222, 0.777778, 0.187500]
	              ]
	elementList = ["Zr" for x in fractCoords]
	outUCell = UCell.UnitCell(lattParams=lattParams, lattAngles=lattAngles, fractCoords=fractCoords, elementList=elementList)
	outUCell.convAngToBohr()
	return outUCell

def getInterstitialPlaneWaveFormationEnergy(structType, interstitialType, relaxType, cellSize):

	paramsToEnergyDict = {("hcp","tetrahedral","unrelaxed","3_3_2"): _getHcpPlaneWaveFormationEnergy_interTetraUnrelaxed332,
	                      ("hcp","octahedral" ,"unrelaxed","3_3_2"): _getHcpPlaneWaveFormationEnergy_interOctaUnrelaxed332,
	                      ("hcp","octahedral", "relaxed_constant_pressure","3_3_2"): _getHcpPlaneWaveFormationEnergy_interOctaRelaxedConstPressure332,
	                      ("hcp","tetrahedral", "relaxed_constant_pressure","3_3_2"): _getHcpPlaneWaveFormationEnergy_interTetraRelaxedConstantPressure332}

	return paramsToEnergyDict[(structType,interstitialType,relaxType,cellSize)]() 


def _getHcpPlaneWaveFormationEnergy_interTetraUnrelaxed332():
	interstitFile = os.path.join(BASE_FOLDER, "interstitial", "unrelaxed", "Zr_hcp_unrelaxed_tetra_inter.castep")
	noInterFile = os.path.join(BASE_FOLDER, "interstitial", "no_inter", "Zr_hcp_no_inter.castep")
	defectEnergy = helpers.getDefectEnergyFromCastepNoDefectAndDefectFiles(noInterFile, interstitFile)
	return defectEnergy	


def _getHcpPlaneWaveFormationEnergy_interOctaUnrelaxed332():
	interstitFile = os.path.join(BASE_FOLDER, "interstitial", "unrelaxed", "Zr_hcp_unrelaxed_octa_inter.castep")
	noInterFile = os.path.join(BASE_FOLDER, "interstitial", "no_inter", "Zr_hcp_no_inter.castep")
	defectEnergy = helpers.getDefectEnergyFromCastepNoDefectAndDefectFiles(noInterFile, interstitFile)
	return defectEnergy


def _getHcpPlaneWaveFormationEnergy_interOctaRelaxedConstPressure332():
	interstitFile = os.path.join(BASE_FOLDER, "interstitial", "relaxed", "constant_p", "Zr_hcp_strain6_0-Ointerstitial.castep")
	noInterFile = os.path.join(BASE_FOLDER, "interstitial", "no_inter", "Zr_hcp_no_inter.castep")
	defectEnergy = helpers.getDefectEnergyFromCastepNoDefectAndDefectFiles(noInterFile, interstitFile)
	return defectEnergy

def _getHcpPlaneWaveFormationEnergy_interTetraRelaxedConstantPressure332():
	interstitFile = os.path.join(BASE_FOLDER, "interstitial", "relaxed", "constant_p", "Zr_hcp_strain6_0-Tinterstitial.castep")
	noInterFile = os.path.join(BASE_FOLDER, "interstitial", "no_inter", "Zr_hcp_no_inter.castep")
	defectEnergy = helpers.getDefectEnergyFromCastepNoDefectAndDefectFiles(noInterFile, interstitFile)
	return defectEnergy


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
		raise NotImplementedError("relaxType = {} not currently implemented".format(relaxType))


def _getVacancyHcpPlaneWaveStruct_relaxedConstantPressure():
	vacFilePath = os.path.join(BASE_FOLDER, "vacancy", "relaxed", "constant_p", "Zr_hcp_strain6_0-vacancy.castep")
	return helpers.getUCellInBohrFromCastepOutFile(vacFilePath)


def getVacancyPlaneWaveEnergy(structType, relaxType, cellSize):
	paramsToEnergyDict = {("hcp", "unrelaxed", "3_3_2"): _getVacancyEnergyHcpUnrelaxed332,
	                      ("hcp", "relaxed_constant_pressure","3_3_2"): _getVacancyEnergyHcpRelaxedConstantPressure332}
	return paramsToEnergyDict[(structType, relaxType, cellSize)]()


def _getVacancyEnergyHcpUnrelaxed332():
	vacFilePath = os.path.join(BASE_FOLDER, "vacancy", "unrelaxed", "3_3_2", "Zr_hcp_unrelaxed_vacancy.castep" )
	bulkFilePath = os.path.join(BASE_FOLDER, "interstitial", "no_inter", "Zr_hcp_no_inter.castep")
	outEnergy = helpers.getDefectEnergyFromCastepNoDefectAndDefectFiles(bulkFilePath, vacFilePath)
	return outEnergy

def _getVacancyEnergyHcpRelaxedConstantPressure332():
	vacFilePath = os.path.join(BASE_FOLDER, "vacancy", "relaxed", "constant_p", "Zr_hcp_strain6_0-vacancy.castep")
	bulkFilePath = os.path.join(BASE_FOLDER, "interstitial", "no_inter", "Zr_hcp_no_inter.castep")
	return helpers.getDefectEnergyFromCastepNoDefectAndDefectFiles(bulkFilePath, vacFilePath)


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

	structTypeToFunct = {"hcpExpt".lower(): getExptStructAsUCell,
	                     "hcpCompr".lower(): _getHcpComprGeomForDos,
	                     "bccCompr".lower(): _getBccComprGeomForDos,
	                     "bcc": _getBccEqmGeomForDos}

	return structTypeToFunct[structType.lower()]()

def _getHcpComprGeomForDos():
	outFile = os.path.join(BASE_FOLDER, "dos", "hcp_compr", "Zr_hcp_compr.castep")
	return helpers.getUCellInBohrFromCastepOutFile(outFile)


def _getBccComprGeomForDos():
	outFile = os.path.join(BASE_FOLDER, "dos", "bcc_compr", "Zr_bcc_compr_spe.castep")
	return helpers.getUCellInBohrFromCastepOutFile(outFile)

def _getBccEqmGeomForDos():
	outFile = os.path.join(BASE_FOLDER, "dos", "bcc_eqm", "Zr_bcc_eqm_spe.castep")
	return helpers.getUCellInBohrFromCastepOutFile(outFile)


def getPlaneWaveDissocSepVsTotalE(inBohr=True):
	outFolder = os.path.join(BASE_FOLDER,"dissoc_curve")
	outFiles = helpers.getCastepOutPathsForFolder(outFolder)
	outEnergies = [parseCastep.parseCastepOutfile(x)["energies"].electronicTotalE for x in outFiles]
	outSeps = [helpers.getDimerSepFromCastepOutFile(x,inBohr=True) for x in outFiles]
	outList = list()
	for sep,e in it.zip_longest(outSeps,outEnergies):
		outList.append( [sep,e] ) 

	return outList


def getPlaneWaveSurfaceParsedFileObject(surfType, nLayers=None, relaxType="unrelaxed"):
	""" Returns a ParsedFile object for castep calculation on a surface. Note the unrelaxed geometries are built from from the castep optimised bulk cell
	
	Args:
		surfType (str): The type of surface; examples are hcp0001 and hcp10m10 
		nLayers (int, optional): The number of surface layers used in the calculation. Default will vary for different surfTypes; but will represent converged structures
		relaxType (str, optional): String denoting the type of relaxation applied. Default="unrelaxed"; other options are "constant_volume" for now

	Returns
		parsedFile (ParsedFile object): This contains the geometry and total energy of the requested structure
	
	"""
	structTypeDefaultNLayers = {"hcp0001":10}
	if nLayers is None:
		nLayers = structTypeDefaultNLayers[surfType]

	structTypeToFunct = {
	                      ("hcp0001" , 10,"constant_volume"  ): _getZrHcp0001PlaneWaveRelaxedParsedFile_10layers
	                     }

	return structTypeToFunct[(surfType,nLayers,relaxType)]()


def _getZrHcp0001PlaneWaveRelaxedParsedFile_10layers():
    refPath = os.path.join(BASE_FOLDER,"surface_energies", "hcp0001", "castep_geom", "relaxed", "nlayers_10_absvac_10_ang", "geom_opt.castep")
    return castepCreator.getParsedFileObjFromCastepOutputFile(refPath)


def getPlaneWaveSurfEnergy(structKey):
	outDict = {"hcp0001": _getSurfaceEnergyHcp0001}
	return outDict[structKey.lower()]()


def _getSurfaceEnergyHcp0001():
	refBaseFolder = os.path.join(BASE_FOLDER,"surface_energies", "hcp0001" )
	bulkFilePath = os.path.join( BASE_FOLDER,"dos", "hcp_expt", "Zr_hcp_expt.castep" )
	surfFile = os.path.join(refBaseFolder, "Zr24.castep")
	return helpers.getCastepRefHcp0001SurfaceEnergyFromSurfAndBulkFilePaths(surfFile,bulkFilePath)
