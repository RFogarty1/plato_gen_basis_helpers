
import os

import plato_pylib.plato.parse_inv_sk as parseInvSk

""" Purpose of this code is to make it easier to access inverse-SK data """



MG_PERFECT_CRYSTALS_DATA = dict()
ZR_PERFECT_CRYSTALS_DATA = dict()

#Despite the data being the same in both Mg folders, they seem to lead to SLIGHTLY different fits
#BASE_MG_DATA_FOLDER = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/testing_models/work_folder/inv_sk/all_terms/mg"
BASE_MG_DATA_FOLDER = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/opt_basis/10el_PP_dorbs/att6_longer_range/create_basis_sets/rc_7pt3/work_folder/inv_sk_calcs/screen1_no_damp"
BASE_ZR_DATA_FOLDER = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/testing_models/work_folder/inv_sk/all_terms/zr"

def registerFunctoToDict(key,regDict):
	def decorate(funct):
		regDict[key.lower()] = funct
		return funct
	return decorate



def getInverseSkPerfectCrystsSingleDataSet(eleKey:str, structType:str, removeXtal=True):
	""" See the relevant _PERFECT_CRYSTALS_DATA dict for available structType keys. """
	if eleKey.lower() == "mg":
		return getMgInverseSkPerfectCrystalsSingleDataset(structType, removeXtal)
	elif eleKey.lower() == "zr":
		return getZrInverseSkPerfectCrystalsSingleDataset(structType, removeXtal)
	raise ValueError("eleKey = {} is an invalid option".format(eleKey))



def getZrInverseSkPerfectCrystalsSingleDataset(structType:str, removeXtal=True):
	""" See Zr_PERFECT_CRYSTALS_DATA.keys() for available structType keys """
	outData = ZR_PERFECT_CRYSTALS_DATA[structType.lower()]()
	outData.removeXtalFieldTerms()
	return outData

def getMgInverseSkPerfectCrystalsSingleDataset(structType:str, removeXtal=True):
	""" See MG_PERFECT_CRYSTALS_DATA.keys() for available structType keys """
	outData = MG_PERFECT_CRYSTALS_DATA[structType.lower()]()
	outData.removeXtalFieldTerms()
	return outData

@registerFunctoToDict("hcp", MG_PERFECT_CRYSTALS_DATA)
def _getInvSkDataPerfectHcpCrystalsMg():
	outFolder = os.path.join(BASE_MG_DATA_FOLDER,"hcp")
	return parseAllInvSkOneFolder(outFolder)

@registerFunctoToDict("bcc", MG_PERFECT_CRYSTALS_DATA)
def _getInvSkDataPerfectHcpCrystalsMg():
	outFolder = os.path.join(BASE_MG_DATA_FOLDER,"bcc")
	return parseAllInvSkOneFolder(outFolder)

@registerFunctoToDict("fcc", MG_PERFECT_CRYSTALS_DATA)
def _getInvSkDataPerfectHcpCrystalsMg():
	outFolder = os.path.join(BASE_MG_DATA_FOLDER, "fcc")
	return parseAllInvSkOneFolder(outFolder)

@registerFunctoToDict("hcp", ZR_PERFECT_CRYSTALS_DATA)
def _getInvSkDataPefectHcpCrystalsZr():
	outFolder = os.path.join(BASE_ZR_DATA_FOLDER,"perfect_crystals", "hcp")
	return parseAllInvSkOneFolder(outFolder)

@registerFunctoToDict("fcc", ZR_PERFECT_CRYSTALS_DATA)
def _getInvSkDataPefectFccCrystalsZr():
	outFolder = os.path.join(BASE_ZR_DATA_FOLDER,"perfect_crystals", "fcc")
	return parseAllInvSkOneFolder(outFolder)

@registerFunctoToDict("bcc", ZR_PERFECT_CRYSTALS_DATA)
def _getInvSkDataPefectBccCrystalsZr():
	outFolder = os.path.join(BASE_ZR_DATA_FOLDER,"perfect_crystals", "bcc")
	return parseAllInvSkOneFolder(outFolder)


def parseAllInvSkOneFolder(inpFolder):
	parsedData = list()
	invSkFileNames = [x for x in os.listdir(inpFolder) if x.endswith("SK.csv")]

	for currFile in invSkFileNames:
		inpPath = os.path.join(inpFolder,currFile)
		parsedData.append( parseInvSk.parseInvSK(inpPath) )

	outData = parsedData[0]
	for x in parsedData:
		outData.addInvSKParsedFileData(x)
	return outData


