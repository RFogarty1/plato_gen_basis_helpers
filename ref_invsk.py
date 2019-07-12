
import os

import plato_pylib.plato.parse_inv_sk as parseInvSk

""" Purpose of this code is to make it easier to access inverse-SK data """



MG_PERFECT_CRYSTALS_DATA = dict()

def registerFunctoToDict(key,regDict):
	def decorate(funct):
		regDict[key.lower()] = funct
		return funct
	return decorate



def getMgInverseSkPerfectCrystalsSingleDataset(structType:str, removeXtal=True):
	""" See MG_PERFECT_CRYSTALS_DATA.keys() for available structType keys """
	outData = MG_PERFECT_CRYSTALS_DATA[structType.lower()]()
	outData.removeXtalFieldTerms()
	return outData

@registerFunctoToDict("hcp", MG_PERFECT_CRYSTALS_DATA)
def _getInvSkDataPerfectHcpCrystals():
	baseFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/opt_basis/10el_PP_dorbs/att6_longer_range/create_basis_sets/rc_7pt3/work_folder/inv_sk_calcs/screen1_no_damp"
	actFolder = os.path.join(baseFolder,"hcp")
	return parseAllInvSkOneFolder(actFolder)

@registerFunctoToDict("bcc", MG_PERFECT_CRYSTALS_DATA)
def _getInvSkDataPerfectHcpCrystals():
	baseFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/opt_basis/10el_PP_dorbs/att6_longer_range/create_basis_sets/rc_7pt3/work_folder/inv_sk_calcs/screen1_no_damp"
	actFolder = os.path.join(baseFolder,"bcc")
	return parseAllInvSkOneFolder(actFolder)


@registerFunctoToDict("fcc", MG_PERFECT_CRYSTALS_DATA)
def _getInvSkDataPerfectHcpCrystals():
	baseFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/Building_Mg_Model/opt_basis/10el_PP_dorbs/att6_longer_range/create_basis_sets/rc_7pt3/work_folder/inv_sk_calcs/screen1_no_damp"
	actFolder = os.path.join(baseFolder,"fcc")
	return parseAllInvSkOneFolder(actFolder)


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


