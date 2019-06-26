
''' My attempt to factor out calculation methods to a single methodStr '''

import os
import pathlib

import plato_pylib.plato.parse_plato_out_files as parsePlato 
import plato_pylib.plato.mod_plato_inp_files as platoInp
import plato_pylib.utils.job_running_functs as jobRun

import evol_test_vs_ref_helpers as evolHelp
import dft_file_helpers as dftHelpers


METHOD_STR_TO_OBJ = dict()
ALL_METHOD_STRS = set()


#Basic decorator factories
def registerMethodStrToObj(key):
	def decorate(funct):
		METHOD_STR_TO_OBJ[key] = funct
		ALL_METHOD_STRS.add(key)
		return funct
	return decorate

#Method class should encapsulate all the things that vary between different methods in plato
class PlatoMethod():

	def __init__(self, optDict, runCommFunction, strDictFromOptDictFunction, gridKwarg=None):
		self.optDict = optDict
		self.runCommFunction = runCommFunction #interface is getRunComm(fileName)
		self.getStrDictFromOptDict = strDictFromOptDictFunction #interface is getStrDictFromOptDict
		self.gridKwarg = "IntegralMeshSpacing".lower() if gridKwarg is None else gridKwarg #Different for dft case annoyingly
#		self.gridKwarg = "IntegralMeshSpacing".lower()

	@property
	def kpts(self):
		return self.optDict["BlochStates".lower()]

	@kpts.setter
	def kpts(self,value):
		self.optDict["BlochStates".lower()] = value

	@property
	def dataSet(self):
		return self.optDict["dataset"]

	@dataSet.setter
	def dataSet(self,value:"str, plato relative path"):
		self.optDict["dataset"] = value

	@property
	def integGrid(self):
		return self.optDict[self.gridKwarg]

	@integGrid.setter
	def integGrid(self,value):
		self.optDict[self.gridKwarg] = value

	def getStrDictWithStruct(self, struct:"ucell obj"):
		strDict = self.getStrDictFromOptDict(self.optDict)
		geomDict = platoInp.getPlatoGeomDictFromUnitCell(struct)
		strDict.update(geomDict)
		return strDict


#Different plato Run Comms which get attached to method objects
def _getRunCommPlatoDft2(inpFilePath):
	runComm = jobRun.pathListToPlatoRunComms([inpFilePath],"dft2")[0]
	return runComm

def _getRunCommPlatoTb1(inpFilePath):
	runComm = jobRun.pathListToPlatoRunComms([inpFilePath],"tb1")[0]
	return runComm

def _getRunCommPlatoDft(inpFilePath):
	runComm = jobRun.pathListToPlatoRunComms([inpFilePath],"dft")[0]
	return runComm



def getPlatoCalcObjFromInpPathAndStrDictAndRunCommFunction(inpFilePath, strDict, runCommFunction):
	strDictWriteFunction = platoInp.writePlatoOutFileFromDict
	fileParser = parsePlato.parsePlatoOutFile_energiesInEv
	return CalcObj.fromEnforcedKwargs(filePath=inpFilePath, strDict=strDict, strDictWriteFunction=strDictWriteFunction,
	                                  fileParser=fileParser, runCommFunction=runCommFunction, inpFileExt=".in",
	                                  outFileExt = ".out")

class CalcObj():
	def __init__(self, filePath, strDict, strDictWriteFunction, fileParser, runCommFunction, inpFileExt, outFileExt):
		self.filePath = os.path.splitext(filePath)[0]
		self.strDict = strDict
		self.inpFileExt = inpFileExt
		self.outFileExt = outFileExt
		self._strDictWriteFunction = strDictWriteFunction
		self._fileParser = fileParser
		self._runCommFunction = runCommFunction

	@classmethod
	def fromEnforcedKwargs(cls, **kwargs):
		kwargs = {k.lower():v for k,v in kwargs.items()}
		reqArgsOrder = ["filePath", "strDict", "strDictWriteFunction", "fileParser",
		                "runCommFunction", "inpFileExt", "outFileExt"]
		argsList = list()
		for x in reqArgsOrder:
			argsList.append( kwargs[x.lower()] )
	
		return cls(*argsList)

	def getRunComm(self):
		return self._runCommFunction(self.filePath + self.inpFileExt)

	def writeFile(self):
		writeFolder = os.path.split(self.filePath)[0]
		pathlib.Path(writeFolder).mkdir(exist_ok=True, parents=True)
		self._strDictWriteFunction(self.filePath + self.inpFileExt, self.strDict)

	def parseOutFile(self):
		return self._fileParser(self.filePath + self.outFileExt)

#Creating PlatoMethod objects
def createPlatoMethodObj(methodStr):
	return METHOD_STR_TO_OBJ[methodStr]()

@registerMethodStrToObj("tb1_2c_only")
def createPlatoMethod_tb1_2c_only():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getStrDict_dft2_tb1_2c_only()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getStrDict_dft2_tb1_2c_only():
	return evolHelp.loadDefaultTb1OptDict()


@registerMethodStrToObj("tb1_2c_noxtalxc")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_2c_noxtalxc()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_2c_noxtalxc():
	optDict = evolHelp.loadDefaultTb1OptDict()
	optDict["CrystalFieldXCWeight".lower()] = 0.0
	return optDict



@registerMethodStrToObj("tb1_vxc_sncorr3_exc_corr1")
def createPlatoMethod_tb1_vxc_sncorr4_exc_corr1():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_vxc_sncorr3_exc_corr1()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_vxc_sncorr3_exc_corr1():
	outOptDict = evolHelp.loadDefaultTb1OptDict()
	modOptsDict = {"excmbcorr": [1, 1.0],
				  "VxcMBCorrXtal".lower(): [3,1.0]}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_2centxc_noxtal")
def createPlatoMethod_dft2_2centxc_noxtal():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2centxc_noxtal()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)


def _getOptDict_dft2_2centxc_noxtal():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():1,
				  "excMbCorr".lower():0,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_2centxc_withxtal")
def createPlatoMethod_dft2_2centxc_noxtal():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2centxc_withxtal()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2centxc_withxtal():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():3,
				  "excMbCorr".lower():0,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_mcweda4_pp")
def createPlatoMethod_dft2_mcweda4_pp():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_mcweda4_pp()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)


def _getOptDict_dft2_mcweda4_pp():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():4,
				  "excMbCorr".lower():0,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict

@registerMethodStrToObj("dft2_mcweda4_pp_uniform")
def createPlatoMethod_dft2_mcweda4_pp_uniform():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_mcweda4_pp_uniform()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)


def _getOptDict_dft2_mcweda4_pp_uniform():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():4,
				  "excMbCorr".lower():1,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_mcweda4_pp_uniform_second_order")
def createPlatoMethod_dft2_mcweda4_pp_uniform_second_order():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_dft2_mcweda4_pp_uniform_second_order()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)


def _getOptDict_dft2_dft2_mcweda4_pp_uniform_second_order():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():4,
				  "excMbCorr".lower():3,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_mcweda4_exact_e0")
def createPlatoMethod_dft2_mcweda4_exact_e0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_mcweda4_exact_e0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_mcweda4_exact_e0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():4,
				  "excMbCorr".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict

@registerMethodStrToObj("dft2_mcweda5_exact_e0")
def createPlatoMethod_dft2_mcweda5_exact_e0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_mcweda5_exact_e0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_mcweda5_exact_e0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():5,
				  "excMbCorr".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_hop_mcweda_xtal_halfmcweda_halfsn2b_exact_e0")
def createPlatoMethod_dft2_hop_mcweda_xtal_halfmcweda_halfsn2b_exact_e0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_hop_mcweda_xtal_halfmcweda_halfsn2b_exact_e0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_hop_mcweda_xtal_halfmcweda_halfsn2b_exact_e0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): [5,2,*[3,4],*[0.5,0.5]], #[flag, numbMethodsToMix, [flags for mixnig methods], [scaling factors] ]
	               "hopxcmethod".lower():2,
				  "excMbCorr".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict

@registerMethodStrToObj("dft2_mcwedahop_exactxtal_exact_e0")
def createPlatoMethod_dft2_mcwedahop_exactxtal_exact_e0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_mcwedahop_exactxtal_exact_e0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_mcwedahop_exactxtal_exact_e0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 0,
	               "hopxcmethod".lower():2,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_2bodyhop_exactxtal_exacte0")
def createPlatoMethod_dft2_2bodyhop_exactxtal_exacte0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bodyhop_exactxtal_exacte0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bodyhop_exactxtal_exacte0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 0,
	               "hopxcmethod".lower():1,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_exacthop_noxtal_exacte0")
def createPlatoMethod_dft2_exacthop_noxtal_exacte0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exacthop_noxtal_exacte0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exacthop_noxtal_exacte0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 1,
	               "hopxcmethod".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_exacthop_2bxtal_exacte0")
def createPlatoMethod_dft2_exacthop_2bxtal_exacte0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exacthop_2bxtal_exacte0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exacthop_2bxtal_exacte0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 2,
	               "hopxcmethod".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_exacthop_2bxtal_sncorr_exacte0")
def createPlatoMethod_dft2_exacthop_2bxtal_sncorr_exacte0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exacthop_2bxtal_sncorr_exacte0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exacthop_2bxtal_sncorr_exacte0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 3,
	               "hopxcmethod".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_exacthop_mcweda_xtal_exacte0")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exacthop_mcweda_xtal_exacte0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exacthop_mcweda_xtal_exacte0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 4,
	               "hopxcmethod".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_2bodyhop_2bxtal_exacte0")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bodyhop_2bxtal_exacte0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bodyhop_2bxtal_exacte0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 2,
	               "hopxcmethod".lower():1,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_2bodyhop_mcwedaxtal_exacte0")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bodyhop_mcwedaxtal_exacte0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bodyhop_mcwedaxtal_exacte0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 4,
	               "hopxcmethod".lower():1,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_exact_e1_pp")
def createPlatoMethod_dft2_exact():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exact_e1_pp()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exact_e1_pp():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():0,
				  "excMbCorr".lower():0,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_2body_e1_exact_e0")
def createPlatoMethod_dft2_exact():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2body_e1_exact_e0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2body_e1_exact_e0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():3,
				  "excMbCorr".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_exact_e1_pp_uniform_e0")
def createPlatoMethod_dft2_exact_e1_pp_uniform_e0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exact_e1_pp_uniform_e0()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exact_e1_pp_uniform_e0():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():0,
				  "excMbCorr".lower():1,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict




@registerMethodStrToObj("dft2_exact")
def createPlatoMethod_dft2_exact():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exact()
	getStrDictFromOptDict = evolHelp.getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exact():
	outOptDict = evolHelp.loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():0,
				  "excMbCorr".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft_plato")
def createPlatoMethod_dft_plato():
	runCommFunct = _getRunCommPlatoDft
	optDict = _getOptDict_dft_plato()
	getStrDictFromOptDict = dftHelpers.getPlatoStrDictFromOptDict_dft
	gridKwarg = "fftGridSpacing".lower()
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict, gridKwarg=gridKwarg)


def _getOptDict_dft_plato():
	outOptDict = dftHelpers.loadDefaultDftOptDict()
	modOptDict = {"xcfunctional": "pbe",
	              "scfflag": 1,
	              "optimizemesh": 0} #Some of these MIGHT be defaults
	outOptDict.update(modOptDict)
	return outOptDict






