
''' My attempt to factor out calculation methods to a single methodStr '''

import os

import plato_pylib.plato.parse_plato_out_files as parsePlato 
import plato_pylib.plato.mod_plato_inp_files as platoInp
import plato_pylib.utils.job_running_functs as jobRun


from .plato_calc_objs import PlatoMethod, CalcObj


METHOD_STR_TO_OBJ = dict()
ALL_METHOD_STRS = set()


#TODO: Would be good to remove these functions eventually(replace with funct partial calls probably)
def getPlatoStrDictFromOptDict_tb1OrTb2(optDict):
	return platoInp.getStrDictFromOptDict(optDict,"dft2")

def getPlatoStrDictFromOptDict_dft(optDict):
	return platoInp.getStrDictFromOptDict(optDict,"dft")

def loadDefaultTb1OptDict():
	return platoInp.getDefOptDict("tb1")

def loadDefaultTb2OptDict():
	return platoInp.getDefOptDict("tb2")

def loadDefaultDftOptDict():
	return platoInp.getDefOptDict("dft")

#Basic decorator factories
def registerMethodStrToObj(key):
	def decorate(funct):
		METHOD_STR_TO_OBJ[key] = funct
		ALL_METHOD_STRS.add(key)
		return funct
	return decorate


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


#Creating PlatoMethod objects
def createPlatoMethodObj(methodStr):
	return METHOD_STR_TO_OBJ[methodStr]()

@registerMethodStrToObj("tb1_2c_only")
def createPlatoMethod_tb1_2c_only():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getStrDict_dft2_tb1_2c_only()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getStrDict_dft2_tb1_2c_only():
	return loadDefaultTb1OptDict()


@registerMethodStrToObj("tb1_2c_noxtalxc")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_2c_noxtalxc()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_2c_noxtalxc():
	optDict = loadDefaultTb1OptDict()
	optDict["CrystalFieldXCWeight".lower()] = 0.0
	return optDict


@registerMethodStrToObj("tb1_2bxtal_plus_sncorr")
def createPlatoMethod_tb1_2bxtal_plus_sncorr():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_2bxtal_plus_sncorr()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_2bxtal_plus_sncorr():
	outOptDict = loadDefaultTb1OptDict()
	modOptsDict = {"VxcMBCorrXtal".lower(): [3,1.0]}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("tb1_snxtal_2bmb")
def createPlatoMethod_tb1_snxtal_2bmb():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_snxtal_2bmb()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_snxtal_2bmb():
	outOptDict = loadDefaultTb1OptDict()
	modOptsDict = {"VxcMBCorrXtal".lower(): [4,1.0],
	               "CrystalFieldXCWeight".lower(): 0}
	outOptDict.update(modOptsDict)
	return outOptDict

@registerMethodStrToObj("tb1_2c_fittedhop")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_2c_fittedHop()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_2c_fittedHop():
	optDict = loadDefaultTb1OptDict()
	optDict["addCorrectingHopFromBdt".lower()] = 1
	return optDict


@registerMethodStrToObj("tb1_2bxtal_plus_sncorr_fittedhop")
def createPlatoMethod_tb1_2bxtal_plus_sncorr_fittedhop():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_2bxtal_plus_sncorr_fittedhop()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_2bxtal_plus_sncorr_fittedhop():
	outOptDict = loadDefaultTb1OptDict()
	modOptsDict = {"VxcMBCorrXtal".lower(): [3,1.0]}
	modOptsDict["addCorrectingHopFromBdt".lower()] = 1
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("tb1_snxtal_2bmb_fittedHop")
def createPlatoMethod_tb1_snxtal_2bmb_fittedHop():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_snxtal_2bmb_fittedHop()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_snxtal_2bmb_fittedHop():
	outOptDict = loadDefaultTb1OptDict()
	modOptsDict = {"VxcMBCorrXtal".lower(): [4,1.0],
	               "CrystalFieldXCWeight".lower(): 0}
	modOptsDict["addCorrectingHopFromBdt".lower()] = 1
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("tb1_snxtal_2bmb_fittedHop_fittedPP")
def createPlatoMethod_tb1_snxtal_2bmb_fittedHop_fittedPP():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_tb1_snxtal_2bmb_fittedHop_fittedPP()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_tb1_snxtal_2bmb_fittedHop_fittedPP():
	outOptDict = loadDefaultTb1OptDict()
	modOptsDict = {"VxcMBCorrXtal".lower(): [4,1.0],
	               "CrystalFieldXCWeight".lower(): 0}
	modOptsDict["addCorrectingHopFromBdt".lower()] = 1
	modOptsDict["addCorrectingPPFromBdt".lower()] = 1
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("tb1_snxtal_2bmb_fittedHop_uniformdens_pp")
def createPlatoMethod_tb1_snxtal_2bmb_fittedHop_uniformdens():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_snxtal_2bmb_fittedHop_uniformdens()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_snxtal_2bmb_fittedHop_uniformdens():
	outOptDict = loadDefaultTb1OptDict()
	modOptsDict = {"VxcMBCorrXtal".lower(): [4,1.0],
	               "CrystalFieldXCWeight".lower(): 0}
	modOptsDict["addCorrectingHopFromBdt".lower()] = 1
	modOptsDict["excmbcorr"] = [1, 1.0]
	outOptDict.update(modOptsDict)
	return outOptDict

@registerMethodStrToObj("tb1_snxtal_2bmb_fittedHop_uniformdens_pp_fittedPP")
def createPlatoMethod_tb1_snxtal_2bmb_fittedHop_uniformdens_pp_fittedPP():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_snxtal_2bmb_fittedHop_uniformdens_pp_fittedPP()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_snxtal_2bmb_fittedHop_uniformdens_pp_fittedPP():
	outOptDict = loadDefaultTb1OptDict()
	modOptsDict = {"VxcMBCorrXtal".lower(): [4,1.0],
	               "CrystalFieldXCWeight".lower(): 0}
	modOptsDict["addCorrectingHopFromBdt".lower()] = 1
	modOptsDict["addCorrectingPPFromBdt".lower()] = 1
	modOptsDict["excmbcorr"] = [1, 1.0]
	outOptDict.update(modOptsDict)
	return outOptDict




@registerMethodStrToObj("tb1_2bxtal_plus_sncorr_uniform_dens_pp")
def createPlatoMethod_tb1_2bxtal_plus_sncorr_uniform_dens_pp():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_2bxtal_plus_sncorr_uniform_dens_pp()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_2bxtal_plus_sncorr_uniform_dens_pp():
	outOptDict = loadDefaultTb1OptDict()
	modOptsDict = {"VxcMBCorrXtal".lower(): [3,1.0]}
	modOptsDict["excmbcorr"] = [1, 1.0]
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("tb1_2bxtal_plus_sncorr_fittedhop_uniform_dens_pp")
def createPlatoMethod_tb1_2bxtal_plus_sncorr_fittedhop_uniform_dens_pp():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_2bxtal_plus_sncorr_fittedhop_uniform_dens_pp()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_2bxtal_plus_sncorr_fittedhop_uniform_dens_pp():
	outOptDict = loadDefaultTb1OptDict()
	modOptsDict = {"VxcMBCorrXtal".lower(): [3,1.0]}
	modOptsDict["addCorrectingHopFromBdt".lower()] = 1
	modOptsDict["excmbcorr"] = [1, 1.0]
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("tb1_2bxtal_plus_sncorr_fittedhop_uniform_dens_pp_fitted_pp")
def createPlatoMethod_tb1_2bxtal_plus_sncorr_fittedhop_uniform_dens_pp_fitted_pp():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_2bxtal_plus_sncorr_fittedhop_uniform_dens_pp_fitted_pp()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_2bxtal_plus_sncorr_fittedhop_uniform_dens_pp_fitted_pp():
	outOptDict = loadDefaultTb1OptDict()
	modOptsDict = {"VxcMBCorrXtal".lower(): [3,1.0]}
	modOptsDict["addCorrectingHopFromBdt".lower()] = 1
	modOptsDict["addCorrectingPPFromBdt".lower()] = 1
	modOptsDict["excmbcorr"] = [1, 1.0]
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("tb1_vxc_sncorr3_exc_corr1")
def createPlatoMethod_tb1_vxc_sncorr4_exc_corr1():
	runCommFunct = _getRunCommPlatoTb1
	optDict = _getOptDict_tb1_vxc_sncorr3_exc_corr1()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_tb1_vxc_sncorr3_exc_corr1():
	outOptDict = loadDefaultTb1OptDict()
	modOptsDict = {"excmbcorr": [1, 1.0],
				  "VxcMBCorrXtal".lower(): [3,1.0]}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_2bhopxc_noxtalxc_pp")
@registerMethodStrToObj("dft2_2centxc_noxtal")
def createPlatoMethod_dft2_2centxc_noxtal():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2centxc_noxtal()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)


def _getOptDict_dft2_2centxc_noxtal():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 1,
	              "hopxcmethod".lower(): 1,
				  "excMbCorr".lower():0,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict

@registerMethodStrToObj("dft2_2bhopxc_xt_exact_pp_exact")
def createPlatoMethod_dft2_2bhopxc_xt_exact_pp_exact():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bhopxc_xt_exact_pp_exact()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bhopxc_xt_exact_pp_exact():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"hopxcmethod".lower(): 1}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_2bhopxc_xt_2b_pp_2b")
def createPlatoMethod_dft2_2bhopxc_xt_2b_pp_2b():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bhopxc_xt_2b_pp_2b()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bhopxc_xt_2b_pp_2b():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"hopxcmethod".lower(): 1,
	               "xtalxcmethod":2,
	               "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_mcwedahop_2bxtalxc_pp")
@registerMethodStrToObj("dft2_2centxc_withxtal")
def createPlatoMethod_dft2_2centxc_withxtal():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2centxc_withxtal()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2centxc_withxtal():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower():2,
	              "hopxcmethod".lower(): 2,
				  "excMbCorr".lower():0,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_mcwedahop_mcwedaxtal_pp")
@registerMethodStrToObj("dft2_mcweda4_pp")
def createPlatoMethod_dft2_mcweda4_pp():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_mcweda4_pp()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)


def _getOptDict_dft2_mcweda4_pp():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 4,
	              "hopxcmethod".lower(): 2,
				  "excMbCorr".lower():0,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict

@registerMethodStrToObj("dft2_mcwedahop_mcwedaxtal_pp_uniform")
@registerMethodStrToObj("dft2_mcweda4_pp_uniform")
def createPlatoMethod_dft2_mcweda4_pp_uniform():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_mcweda4_pp_uniform()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)


def _getOptDict_dft2_mcweda4_pp_uniform():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 4,
	              "hopxcmethod".lower(): 2,
				  "excMbCorr".lower():1,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict

@registerMethodStrToObj("dft2_mcwedahop_mcwedaxtal_pp_fit_uniform_scf_monopoles")
def _unused():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_mcwedahop_mcwedaxtal_pp_fit_uniform()
	optDict["scfflag"] = 1
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

@registerMethodStrToObj("dft2_mcwedahop_mcwedaxtal_pp_fit_uniform_scf_dipoles")
def _unused():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_mcwedahop_mcwedaxtal_pp_fit_uniform()
	optDict["scfflag"] = 2
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

@registerMethodStrToObj("dft2_mcwedahop_mcwedaxtal_pp_fit_uniform")
def createPlatoMethod_dft2_mcweda4_pp_uniform():
    runCommFunct = _getRunCommPlatoDft2
    optDict = _getOptDict_dft2_mcwedahop_mcwedaxtal_pp_fit_uniform()
    getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
    return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)


def _getOptDict_dft2_mcwedahop_mcwedaxtal_pp_fit_uniform():
    outOptDict = loadDefaultTb2OptDict()
    modOptsDict = {"xtalxcmethod".lower(): 4,
                  "hopxcmethod".lower(): 2,
                  "excMbCorr".lower():1,
                  "e0method":1,
	              "addCorrectingPPFromBdt".lower():1}
    outOptDict.update(modOptsDict)
    return outOptDict

@registerMethodStrToObj("dft2_mcwedahop_mcwedaxtal_pp_uniform_second_order")
@registerMethodStrToObj("dft2_mcweda4_pp_uniform_second_order")
def createPlatoMethod_dft2_mcweda4_pp_uniform_second_order():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_dft2_mcweda4_pp_uniform_second_order()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)


def _getOptDict_dft2_dft2_mcweda4_pp_uniform_second_order():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 4,
	              "hopxcmethod".lower(): 2,
				  "excMbCorr".lower():3,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_mcwedahop_mcwedaxtal_exact_e0")
@registerMethodStrToObj("dft2_mcweda4_exact_e0")
def createPlatoMethod_dft2_mcweda4_exact_e0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_mcweda4_exact_e0()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_mcweda4_exact_e0():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 4,
	               "hopxcmethod".lower(): 2,
				  "excMbCorr".lower():0,
				  "e0method":0}

	outOptDict.update(modOptsDict)
	return outOptDict

@registerMethodStrToObj("dft2_mcwedahop_xtalsncorrd_exact_e0")
@registerMethodStrToObj("dft2_mcweda5_exact_e0")
def createPlatoMethod_dft2_mcweda5_exact_e0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_mcweda5_exact_e0()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_mcweda5_exact_e0():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():5,
				  "excMbCorr".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_hop_mcweda_xtal_halfmcweda_halfsn2b_exact_e0")
def createPlatoMethod_dft2_hop_mcweda_xtal_halfmcweda_halfsn2b_exact_e0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_hop_mcweda_xtal_halfmcweda_halfsn2b_exact_e0()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_hop_mcweda_xtal_halfmcweda_halfsn2b_exact_e0():
	outOptDict = loadDefaultTb2OptDict()
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
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_mcwedahop_exactxtal_exact_e0():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 0,
	               "hopxcmethod".lower():2,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_2bodyhop_exactxtal_exacte0")
def createPlatoMethod_dft2_2bodyhop_exactxtal_exacte0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bodyhop_exactxtal_exacte0()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bodyhop_exactxtal_exacte0():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 0,
	               "hopxcmethod".lower():1,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_exacthop_noxtal_exacte0")
def createPlatoMethod_dft2_exacthop_noxtal_exacte0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exacthop_noxtal_exacte0()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exacthop_noxtal_exacte0():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 1,
	               "hopxcmethod".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_exacthop_2bxtal_exacte0")
def createPlatoMethod_dft2_exacthop_2bxtal_exacte0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exacthop_2bxtal_exacte0()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exacthop_2bxtal_exacte0():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 2,
	               "hopxcmethod".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_exacthop_2bxtal_sncorr_exacte0")
def createPlatoMethod_dft2_exacthop_2bxtal_sncorr_exacte0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exacthop_2bxtal_sncorr_exacte0()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exacthop_2bxtal_sncorr_exacte0():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 3,
	               "hopxcmethod".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_exacthop_mcweda_xtal_exacte0")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exacthop_mcweda_xtal_exacte0()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exacthop_mcweda_xtal_exacte0():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 4,
	               "hopxcmethod".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_e1_2bxc_exacte0")
@registerMethodStrToObj("dft2_2bodyhop_2bxtal_exacte0")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bodyhop_2bxtal_exacte0()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bodyhop_2bxtal_exacte0():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 2,
	               "hopxcmethod".lower():1,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_2bodyhop_mcwedaxtal_exacte0")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bodyhop_mcwedaxtal_exacte0()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bodyhop_mcwedaxtal_exacte0():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalxcmethod".lower(): 4,
	               "hopxcmethod".lower():1,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_2bodyhopvnl_else_exact")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bodyhopvnl_else_exact()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bodyhopvnl_else_exact():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"hopVnlMethod".lower():1,
				   "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_2bodyhopvna_else_exact")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bodyhopvna_else_exact()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bodyhopvna_else_exact():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"hopVnaMethod".lower():1,
				   "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_no_xtal_at_all_else_exact")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_no_xtal_at_all_else_exact()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_no_xtal_at_all_else_exact():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"xtalVnlMethod".lower():1,
	               "xtalVnaMethod".lower():1,
	               "xtalXcMethod".lower():1,
				   "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_2bodyhopvna_2bodyhopvnl_else_exact")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bodyhopvna_2bodyhopvnl_else_exact()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bodyhopvna_2bodyhopvnl_else_exact():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"hopVnlMethod".lower():1,
	               "hopVnaMethod".lower():1,
				   "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_2bodyhop_all_else_exact")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bodyhop_all_else_exact()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bodyhop_all_else_exact():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"hopVnlMethod".lower():1,
	               "hopVnaMethod".lower():1,
	               "hopXcMethod".lower():1,
				   "e0method".lower():0}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_2bodyhop_2bodyxtal_pp")
def createPlatoMethod_():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2bodyhop_2bodyxtal_pp()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2bodyhop_2bodyxtal_pp():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"hopVnlMethod".lower():1,
	               "hopVnaMethod".lower():1,
	               "hopXcMethod".lower():1,
				   "e0method".lower():1}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_exact_e1_pp")
def createPlatoMethod_dft2_exact():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exact_e1_pp()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exact_e1_pp():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():0,
				  "excMbCorr".lower():0,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_2body_e1_exact_e0")
def createPlatoMethod_dft2_exact():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_2body_e1_exact_e0()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_2body_e1_exact_e0():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():3,
				  "excMbCorr".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_exact_e1_pp_uniform_e0")
def createPlatoMethod_dft2_exact_e1_pp_uniform_e0():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exact_e1_pp_uniform_e0()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exact_e1_pp_uniform_e0():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():0,
				  "excMbCorr".lower():1,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_exact_e1_pp_uniform_second_order")
def createPlatoMethod_dft2_exact_e1_pp_uniform_second_order():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exact_e1_pp_uniform_second_order()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exact_e1_pp_uniform_second_order():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"excMbCorr".lower():3,
				  "e0method":1}
	outOptDict.update(modOptsDict)
	return outOptDict


#@registerMethodStrToObj("dft2_minimal_for_debuggingA")
#def createPlatoMethod_dft2_minimal_for_debuggingA():
#	runCommFunct = _getRunCommPlatoDft2
#	optDict = _getOptDict_dft2_minimal_for_debuggingA()
#	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
#	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)
#
#def _getOptDict_dft2_minimal_for_debuggingA():
#	outOptDict = loadDefaultTb2OptDict()
#	modOptsDict = {"XtalXCMethod".lower():-1,
#				  "HopXCMethod".lower():-1,
#	              "XtalVnaMethod".lower():-1,
#	              "HopVnaMethod".lower():-1,
#	              "XtalVnlMethod".lower():-1,
#	              "HopVnlMethod".lower():-1}
#	outOptDict.update(modOptsDict)
#	return outOptDict


@registerMethodStrToObj("dft2_exact_scf_monopoles")
def _():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exact()
	optDict["scfflag"] = 1
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

@registerMethodStrToObj("dft2_exact_scf_dipoles")
def _():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exact()
	optDict["scfflag"] = 2
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)


@registerMethodStrToObj("dft2_exact")
def createPlatoMethod_dft2_exact():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_exact()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_exact():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():0,
				  "excMbCorr".lower():0,
				  "e0method":0}
	outOptDict.update(modOptsDict)
	return outOptDict


@registerMethodStrToObj("dft2_lda")
def createPlatoMethod_dft2_lda():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_lda()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_lda():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():0,
				  "excMbCorr".lower():0,
				  "e0method":0,
	              "XCFunctional".lower():"lda"}
	outOptDict.update(modOptsDict)
	return outOptDict



@registerMethodStrToObj("dft2_lda_1body_pbe_corr")
def createPlatoMethod_dft2_lda_1body_pbe_corr():
	runCommFunct = _getRunCommPlatoDft2
	optDict = _getOptDict_dft2_lda_1body_pbe_corr()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_tb1OrTb2
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict)

def _getOptDict_dft2_lda_1body_pbe_corr():
	outOptDict = loadDefaultTb2OptDict()
	modOptsDict = {"McWedaXcFlag".lower():0,
				  "excMbCorr".lower():0,
				  "e0method":0,
	              "XCFunctional".lower():"lda",
	              "E0NonLocalXcCorr".lower(): 1,
	              "e1XtalNonLocalXcCorr".lower():1,
	              "e1HopNonLocalXcCorr".lower():1}
	outOptDict.update(modOptsDict)
	return outOptDict




@registerMethodStrToObj("dft_plato")
def createPlatoMethod_dft_plato():
	runCommFunct = _getRunCommPlatoDft
	optDict = _getOptDict_dft_plato()
	getStrDictFromOptDict = getPlatoStrDictFromOptDict_dft
	gridKwarg = "fftGridSpacing".lower()
	return PlatoMethod(optDict, runCommFunct, getStrDictFromOptDict, gridKwarg=gridKwarg)


def _getOptDict_dft_plato():
	outOptDict = loadDefaultDftOptDict()
	modOptDict = {"xcfunctional": "pbe",
	              "scfflag": 1,
	              "optimizemesh": 0} #Some of these MIGHT be defaults
	outOptDict.update(modOptDict)
	return outOptDict






