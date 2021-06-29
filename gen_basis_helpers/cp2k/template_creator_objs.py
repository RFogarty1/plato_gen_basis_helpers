
import plato_pylib.shared.unit_convs as uConvHelp

from . import collective_vars as colVarHelp
from . import cp2k_creator as cp2kCreatorHelp
from . import cp2k_misc_objs as cp2kMiscObjHelp
from . import cp2k_md_options as cp2kMdOptsHelp
from ..shared import register_key_decorator as regKeyDeco


_CREATOR_MOD_OBJ_DICT = dict()


registerCreatorModDeco = regKeyDeco.RegisterKeyValDecorator(_CREATOR_MOD_OBJ_DICT, forceKeysToCase="lower")



#Utility functions here
def getRegisteredTemplateCreatorObjStrings():
	return _CREATOR_MOD_OBJ_DICT.keys()

#Maybe have an option to pass in a list; apply modifications in turn
def getTemplateCreatorFromString(inpStr):
	""" Gets a CP2KCalcObjFactoryStandard object with various options preset from a string
	
	Args:
		inpStr: (Str or iter of str, case-insensitive) See getRegisteredTemplateCreatorObjStrings for options. If an iter of strings, then templates will be applied in order to a blank CP2KCalcObjFactoryStandard object (later can override earlier)
			 
	Returns
		outObj: (CP2KCalcObjFactoryStandard object) With some options preset
 
	"""
	outObj = cp2kCreatorHelp.CP2KCalcObjFactoryStandard()
	if isinstance(inpStr,str):
		_CREATOR_MOD_OBJ_DICT[inpStr.lower()](outObj)
	else:
		for currStr in inpStr:
			_CREATOR_MOD_OBJ_DICT[currStr.lower()](outObj)
	return outObj




#-------------------->Standard mod functions below<-------------------------

#This ones so that we have some stuff set; mainly to remind me what needs chaning 
@registerCreatorModDeco("simple_generic")
def _unused(inpObj):
	kwargDict = {"absGridCutoff":800,"addedMOs":10, "epsDef":1e-10, "epsScf":1e-6,
	             "kPts":"None", "methodStr":"spe_standard", "relGridCutoff":500,
	             "xcFunctional":"PBE", "useSmearing":True, "walltime":3600*71}
	inpObj._updateAttrsFromKwargs(**kwargDict)

@registerCreatorModDeco("metal_ot_spe")
def _unused(inpObj):
	kwargDict = {"addedMOs":0, "epsDef":1e-12, "epsScf":1e-6, "kPts":"None",
	             "methodStr":"spe_standard", "scfMaxIters":20, "scfOuterMaxIters":20,
	             "scfDiagOn":False, "scfOTEnergyGap":1e-4, "scfOTMinimizer":"diis",
	             "scfOTPreconditioner":"FULL_SINGLE_INVERSE", "scfOuterEps":1e-6,
	             "scfOTStepsize":0.02, "useSmearing":False, "walltime":3600*71}
	inpObj._updateAttrsFromKwargs(**kwargDict)


@registerCreatorModDeco("optb88-vdw-drsll")
def _unused(inpObj):
	dispObj = cp2kMiscObjHelp.NonLocalDispersionsCorrOptsCP2K(corrType="DRSLL")
	kwargDict = {"nonLocalDisp":dispObj,"xcFunctional":"optb88_PW92"}
	inpObj._updateAttrsFromKwargs(**kwargDict)


@registerCreatorModDeco("high_frict_langevin_md")
def _unused(inpObj):
	currKwargs = {"timeStep": 0.5, "nSteps": 50000, "ensemble": "Langevin",
	              "temperature": 300, "printKindTemp": True}
	stdMdOpts = cp2kMdOptsHelp.MolDynamicsOptsCP2KStandard(**currKwargs)
	thermostatOpts = cp2kMdOptsHelp.LangevinThermostatOpts(gamma=0.01)

	outKwargDict = {"mdOpts":stdMdOpts, "thermostatOpts":thermostatOpts, "runType":"md"}
	inpObj._updateAttrsFromKwargs(**outKwargDict)	


@registerCreatorModDeco("metadyn_aggressive_single_colvar")
def _unused(inpObj):
	metaVar = colVarHelp.MetaVarStandard()
	evToHa = 0.5*uConvHelp.EV_TO_RYD
	currKwargs = {"ntHills":20, "hillHeight":0.2*evToHa}
	metaObj = cp2kMdOptsHelp.MetaDynamicsOptsCP2KStandard([metaVar], **currKwargs)
	kwargDict = {"metaDynOpts":metaObj}
	inpObj._updateAttrsFromKwargs(**kwargDict)

@registerCreatorModDeco("hirshfeld_iterative_standard")
def _unused(inpObj):
	kwargDict = {"hirshfeld_selfConsistent":True, "hirshfeld_shapeFunction":"density"}
	inpObj._updateAttrsFromKwargs(**kwargDict)




