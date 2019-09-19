
import copy
import os

import plato_pylib.plato.mod_plato_inp_files as modInp
import plato_fit_integrals.initialise.create_ecurve_workflows as eCurveWorkFlows

from . import convergers as convergers
from . import integ_grid_conv as gridConv

''' Module is to help converging the total energy of single structures '''


POSSIBLE_CONV_TYPES = ["angularGridConv","fftGridConv", "kPtConv", "radialGridConv"]


class ElementConfigEqmStructs():
	"""Static class (no args in init) used to get default grid meshes/k points etc. for running convergence.

	User is expected to subclass this and override the relevant methods

	"""

	def getConvVariableValues(self, structType, varyType):
		''' Get a list of parameters which we vary (e.g. the k-points we use to look at convergence)'''
		raise NotImplementedError
	
	def getDefaultMeshSpacing(self, structStr, varyType):
		''' Get the default mesh parameters '''
		raise NotImplementedError
		
	def getMeshSpacingKwarg(self, varyType):
		''' Get the mesh keyword for varytype. Used to differentiate between grid and code types '''
		if (varyType == "fftGridConv"):
			return "fftGridSpacing".lower()
		else:
			return "IntegralMeshSpacing".lower()
		
	def getPlatoCodeFromVaryType(self,varyType):
		if varyType=="fftGridConv":
			return "dft"
		else:
			return "dft2"

	def getDefKPoints(self,structStr):
		raise NotImplementedError



#TODO: Implement eleAndStructKeyToUCell
class TotalEnergyConvRunnerOptions():
	"""Class which encapsulates all the options (except converged paramter type) for running convergence calculations 
	for total energy for a single structure. NOTE: Recommended to use the fromEnforcedKwargs initialisation function

	Attributes:
		eleConfigStructs (dict): Keys are elements(lower case), values are ElementConfigEqmStructs objects
		elements (iter of str): list of the keys for eleConfigStructs
		structKeys (iter of str): list of keys which map to structures
		eleAndStructKeyToUCell (function): Function which accepts element and structKey, and returns a UnitCell object. Callable class recommended
		refDataObj (SimpleNamespace): Attrs are element symbols (element.capitalize()), values are RefElementalDataBase objects
		baseFolder (str): base workFolder. Recommened os.path.join("work_folder", "convergence_calcs")
		runJobsDict (dict): keys are convergence types, [see VALID_VARY_TYPES] while values are True/False (determines whether to run jobs for those varyTypes)
		eTypeAll (str): String denoting the type of energy to use
	"""

	def __init__(self, **kwargs):
		self.eleConfigStructs = kwargs.get( "eleConfigStructs", None )
		self.elements = kwargs.get( "elements", None )
		self.structKeys = kwargs.get( "structKeys", None )
		self.eleAndStructKeyToUCell = kwargs.get( "eleAndStructKeyToUCell", None )
		self.refDataObj = kwargs.get( "refDataObj", None ) #SimpleNamespace, keys are elements. Values are refEleObjs.RefElementalDataBase objects
		self.baseFolder = kwargs.get( "baseFolder", None )
		self.runJobsDict = kwargs.get( "runJobsDict", None ) #Keys are convType keywords, values are True/False for whether to run the jobs or not
		self.eTypeAll = kwargs.get( "eTypeAll", None ) #String denoting the type of energy to pull from files
		
	@classmethod
	def fromEnforcedKwargs(cls,**kwargs):
		enforcedKwargs = ("eleConfigStructs", "elements", "structKeys", "refDataObj", "baseFolder",
						  "runJobsDict", "eTypeAll", "eleAndStructKeyToUCell")
		for key in enforcedKwargs:
			assert key in kwargs.keys(), "{} is a required keyword argument".format(key)
		return cls(**kwargs)





#THIS is the main funciton that should be called
def createJobRunnersOneConvType(convRunnerOpts:"TotalEnergyConvRunnerOptions obj", varyType, labelExt=None):
        allRunners = list()

        if labelExt is None:
            labelExt = ""
            
        #Create the actual individual runners
        for ele in convRunnerOpts.elements:
            currRunner = _createJobRunnerSingleElementAndConvType(ele, convRunnerOpts,varyType,labelExt)
            allRunners.append(currRunner)
        outRunner = convergers.PropConvJobRunComposite(allRunners)
        
        return outRunner
        
def _createJobRunnerSingleElementAndConvType(element, convRunnerOpts, varyType, labelExt):
    allRunners = list()
    for struct in convRunnerOpts.structKeys:
        currRunner = _createJobRunnerSingleElementConvStructType(element, convRunnerOpts, struct, varyType, labelExt)
        allRunners.append(currRunner)
    outRunner = convergers.PropConvJobRunComposite(allRunners)
    return outRunner




def _createJobRunnerSingleElementConvStructType(element, convRunnerOpts, structKey,varyType, labelExt):
    paramToOptDictMapper = _getVaryTypeParamToOptDictMapper(varyType)
    uCell = convRunnerOpts.eleAndStructKeyToUCell(element,structKey)
    label = "{}_{}_{}{}".format(element,varyType,structKey,labelExt)
    workFolder = os.path.join(convRunnerOpts.baseFolder,element,structKey,varyType)
    platoCode = convRunnerOpts.eleConfigStructs[element].getPlatoCodeFromVaryType(varyType)
    
    baseModOptsDict = _getBaseModOptsDictFromVaryTypeEleConfigAndStructKey(varyType, convRunnerOpts.eleConfigStructs[element], structKey)
    
    baseModOptsDict["dataset"] = _getDataSetFromPlatoCodeAndElementAndRefDataStruct(platoCode,element,convRunnerOpts.refDataObj)
    varyParams = convRunnerOpts.eleConfigStructs[element].getConvVariableValues(structKey,varyType)
    
    varyParamsSingleVals = _getSingleValRepForVaryParams(varyParams)
    
    #Create each workflow
    allWorkFlows = list()
    for vParam,vSingleVal in zip(varyParams,varyParamsSingleVals):
        modOptsDict = copy.deepcopy(baseModOptsDict)
        paramToOptDictMapper.modOptDictWithVariableParam(modOptsDict,vParam)
        
        varyParamStr = "{:.3f}".format(vSingleVal).replace(".","pt")
        currWorkFolder = os.path.join(workFolder, varyParamStr )
        outWorkFlow = eCurveWorkFlows.CreateStructEnergiesWorkFlow([uCell], modOptsDict, currWorkFolder, platoCode,
                                                                     varyType=None, eType=convRunnerOpts.eTypeAll)()
        _modECurveWorkFlowToReturnFirstListValue(outWorkFlow)
        _modECurveWorkFlowToReturnValueInEv(outWorkFlow)
        allWorkFlows.append(outWorkFlow)
        
    shouldWeRunCalcs = convRunnerOpts.runJobsDict[varyType]
    outRunner = convergers.PropConvJobRunnerStandard(allWorkFlows, varyParamsSingleVals, label,shouldWeRunCalcs=shouldWeRunCalcs)
    
    return outRunner

def _getSingleValRepForVaryParams(varyParams):
    ''' We need our params that we vary to be a single number, but sometimes (e.g k points) the reality 
    is that we need to originally represent them as an iterable. This handles that by taking the first value
    of an iterable (if it IS iterable), and else just using the original value
    '''
    firstVal = varyParams[0]
    #Should return if individual params arent iterables
    try:
        iter(firstVal)
    except TypeError:
        return list(varyParams)
    
    #Should only run for iterables
    return [x[0] for x in varyParams]
    
def _getVaryTypeParamToOptDictMapper(varyType):
    if varyType == "angularGridConv":
        outVal = gridConv.createDft2AngularSpacingOptDictMapper()
    elif varyType == "fftGridConv":
        outVal = gridConv.createDftFFTSpacingOptDictMapper()
    elif varyType == "kPtConv":
        outVal = gridConv.createKPointsMPGridOptDictMapper()
    elif varyType == "radialGridConv":
        outVal = gridConv.createDft2RadialSpacingOptDictMapper()
    else:
        raise ValueError("varyType={} is wrong".format(varyType))
    return outVal
        

def _getBaseModOptsDictFromVaryTypeEleConfigAndStructKey(varyType,eleConfig,structKey):
    platoCode = eleConfig.getPlatoCodeFromVaryType(varyType)
    kPoints = eleConfig.getDefKPoints(structKey)
    gridParams = eleConfig.getDefaultMeshSpacing(structKey, varyType)
    gridKwarg = eleConfig.getMeshSpacingKwarg(varyType)
    modOptsDict = modInp.getDefOptDict(platoCode) #Need full dict or error thrown by paramToOptDictMapper 
    modOptsDict["blochstates"] = kPoints
    modOptsDict[gridKwarg] = gridParams
    return modOptsDict


def _getDataSetFromPlatoCodeAndElementAndRefDataStruct(platoCode, element,refDataObj):
    eleData = getattr(refDataObj,element.capitalize())
    if platoCode == "dft2":
        return eleData.modelFiles.dft2PlatoPath
    elif platoCode == "dft":
        return eleData.modelFiles.dftPlatoPath #TODO: modelFiles is bugged for this; fix it
    else:
        raise ValueError("platoCode = {} isnt allowed".format(platoCode))

        
def _modECurveWorkFlowToReturnFirstListValue(workFlow):
    f = workFlow.run
    def runThenCollapseListIntoVal():
        f()
        assert len(workFlow.namespaceAttrs) == 1
        setList = getattr(workFlow.output, workFlow.namespaceAttrs[0])
        singleValue = setList[0]
        setattr(workFlow.output, workFlow.namespaceAttrs[0], singleValue)
        return None
    workFlow.run = runThenCollapseListIntoVal
    return None

#NOTE: assumed we've already applied _modECurveWorkFlowToReturnFirstListValue
def _modECurveWorkFlowToReturnValueInEv(workFlow):
    f = workFlow.run
    rydToEV = 13.6056980659
    def convRunOutputFromRydToEv():
        f()
        assert len(workFlow.namespaceAttrs) == 1
        rydVal = getattr(workFlow.output, workFlow.namespaceAttrs[0])
        evVal = rydVal*rydToEV
        setattr(workFlow.output, workFlow.namespaceAttrs[0], evVal)
    workFlow.run = convRunOutputFromRydToEv

