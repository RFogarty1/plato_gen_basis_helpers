
import copy
import os
from collections import OrderedDict

import plato_fit_integrals.initialise.create_eos_workflows as createEosWorkflow
import plato_fit_integrals.core.workflow_coordinator as wFlowCoordinator

from ..shared import calc_methods as calcMethods




def createElementInfoStructs(elements, structKeys, refStructDataPureElements, refConvDataPureElements, labelExt=None):
	#Step 1 = create ALL possible elemental structs
	outDict = OrderedDict()
	for eleKey in elements:
		outDict[eleKey] = _createElementalInfoStructForOneKey(eleKey, refStructDataPureElements, refConvDataPureElements, structKeys)
	return outDict


def _createElementalInfoStructForOneKey(eleKey, refStructDataPureElements, refConvDataPureElements, structKeys):
	relStructData = getattr(refStructDataPureElements,eleKey.capitalize())
	relConvData = getattr(refConvDataPureElements, eleKey.capitalize())
	modelFolders = copy.deepcopy( getattr(relStructData, "modelFiles") )
	structDict = _getStructDictFromReferenceElementAndConvData(relStructData, relConvData, structKeys)
	return ElementalInfo(modelFolders=modelFolders, structDict=structDict, label=eleKey.lower())

	
def _getStructDictFromReferenceElementAndConvData(refEleData, refConvData, structKeys):
	outDict = OrderedDict()
	for structKey in structKeys:
		currStructs = refEleData.getStructsForEos(structKey)
		currKptsInfo = refConvData.kptGridVals.getKptsPrimCell(structKey)
		dft2GridInfo = refConvData.integGridVals.getPrimCellDft2AngularGrid(structKey)
		fftGridInfo = refConvData.integGridVals.getPrimCellDftFFTGrid(structKey)
		outDict[structKey] = StructSetInfo(geoms=currStructs, kPts=currKptsInfo,
										   dft2Grid=dft2GridInfo, fftGrid=fftGridInfo)
	return outDict





class ElementalInfo():
	""" Holds Info that varies between elements tested. Such as the model data folder"""
	def __init__(self,modelFolders=None, structDict=None, label=None):
		self.modelFolders = modelFolders
		self.structDict = structDict #Contains StructSetInfo objects
		self.label = label #The actual element

		reqArgs = [modelFolders,structDict,label]
		if any([x is None for x in reqArgs]):
			raise ValueError("Required argument missing")
		
		
class StructSetInfo():
	""" Holds calculation info thats relevent to a single set of structures (e.g. bcc)"""
	def __init__(self,geoms=None, kPts=None, dft2Grid=None, fftGrid=None):
		self.kPts = kPts #1 value; assumed same k-points for a set of structures here
		self.dft2GridVals = dft2Grid
		self.fftGridVals = fftGrid
		self.geoms = geoms #geometries
		reqArgs = [geoms,kPts,dft2Grid,fftGrid]
		if any([x is None for x in reqArgs]):
			raise ValueError("One of the arguments was not set")















def createWorkFlowCoordinator(eleInfoStructs, methodStrs, eos, baseFolder, nCores):
	allWorkFlows = list()
	for eleStr in eleInfoStructs.keys():
		currEleInfo = eleInfoStructs[eleStr]
		for methodStr in methodStrs:
			currWorkFlows = _createWorkFlowsFromElementInfoAndMethodStr(currEleInfo, methodStr, eos, baseFolder)
			allWorkFlows.extend(currWorkFlows)
		
	wFlowCoord = wFlowCoordinator.WorkFlowCoordinator(allWorkFlows, nCores=nCores)
	wFlowCoord.quietPreShellComms = False
	return wFlowCoord


def _createWorkFlowsFromElementInfoAndMethodStr(elementInfo,methodStr, eos, baseFolder):
	allWorkFlows = list()
	for structKey in elementInfo.structDict.keys():
		currFactory = _createWorkFlowFactoryFromElementInfoAndStructKeyAndMethodStr(elementInfo,structKey, methodStr, eos, baseFolder)
		allWorkFlows.append( currFactory() )
	return allWorkFlows

def _createWorkFlowFactoryFromElementInfoAndStructKeyAndMethodStr(elementInfo, structKey, methodStr, eos, baseFolder):
	currStructSetInfo = elementInfo.structDict[structKey]
	#Sort grid conv based on the method str
	if methodStr.startswith("dft2_"):
		gridConv = currStructSetInfo.dft2GridVals
		modelFolder = elementInfo.modelFolders.dft2PlatoPath
	elif methodStr.startswith("dft_"):
		gridConv = currStructSetInfo.fftGridVals
		modelFolder = elementInfo.modelFolders.dftModel
	elif methodStr.startswith("tb1_"):
		gridConv = currStructSetInfo.dft2GridVals
		modelFolder = elementInfo.modelFolders.tb1PlatoPath
	else:
		raise ValueError("Cant determine grid spacing for methodStr {}".format(methodStr))
	
	geoms = currStructSetInfo.geoms
	kPts = currStructSetInfo.kPts
	workFolder = os.path.join(baseFolder, elementInfo.label,structKey)
	
	return WorkFlowFactory(modelFolder=modelFolder, methodStr=methodStr, kPts=kPts, gridConv=gridConv,
						   structName=structKey, workFolder=workFolder, structList=geoms,
						   element=elementInfo.label, eos=eos)





class WorkFlowFactory():

	def __init__(self, modelFolder=None, methodStr=None, kPts=None, gridConv=None, structName=None,
				workFolder=None, structList=None, element=None, eos=None):
		self.modelFolder = modelFolder
		self.methodStr = methodStr
		self.kPts = kPts
		self.gridConv = gridConv
		self.structName = structName #e.g. bcc/fcc/hcp
		self.workFolder = workFolder
		self.structList = structList #List of u-cell structures
		self.element = element
		self.eos = eos
		
		reqArgs = [modelFolder, methodStr, kPts, gridConv, structName, workFolder, structList, element, eos]
		if any([x is None for x in reqArgs]):
			raise ValueError("Required argument is missing")
		
	@property
	def _platoCode(self):
		if self.methodStr.startswith("dft2_"):
			return "dft2"
		elif self.methodStr.startswith("dft_"):
			return "dft"
		elif self.methodStr.startswith("tb1_"):
			return "tb1"
		else:
			raise ValueError("self.methodStr = {} is an invalid value".format(self.methodStr))
	
	def __call__(self):
		uniqueKey = self.element + "_" + self.structName + "_" + self.methodStr
		structDict = {uniqueKey:self.structList}
		workFolder = os.path.join(self.workFolder,self.methodStr)
		
		#Get the opt-dict
		platoMethod = calcMethods.createPlatoMethodObj(self.methodStr)
		platoMethod.kpts = self.kPts
		platoMethod.dataSet = self.modelFolder
		platoMethod.integGrid = self.gridConv
		baseOptDict = platoMethod.optDict
		outDict = {uniqueKey:baseOptDict}
		
		#TODO: THIS BIT OF CODE - finish making the factory object
		
		factoryInstance = createEosWorkflow.CreateEosWorkFlow(structDict, outDict, workFolder,
															  self._platoCode, eosModel=self.eos, varyType=None)
		objInstance = factoryInstance()
		objInstance.element = self.element
		objInstance.structKey = self.structName
		objInstance.methodKey = self.methodStr
		
		return objInstance

