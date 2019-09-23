

from . import multi_cryst_eos as multCrystEos


def createSingleCrystEosFromWorkFlowCoord(wFlowCoord):
	allSingleCrysts = list()
	for wFlow in wFlowCoord._workFlows:
		currStruct = wFlow.structKey
		currElement = wFlow.element
		currFitDict = wFlow.extraOutput.full_eos
		currSingleCryst = multCrystEos.SingleCrystEosResult(v0=currFitDict["v0"], b0=currFitDict["b0"], e0=currFitDict["e0"], fitData=currFitDict["fitdata"],
											   actData = currFitDict["data"], structLabel = currStruct, elementLabel=currElement,
											  methodLabel=wFlow.methodKey)
		allSingleCrysts.append(currSingleCryst)
	return allSingleCrysts

def getMultiCrystEosResultsFromSingleCrysts(allSingleCrysts):
	singleWorkFlowsGroupedByMethodAndElement = _groupSingleCrystEosByMethodAndElement(allSingleCrysts)
	multiEosGrouped = _groupSingleCrystEosIntoMultiCrysts(singleWorkFlowsGroupedByMethodAndElement)
	outList = _groupMultiCrystsByElement(multiEosGrouped)
	return outList


def getMultiCrystEosResultsFromWorkFlowCoord(wFlowCoord):
	#Loop over the individual workflows directly. Use the labels i inserted into them
	singleWorkFlows = _createSingleCrystEosFromWorkFlowCoord(wFlowCoord)

	

def _groupSingleCrystEosByMethodAndElement(singleCrystEosList):
	allMethods = set( [x.methodLabel for x in singleCrystEosList] )
	allElements = set( [x.elementLabel for x in singleCrystEosList] )
	groupedObjs = list()
	print("allMethods = {}".format(allMethods))
	print("allElements = {}".format(allElements))
	for mKey in allMethods:
		for eleKey in allElements:
			currList = list()
			for currObj in singleCrystEosList:
				if (currObj.methodLabel == mKey) and (currObj.elementLabel == eleKey):
					currList.append( currObj )
			groupedObjs.append(currList)
		
	return groupedObjs
	
def _groupSingleCrystEosIntoMultiCrysts(singleCrystEosList):
	multiCrystObjs = list()
	for crystSet in singleCrystEosList:
		currSet = multCrystEos.MultiCrystEosResult( crystSet, crystSet[0].methodLabel)
		multiCrystObjs.append(currSet)
	
	return multiCrystObjs

def _groupMultiCrystsByElement(multiCrystsList):
	#Step 1 = get ALL ELEMENTS. TODO: ACTUALLY FINISH THIS
	allElements = list()
	for x in multiCrystsList:
		allElements.append(x.elementLabel)
	allElements = list(set(allElements))
	
	#Step 2 = group them into lists
	allObjInpLists = list()
	for currEle in allElements:
		currList = list()
		for x in multiCrystsList:
			if x.elementLabel == currEle:
				currList.append(x)
		allObjInpLists.append(currList)

	#Step 3 = create objects from the grouped lists
	allObjs = [multCrystEos.GroupedMultiCrystEosForOneElement(x) for x in allObjInpLists]
	
	return allObjs



def getPWSingleCrystEosForRelevantElementsAndStructs(elements, structKeys:list, refDataAllElements, eos):
	allDict = dict()
	allDict["mg"] = _getPWSingleCrystEosMg(structKeys, refDataAllElements, eos)
	allDict["zr"] = _getPWSingleCrystEosZr(structKeys, refDataAllElements, eos)
	outList = list()
	for x in allDict.keys():
		if x in elements:
			outList.extend(allDict[x])
	return outList
	
	
def _getPWSingleCrystEosMg(structKeys,refDataAllElements, eos):
	infoStruct = refDataAllElements.Mg
	elementLabel = "mg"
	return _getPlaneWaveObjListFromInfoStructAndElementLabel(structKeys,infoStruct,elementLabel, eos)

def _getPWSingleCrystEosZr(structKeys,refDataAllElements, eos):
	infoStruct = refDataAllElements.Zr
	elementLabel = "zr"
	return _getPlaneWaveObjListFromInfoStructAndElementLabel(structKeys,infoStruct,elementLabel, eos)


def _getPlaneWaveObjListFromInfoStructAndElementLabel(structKeys,infoStruct,elementLabel, eos):
	methodLabel = "plane-wave"
	outList = list()
	for structKey in structKeys:
		eosDict = infoStruct.getEosFitDict(structKey,eos=eos)
		v0, b0, e0 = eosDict["v0"], eosDict["b0"], eosDict["e0"]
		fitData, actData = eosDict["fitdata"], eosDict["data"]
		currObj = multCrystEos.SingleCrystEosResult(v0=v0, e0=e0, b0=b0, fitData=fitData, actData=actData,
									   structLabel=structKey, elementLabel=elementLabel, methodLabel=methodLabel)
		outList.append(currObj)
	return outList

