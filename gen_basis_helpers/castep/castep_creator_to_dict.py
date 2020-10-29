from ..shared import method_objs as methObjsHelp

def getSimpleCreatorObjToDictMapObj():
	outFuncts = list()
	outFuncts.append(getSelectedBasicInfoDictFromCreatorObj)
	return methObjsHelp.GetDictFromCreatorObj(outFuncts)


def getSelectedBasicInfoDictFromCreatorObj(creatorObj):
	outDict = dict()
	attrsToGet = ["charge", "cutoffEnergy", "kPts"]
	for currAttr in attrsToGet:
		currVal = getattr(creatorObj, currAttr)
		if currVal is not None:
			outDict[currAttr] = currVal
	return outDict


