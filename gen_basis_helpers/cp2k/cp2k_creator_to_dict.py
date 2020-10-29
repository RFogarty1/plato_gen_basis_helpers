
from ..shared import method_objs as methObjsHelp

def getSimpleCreatorObjToDictMapObj():
	outFuncts = list()
	outFuncts.append(getSelectedBasicInfoDictFromCreatorObj)
	outFuncts.append(_getSelectedDictsFromCreatorObj)
	outFuncts.append(_getGeomConstraintInfoFromCreatorObj)
	return methObjsHelp.GetDictFromCreatorObj(outFuncts)


def getSelectedBasicInfoDictFromCreatorObj(creatorObj):
	outDict = dict()
	attrsToGet = ["absGridCutoff", "addedMOs", "charge", "kPts","relGridCutoff", "xcFunctional"]
	for currAttr in attrsToGet:
		currVal = getattr(creatorObj, currAttr)
		if currVal is not None:
			outDict[currAttr] = currVal
	return outDict

def _getSelectedDictsFromCreatorObj(creatorObj):
	outDict = dict()
	attrsToGet = ["grimmeDisp", "nonLocalDisp"]
	for currAttr in attrsToGet:
		currVal = getattr(creatorObj, currAttr)
		if currVal is not None:
			outDict[currAttr] = currVal.toDict()
	return outDict

def _getGeomConstraintInfoFromCreatorObj(creatorObj):
	outDict = dict()
	if creatorObj.geomConstraints is None:
		pass
	else:
		outDict["cell_constraints"] = creatorObj.geomConstraints.cellConstraints.toDict()
	return outDict







