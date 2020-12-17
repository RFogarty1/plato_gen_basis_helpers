
from ..shared import method_objs as methObjsHelp
from . import cp2k_basis_obj as basisObjHelp

def getSimpleCreatorObjToDictMapObj():
	outFuncts = list()
	outFuncts.append(getSelectedBasicInfoDictFromCreatorObj)
	outFuncts.append(_getSelectedDictsFromCreatorObj)
	outFuncts.append(_getGeomConstraintInfoFromCreatorObj)
	outFuncts.append(_getBasisObjsInfoFromCreatorObj)
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
	attrsToGet = ["grimmeDisp", "nonLocalDisp", "surfaceDipoleCorr"]
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


def _getBasisObjsInfoFromCreatorObj(creatorObj):
	outDict = dict()
	if creatorObj.basisObjs is None:
		pass
	else:
		outDict["basisObjs"] = getDictFromBasisObjList( creatorObj.basisObjs )

	return outDict


def getDictFromBasisObjList(basisObjs):
	""" Gets a convenient (in terms of searching) dictionary representation of an iter of CP2KBasisObjStandard. Purpose is to dump basis set information to the database
	
	Args:
		basisObjs: (iter of CP2KBasisObjStandard) 
			 
	Returns
		 outDict: (dict). Keys are .kind values of basisObjs Values are dicts containing all info on the basis set assigned to that kind
 
	"""
	outDict = dict()
	for obj in basisObjs:
		currKey = obj.kind
		currDict = obj.toDict()
		outDict[currKey] = currDict
	return outDict

def getBasisObjListFromDict(basisObjDict):
	""" Gets a list of CP2KBasisObjStandard objects from a dictionary representation used for the database
	
	Args:
		basisObjDict: (dict) This is whats returned from getDictFromBasisObjList; it should basically never be generated any other way
			 
	Returns
		basisObjList: (list of CP2KBasisObjStandard objs) Each object represents the basis set information for one "kind" (usually element, but can have something like Mg_surface to use diff basis for diff "kinds" of atom which are the same element)
 
	"""
	outList = list()
	for key in basisObjDict.keys():
		currObj = basisObjHelp.CP2KBasisObjStandard.fromDict(basisObjDict[key])
		outList.append(currObj)
	return outList

