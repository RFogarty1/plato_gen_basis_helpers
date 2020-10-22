

def getSimpleCreatorObjToDictMapObj():
	outFuncts = list()
	outFuncts.append(getSelectedBasicInfoDictFromCreatorObj)
	outFuncts.append(_getSelectedDictsFromCreatorObj)
	outFuncts.append(_getGeomConstraintInfoFromCreatorObj)
	return GetDictFromCreatorObj(outFuncts)

class GetDictFromCreatorObj():
	""" Callable class for getting a dictionary object from a BaseCP2KCalcObjFactory object. Callable interface is f(BaseCP2KCalcObjFactory)->dict

	"""

	def __init__(self, creatorToDictFuncts, modFinalDictFuncts=None):
		""" Initialiser
		
		Args:
			creatorToDictFuncts: (iter of functions) Each function takes a BaseCP2KCalcObjFactory instance and returns a dictionary
			modFinalDictFuncts: (iter of functions, Optional) Each function modifies a dictionary in place

		Important notes (for args):
			The creatorToDictFuncts are each called first IN ORDER. The output dict updates with EACH of the dicts produced by these functions, so if two return the same key then the second will overwrite the first. The modFinalDictFuncts are run AFTER the creatorToDictFuncts, so they can do things lke remove certain keys if required.
				 
		"""
		self.creatorToDictFuncts = list(creatorToDictFuncts)
		self.modFinalDicts = list() if modFinalDictFuncts is None else list(modFinalDictFuncts)


	def _getOutputDict(self, creator):
		outDict = dict()
		for x in self.creatorToDictFuncts:
			outDict.update( x(creator) )

		for x in self.modFinalDicts:
			x(outDict)

		return outDict

	def __call__(self, creator):
		return self._getOutputDict(creator)

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







