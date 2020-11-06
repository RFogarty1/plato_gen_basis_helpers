
import copy

def createFunctToGetStandardInpCreatorFromOverideDicts(templateObjCreator, templateStdInpCreator):
	""" Gets a function for getting standard creator objects from overide dictionaries. See CallableGetStandardInpCreatorFromOverideDictsStandard.create for details
	
	Args:
		templateObjCreator: (CalcMethodFactoryBase)
		templateStdInpCreator: (CreatorWithResetableKwargsTemplate) where .create leads to a StandardInputObj object

	Returns
		outFunct: (CallableGetStandardInpCreatorFromOverideDictsStandard object) This is a callable instance
	
	"""
	return CallableGetStandardInpCreatorFromOverideDictsStandard(templateObjCreator, templateStdInpCreator)


class CallableGetStandardInpCreatorFromOverideDictsStandard():
	""" Callable class for getting a standardInput creator from two override dictionaries. See self._create for the interface

	"""


	def __init__(self, templateObjCreator, templateStdInpCreator):
		""" Initializer
		
		Args:
			templateObjCreator: (CalcMethodFactoryBase)
			templateStdInpCreator: (CreatorWithResetableKwargsTemplate) where .create leads to a StandardInputObj object
				
		"""
		self.templateObjCreator = templateObjCreator
		self.templateStdInpCreator = templateStdInpCreator

	def _create(self, objCreatorOverideDict, stdInpCreatorOverideDict):
		""" Create the stdInput creator. This assumes the creator object will have a .baseCreator attribute for the calcObjCreator object to be set to
		
		Args:
			objCreatorOverideDict: (dict) Each key,val pair represents the attribute/value to set on self.templateObjCreator
			stdInpCreatorOverideDict: (dict) Each key,val pair represents the attribute/value to set on self.templateStdInpCreator

		Returns
			stdInpCreator: Has a .create method which returns a StandardInputObj
		
		"""
		baseCalcObjCreator = getCreatorWithOveridenDict(self.templateObjCreator, objCreatorOverideDict)
		outStdInpCreator = copy.deepcopy(self.templateStdInpCreator)
		outStdInpCreator.baseCreator = baseCalcObjCreator
		outStdInpCreator = getCreatorWithOveridenDict(outStdInpCreator, stdInpCreatorOverideDict)
		return outStdInpCreator

	def __call__(self, objCreatorOverideDict, stdInpCreatorOverideDict):
		return self._create(objCreatorOverideDict, stdInpCreatorOverideDict)


def getCreatorWithOveridenDict(inpCreator, inpDict):
    outCreator = copy.deepcopy(inpCreator)
    for key in inpDict.keys():
        if key in inpCreator.registeredKwargs:
            setattr(outCreator, key, inpDict[key])
        else:
            raise ValueError("{} is an invalid key".format(key))
    return outCreator



