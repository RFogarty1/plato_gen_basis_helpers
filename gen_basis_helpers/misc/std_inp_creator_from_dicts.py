
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

def createFunctToGetAbsGridConvStdInpCreatorFromOverideDicts(absConvVals, fixedRel, calcObjCreatorTemplate, stdInpCreatorTemplate):
	""" Gets a function for creating standard input creator objects from overide dictionaries
	
	Args:
		absConvVals: (iter of ints) Values of absolute grid convergence parameters
		fixedRel: (int) The value to set relative grid convergence to
		calcObjCreatorTemplate: (CalcMethodFactoryBase)
		stdInpCreatorTemplate: (CreatorWithResetableKwargsTemplate) Specifically should really be ParsedFileObjsForMultiGeomsStandardInputCreator

	Returns
		outFunct: (CallableGetStandardInpCreatorFromOverideDictsStandard object) This is a callable instance with f(objCreatorOverideDict, stdInpCreatorOverideDict)
	
	"""
	args = [calcObjCreatorTemplate, stdInpCreatorTemplate, "relGridCutoff", "absGridCutoff", fixedRel, absConvVals]
	return CallableStandardInpCreatorFromOverideDictsForCP2KGridConv(*args)


def createFunctToGetRelGridConvStdInpCreatorFromOverideDicts(relConvVals, fixedAbs, calcObjCreatorTemplate, stdInpCreatorTemplate):
	""" Gets a function for creating standard input creator objects frmo overide dictionaries
	
	Args:
		relConvVals: (iter of ints) Values of relative grid convergence parameters
		fixedAbs: (int) The value to set absolute grid convergence to
		calcObjCreatorTemplate: (CalcMethodFactoryBase)
		stdInpCreatorTemplate: (CreatorWithResetableKwargsTemplate) Specifically should really be ParsedFileObjsForMultiGeomsStandardInputCreator

	Returns
		outFunct: (CallableGetStandardInpCreatorFromOverideDictsStandard object) This is a callable instance with f(objCreatorOverideDict, stdInpCreatorOverideDict)
	
	"""
	args = [calcObjCreatorTemplate, stdInpCreatorTemplate, "absGridCutoff", "relGridCutoff", fixedAbs, relConvVals]
	return CallableStandardInpCreatorFromOverideDictsForCP2KGridConv(*args)

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


class CallableStandardInpCreatorFromOverideDictsForCP2KGridConv():

	def __init__(self, templateObjCreator, templateStdInpCreator, fixedAttr, varyAttr, fixedVal, varyVals):
		""" Initializer
		
		Args:
			templateObjCreator: (CalcMethodFactoryBase)
			templateStdInpCreator: (CreatorWithResetableKwargsTemplate) Specifically should really be ParsedFileObjsForMultiGeomsStandardInputCreator
			fixedAttr: (str) Either absGridCutoff or relGridCutoff (Depending which you want to keep fixed)
			varyAttr: (str) Either absGridCutoff or relGridCutoff (Depending which varies)
			fixedVal: (int) Cutoff value for fixedAttr
			varyVals: (iter of ints) Cutoff values for the varying attribute

		"""
		self.templateObjCreator = templateObjCreator
		self.templateStdInpCreator = templateStdInpCreator
		self.fixedAttr = fixedAttr
		self.varyAttr = varyAttr
		self.fixedVal = fixedVal
		self.varyVals = list(varyVals)

	def _create(self, objCreatorOverideDict, stdInpCreatorOverideDict):
		""" Create the stdInput creator. This assumes the creator object will have a .baseCreator attribute for the calcObjCreator object to be set to
		
		Args:
			objCreatorOverideDict: (dict) Each key,val pair represents the attribute/value to set on self.templateObjCreator
			stdInpCreatorOverideDict: (dict) Each key,val pair represents the attribute/value to set on self.templateStdInpCreator

		Returns
			stdInpCreator: Has a .create method which returns a StandardInputObj
		
		"""
		outObjs = list()
		for varyVal in self.varyVals:
			objCreator = getCreatorWithOveridenDict(self.templateObjCreator, objCreatorOverideDict)		
			setattr(objCreator, self.fixedAttr, self.fixedVal)
			setattr(objCreator, self.varyAttr, varyVal)
			stdInpCreator = getCreatorWithOveridenDict(self.templateStdInpCreator, stdInpCreatorOverideDict)
			absConv, relConv = self._getAbsConvAndRelConv(self.fixedVal, varyVal)
			stdInpCreator.baseFileName = "grid_convs_a{}_r{}".format(absConv,relConv)
			stdInpCreator.baseCreator = objCreator
			outObjs.append(stdInpCreator)
		return outObjs

	def _getAbsConvAndRelConv(self, fixedVal, varyVal):
		if self.fixedAttr=="absGridCutoff" and self.varyAttr=="relGridCutoff":
			absConv, relConv = fixedVal, varyVal
		elif self.fixedAttr=="relGridCutoff" and self.varyAttr=="absGridCutoff":
			absConv, relConv = varyVal, fixedVal
		else:
			raise ValueError("self.fixedAttr={} and self.varyAttr={} are invalid values".format(self.fixedAttr, self.varyAttr))
		return absConv, relConv

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



