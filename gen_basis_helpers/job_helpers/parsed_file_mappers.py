import types

class StandardInptXYMapperForParsedFileWorkflows():
	
	def __init__(self, xFunct, yFunct):
		self.xFunct = xFunct
		self.yFunct = yFunct

	def _getXFunctVals(self, stdInputObj):
		return self.xFunct(stdInputObj)
	
	def _getYFunctVals(self, stdInputObj):
		return self.yFunct(stdInputObj)
	
	def __call__(self, stdInputObj):
		stdInputObj.workflow.run()
		output = types.SimpleNamespace()
		output.xVals = self._getXFunctVals(stdInputObj)
		output.yVals = self._getYFunctVals(stdInputObj)
		return output

class StandardInptMapperForVolumesVsDiagCondNumber(StandardInptXYMapperForParsedFileWorkflows):
	
	def __init__(self, divVolsByConstant=None, ignoreNoneValues=False):
		self.divVolByConstant = divVolsByConstant if divVolsByConstant is not None else 1
		self.ignoreNoneValues = ignoreNoneValues
	
	def _getXFunctVals(self, stdInputObj):
		return getVolumesFromStdInputObj(stdInputObj, divVolByConstant=self.divVolByConstant, ignoreNoneValues=self.ignoreNoneValues)
	
	def _getYFunctVals(self,stdInputObj):
		return getDiagConditionNumberFromStdInputObj(stdInputObj, ignoreNoneValues=self.ignoreNoneValues)


class StandardInptMapperForVolumeVsEnergy(StandardInptXYMapperForParsedFileWorkflows):

	def __init__(self, divVolsByConstant=None, ignoreNoneValues=False):
		self.divVolByConstant = divVolsByConstant if divVolsByConstant is not None else 1
		self.ignoreNoneValues = ignoreNoneValues

	def _getXFunctVals(self, stdInputObj):
		return getVolumesFromStdInputObj(stdInputObj, divVolByConstant=self.divVolByConstant, ignoreNoneValues=self.ignoreNoneValues)

	def _getYFunctVals(self,stdInputObj):
		return getTotalEnergiesFromStdInputObj(stdInputObj, ignoreNoneValues=self.ignoreNoneValues)


#TODO: Factor out the duplication between these functions
def getVolumesFromStdInputObj(stdInputObj, divVolByConstant=None, ignoreNoneValues=False):
	if ignoreNoneValues:
		allGeoms = list()
		for x in stdInputObj.workflow.output:
			if x.parsedFile is not None:
				allGeoms.append(x.parsedFile.unitCell)	
	else:
		allGeoms = [x.parsedFile.unitCell for x in stdInputObj.workflow.output]

	allVolumes = [x.volume/divVolByConstant for x in allGeoms]
	return allVolumes


def getDiagConditionNumberFromStdInputObj(stdInputObj, ignoreNoneValues=False):
	if ignoreNoneValues:
		allVals = list()
		for x in stdInputObj.workflow.output:
			if x.parsedFile is not None:
				allVals.append(x.parsedFile.overlap_condition_number.diag.twoNorm)
	else:
		allVals = [x.parsedFile.overlap_condition_number.diag.twoNorm for x in stdInputObj.workflow.output]

	return allVals


def getTotalEnergiesFromStdInputObj(stdInputObj, ignoreNoneValues=False):
	if ignoreNoneValues:
		allVals = list()
		for x in stdInputObj.workflow.output:
			if x.parsedFile is not None:
				allVals.append(x.parsedFile.energies.electronicTotalE)
	else:
		allVals = [x.parsedFile.energies.electronicTotalE for x in stdInputObj.workflow.output]

	return allVals










