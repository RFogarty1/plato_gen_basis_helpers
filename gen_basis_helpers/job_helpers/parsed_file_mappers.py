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
	
	def __init__(self, divVolsByConstant=None):
		self.divVolConstant = divVolsByConstant if divVolsByConstant is not None else 1
	
	def _getXFunctVals(self, stdInputObj):
		allGeoms = [x.parsedFile.unitCell for x in stdInputObj.workflow.output]
		allVolumes = [x.volume/self.divVolConstant for x in allGeoms]
		return allVolumes
	
	def _getYFunctVals(self,stdInputObj):
		allVals = [x.parsedFile.overlap_condition_number.diag.twoNorm for x in stdInputObj.workflow.output]
		return allVals

