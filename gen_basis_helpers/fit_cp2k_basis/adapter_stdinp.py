
import types

from ..shared import calc_runners as calc_runners

from . import core

class GenericMapFunctionForGettingObjFunctFromStdInp(calc_runners.StandardMapFunction):

	def __init__(self, rawOutputToObjFunct):
		""" Initializer
		
		Args:
			rawOutputToObjFunct: function with call(workflow.output)->objFunctVal interface
				 
		"""
		self.rawOutputToObjFunct = rawOutputToObjFunct 

	def _getObjFunctVal(self, workflowOutput):
		return self.rawOutputToObjFunct(workflowOutput)

	def __call__(self, stdInpObj):
		stdInpObj.workflow.run()
		objFunctVal = self._getObjFunctVal(stdInpObj.workflow.output)
		outObj = types.SimpleNamespace(objFunct=objFunctVal, rawData=stdInpObj.workflow.output)
		return outObj


class WorkflowOutputToObjFunctValStandard():

	def _getValuesFromOutput(self, output):
		raise NotImplementedError("")

	def _getObjFunctionValueFromOutputValues(self, outputVals):
		return self._getObjFunctValFromTargValsAndActVals(self.targVals, outputVals)

	def _getObjFunctValFromTargValsAndActVals(self, targVals, actVals):
		raise NotImplementedError("")

	def __call__(self, output):
		outputVals = self._getValuesFromOutput(output)
		return self._getObjFunctionValueFromOutputValues(outputVals)

