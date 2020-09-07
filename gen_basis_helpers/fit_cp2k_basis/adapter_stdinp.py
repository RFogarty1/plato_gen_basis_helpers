
import types

from ..shared import calc_runners as calcRunners

from . import core

#TODO: Add ability to put weights into the map function


def getAdaptedStdInptFromWorkflow(workflow, outputToObjFunct, label=None):
	""" Gets an AdaptedStandardInput object when given a workflow and a function to map its output to an objective function value
	
	Args:
		workflow: (BaseWorkflow object) Used to encapsulate instructions for the calculating something
		outputToObjFunct: (WorkflowOutputToObjFunctValStandard) This maps the output from workflow into an objective function value
			 
	Returns
		stdInp: (AdaptedStandardInput object) A version of () object which creates an output object in the format required for fitting (the overall objective function knows where to find the contributions of this calculation)
 
	"""
	mapFunction = GenericMapFunctionForGettingObjFunctFromStdInp(outputToObjFunct)
	stdInpObj = calcRunners.StandardInputObj( workflow, label, mapFunction=mapFunction )
	return stdInpObj

class GenericMapFunctionForGettingObjFunctFromStdInp(calcRunners.StandardMapFunction):

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

