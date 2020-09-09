
from . import core as coreFit
from . import adapter_stdinp as stdInpAdapterHelp

import types

class ParsedFileWorkflowToCondNumberObjFunction(stdInpAdapterHelp.WorkflowOutputToObjFunctValStandard):

	def __init__(self, targVals, targAndActValsToObjVal, applyFunctToCondNumbs=None):
		""" Initializer
		
		Args:
			targVals: (iter of float) Target values for the overlap workflow; will usually be length 1
			targAndActValsToObjVal: f(targVals,actVals)->objFunctVal
			applyFunctToCondNumbs: f(condNumb)->y If this is set then the function is applied to the condition number. Original purpose is to allow log10(cond-number) to be used
		"""
		self.targVals = targVals
		self.targAndActValsToObjVal = targAndActValsToObjVal
		self.applyFunctToCondNumbs = applyFunctToCondNumbs if applyFunctToCondNumbs is not None else lambda x:x

	def _getValuesFromOutput(self, output):
		condVals = [x.parsedFile.overlap_condition_number.diag.twoNorm for x in output]
		outVals = [self.applyFunctToCondNumbs(x) for x in condVals]
		return outVals

	def _getObjFunctValFromTargValsAndActVals(self, targVals, actVals):
		return self.targAndActValsToObjVal(targVals,actVals)


class ParsedFileWorkflowToAtomEnergyObjFunction(stdInpAdapterHelp.WorkflowOutputToObjFunctValStandard):

	def __init__(self, targVals, targAndActValsToObjVal):
		""" Initializer
		
		Args:
			targVals: (iter of float) Target values for the overlap workflow; will usually be length 1
			targAndActValsToObjVal: f(targVals,actVals)->objFunctVal
		"""
		self.targVals = targVals
		self.targAndActValsToObjVal = targAndActValsToObjVal

	def _getValuesFromOutput(self, output):
		electronicE = [x.parsedFile.energies.electronicTotalE for x in output]
		return electronicE

	def _getObjFunctValFromTargValsAndActVals(self, targVals, actVals):
		return self.targAndActValsToObjVal(targVals,actVals)





