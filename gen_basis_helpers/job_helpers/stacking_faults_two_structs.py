

import itertools as it

from . import standard_template_obj as stdTemplate
from ..shared import calc_runners as calcRunners
from ..workflows import stacking_fault_workflows as stackFaultFlow

class StandardInputCreatorStackFaultTwoStructs(stdTemplate.StandardInputCreatorTemplateBase):

	registeredKwargs = set(stdTemplate.StandardInputCreatorTemplateBase.registeredKwargs)

	registeredKwargs.add("baseCreatorPerfectStruct")
	registeredKwargs.add("baseCreatorStackFaultStruct")
	registeredKwargs.add("perfectStructGeom")
	registeredKwargs.add("stackFaultGeom")

	def _createFromSelf(self):
		workFlow = self._createWorkflow()
		label = self.label
		outObj = calcRunners.StandardInputObj(workFlow, label)
		return outObj


	def _createWorkflow(self):
		perfCalcStruct = self._getCalcObjPerfectStruct()
		stackFaultStruct = self._getCalcObjStackFaultStruct()
		outWorkflow = stackFaultFlow.StackingFaultWorkflowTwoStructs(perfCalcStruct,stackFaultStruct)
		return outWorkflow

	def _getCalcObjStackFaultStruct(self):
		kwargDict = self.__getKwargDictForModdingStackFaultCreator()
		return self.baseCreatorStackFaultStruct(**kwargDict)

	def _getCalcObjPerfectStruct(self):
		kwargDict = self._getKwargDictForModdingPerfectStructCreator()
		return self.baseCreatorPerfectStruct.create(**kwargDict)

	def _getKwargDictForModdingPerfectStructCreator(self):
		outDict = dict()
		outDict["geom"] = self.perfectStructGeom
		outDict["fileName"] = "no_fault_struct"
		outDict["workFolder"] = self.outFolder
		return outDict

	def _getKwargDictForModdingStackFaultCreator(self):
		outDict = dict()
		outDict["geom"] = self.stackFaultGeom
		outDict["fileName"] = "stack_fault_geom"
		outDict["workFolder"] = self.outFolder
		return outDict


