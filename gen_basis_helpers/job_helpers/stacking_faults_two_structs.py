

import itertools as it
import types

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
		kwargDict = self._getKwargDictForModdingStackFaultCreator()
		return self.baseCreatorStackFaultStruct.create(**kwargDict)

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



class MapStackingFaultTwoStructsToUsefulFormatStandard():

	def __init__(self, stackFaultEnergyFmt = "{:.4f}", xVal="methodStr", multStackFaultByFactor=1):
		self.xVal = xVal
		self.stackFaultEnergyFmt = stackFaultEnergyFmt
		self.multStackFaultByFactor = multStackFaultByFactor

	def _getTableData(self, stdInputObj):
		xVal = self._getXValFromStdInp(stdInputObj)
		yVal = self.stackFaultEnergyFmt.format( stdInputObj.workflow.output[0].stackFaultEnergy*self.multStackFaultByFactor )
		return [xVal,yVal]

	def _getXValFromStdInp(self, stdInputObj):
		if self.xVal == "methodStr":
			outVal = stdInputObj.label[0].methodKey
		else:
			raise ValueError("{} is an invalid value for xVal".format(self.xVal))
		return outVal

	def __call__(self, stdInputObj):
		stdInputObj.workflow.run()
		assert len(stdInputObj.workflow.output)==1
		assert len(stdInputObj.label)==1
		output = types.SimpleNamespace(tableData=None)
		output.tableData = self._getTableData(stdInputObj)
		return output		

