
from . import standard_template_obj as stdTemplate
from ..shared import calc_runners as calcRunners
from ..workflows import base_flow as baseFlow
from ..workflows import eos_workflow as eosFlowHelp
from ..shared import label_objs as labelHelp

class StandardInputCreatorEosCurves(stdTemplate.StandardInputCreatorTemplateBase):

	registeredKwargs = set(stdTemplate.StandardInputCreatorTemplateBase.registeredKwargs)

	registeredKwargs.add("structStrs") #Iter of strings
	registeredKwargs.add("baseCreators") #Single creator or iter of strings same length and structStrs
	registeredKwargs.add("geoms") #Iter of iter of geoms (plato_pylib UnitCell); length of top level needs to be the same as structStrs
	registeredKwargs.add("eosFitFunctions") #Single StandardEosFitFunction OR iter of them (to use different functs for each struct)

	def _createFromSelf(self):
		outWorkflow = self._createWorkflow()
		outLabel = labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)
		creatorObj = calcRunners.StandardInputObj(outWorkflow, outLabel)
		return creatorObj.create()

	def _createWorkflow(self):
		allWorkflows = list()
		for structStr in self.structStrs:
			currWorkflow = self._getWorkflowForStructStr(structStr)
			allWorkflows.append(currWorkflow)
		outWorkflow = baseFlow.StandardLabelledWorkflowComposite(allWorkflows)
		return outWorkflow

	def _getWorkflowForStructStr(self, structStr):
		structIdx = self.structStrs.index(structStr)
		calcObjs = self._getCalcObjsForStructStr(structStr)
		outLabel = labelHelp.StandardLabel(eleKey=self.eleKey, structKey=structStr, methodKey=self.methodKey)
		fitObj = self._getIterForAttr("eosFitFunctions")[structIdx]
		outFlow = eosFlowHelp.EosWorkflow(calcObjs, fitObj, outLabel)
		return outFlow

	def _getCalcObjsForStructStr(self, structStr):
		structIdx = self.structStrs.index(structStr)
		allGeoms = self.geoms[structIdx]
		creatorObj = self._getIterForAttr("baseCreators")[structIdx]
		outObjs = list()
		for x in allGeoms:
			currObj = creatorObj.create(geom=x)
			outObjs.append(currObj)
		return outObjs

	def _getIterForAttr(self, attr):
		attrVal = getattr(self,attr)
		try:
			outIter = [x for x in iter(attrVal)]
			assert len(outIter)==len(self.structStrs)
		except TypeError:
			outIter = [attrVal for x in self.structStrs]
		return outIter


