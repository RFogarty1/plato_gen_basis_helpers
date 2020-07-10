
import itertools as it

from . import standard_template_obj as stdTemplate
from ..shared import calc_runners as calcRunners
from ..workflows import stacking_fault_workflows as stackFaultFlow

class StandardInputCreatorTemplate(stdTemplate.StandardInputCreatorTemplateBase):

	registeredKwargs = set(stdTemplate.StandardInputCreatorTemplateBase.registeredKwargs)

	registeredKwargs.add("baseCreator")
	registeredKwargs.add("perfectCellGeom")
	registeredKwargs.add("dispZeroGeom") #Often the same as perfectCellGeom
	registeredKwargs.add("dispVals")
	registeredKwargs.add("fitterObj")
	registeredKwargs.add("stackingFaultGeomGenerator") #BaseStackingFaultGeomGenerator object

	def _createFromSelf(self):
		outWorkflow = self._createWorkflow()
		outLabel = self.label
		return calcRunners.StandardInputObj(outWorkflow, outLabel)

	def _createWorkflow(self):
		perfectCellCalcObj = self._createPerfectCalcObj()
		dispObjs = self._createAllDispCalcObjs()
		return stackFaultFlow.StackingFaultWorkflow(perfectCellCalcObj, dispObjs, self.dispVals, self.fitterObj)

	def _createPerfectCalcObj(self):
		kwargDict = self._getKwargDictForModdingPerfectCellCreatorObj()
		return self.baseCreator.create(**kwargDict)

	def _createAllDispGeoms(self):
		outGeoms = list()
		dispFunct = self.stackingFaultGeomGenerator.getGeomForGivenDisplacement
		for disp in self.dispVals:
			currGeom = dispFunct(self.dispZeroGeom, disp)
			outGeoms.append(currGeom)
		return outGeoms

	def _createAllDispCalcObjs(self):
		outObjs = list()
		allGeoms = self._createAllDispGeoms()
		for geom,dVal in it.zip_longest(allGeoms,self.dispVals):
			currKwargDict = self._getKwargDictForModdingDispCreatorObj(geom,dVal)
			currObj = self.baseCreator.create(**currKwargDict)
			outObjs.append(currObj)
		return outObjs

	def _getKwargDictForModdingPerfectCellCreatorObj(self):
		outDict = dict()
		outDict["workFolder"] = self.outFolder
		outDict["fileName"] = "perfect_cell"
		outDict["geom"] = self.perfectCellGeom
		return outDict

	def _getKwargDictForModdingDispCreatorObj(self, geom, dispVal):
		outDict = dict()
		dispStr = str(dispVal).replace(".","pt")
		outDict["workFolder"] = self.outFolder
		outDict["fileName"] = "disp_val_{}".format(dispStr)
		outDict["geom"] = geom
		return outDict
