
import os

from ..shared import calc_runners as calcRunners
from ..shared import creator_resetable_kwargs as baseCreator
from ..shared import label_objs as labelHelp
from ..workflows import self_point_defects as defectFlow

class CodeSpecificStandardInputCreatorTemplate(baseCreator.CreatorWithResetableKwargsTemplate):

	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs)

	registeredKwargs.add("baseWorkFolder")
	registeredKwargs.add("eleKey")
	registeredKwargs.add("methodKey")
	registeredKwargs.add("structKey")
	registeredKwargs.add("bulkGeom")
	registeredKwargs.add("kPtsBulk")
	registeredKwargs.add("defectGeom")
	registeredKwargs.add("kPtsDefect")

	def _createFromSelf(self):
		bulkObj = self._getBulkCreator().create()
		defectObj = self._getDefectCreator().create()
		labelObj = labelHelp.StandardLabel(eleKey=self.eleKey, methodKey=self.methodKey,
		                                   structKey=self.structKey)
		workflow = defectFlow.SelfPointDefectWorkflow(defectObj,bulkObj)
		return calcRunners.StandardInputObj(workflow, labelObj)

	# This should be overwritten by code-specific versions
	def _createCalcObjCreatorBulk(self):
		""" Return a CalcMethodFactoryBase instance with no kPts or geom present
		"""
		raise NotImplementedError("")

	def _createCalcObjCreatorDefect(self):
		""" Return a CalcMethodFactoryBase instance with no kPts or geom present
		"""
		raise NotImplementedError("")


	def _getDefectCreator(self):
		outCreator = self._createCalcObjCreatorDefect()
		outCreator.workFolder = self.outFolder
		outCreator.fileName = "defect_calc"
		outCreator.kPts = self.kPtsDefect
		outCreator.geom = self.defectGeom
		return outCreator

	def _getBulkCreator(self):
		outCreator = self._createCalcObjCreatorBulk()
		outCreator.workFolder = self.outFolder
		outCreator.fileName = "bulk_calc"
		outCreator.kPts = self.kPtsBulk
		outCreator.geom = self.bulkGeom
		return outCreator

	@property
	def outFolder(self):
		return os.path.join(self.baseWorkFolder, self.eleKey, self.structKey, self.methodKey)
