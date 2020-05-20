
import copy
import os
import types

from . import standard_template_obj as stdTemplate

from ..shared import calc_runners as calcRunners
from ..workflows import self_point_defects as defectFlow

class CodeSpecificStandardInputCreatorTemplate(stdTemplate.StandardInputCreatorTemplateBase):

	registeredKwargs = set(stdTemplate.StandardInputCreatorTemplateBase.registeredKwargs)

	registeredKwargs.add("bulkGeom")
	registeredKwargs.add("kPtsBulk")
	registeredKwargs.add("defectGeom")
	registeredKwargs.add("kPtsDefect")
	registeredKwargs.add("baseCreatorBulk")
	registeredKwargs.add("baseCreatorDefect")


	def _createFromSelf(self):
		bulkObj = self._getBulkCreator().create()
		defectObj = self._getDefectCreator().create()
		workflow = defectFlow.SelfPointDefectWorkflow(defectObj,bulkObj)
		return calcRunners.StandardInputObj(workflow, self.label)

	# This should be overwritten by code-specific versions
	def _createCalcObjCreatorBulk(self):
		""" Return a CalcMethodFactoryBase instance with no kPts or geom present
		"""
		if self.baseCreatorBulk is not None:
			return copy.deepcopy(self.baseCreatorBulk)
		else:
			raise ValueError("baseCreatorBulk not set")

	def _createCalcObjCreatorDefect(self):
		""" Return a CalcMethodFactoryBase instance with no kPts or geom present
		"""
		if self.baseCreatorDefect is not None:
			return copy.deepcopy(self.baseCreatorDefect)
		else:
			raise ValueError("baseCreatorDefect not set")


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


class MapSelfDefectWorkflowOutputToUsefulFormatStandard():
	"""Callable class used to transform SelfPointDefectWorkflow output into a better format for tabulating

		The callable interface takes a non-composite StandardInput object as the sole argument
	"""

	def __init__(self, defectEFormat="{:.2f}", ePerAtomFmtStr="{:.2f}"):
		self.defectEFormat = defectEFormat
		self.ePerAtomFmtStr = ePerAtomFmtStr

	def _getTableData(self, stdInputObj):
		methKey = stdInputObj.label[0].methodKey
		defectE = self.defectEFormat.format(stdInputObj.workflow.output[0].defectE)
		return [methKey, defectE]

	def _getTableDataWithEPerAtom(self,stdInputObj):
		outTable = self._getTableData(stdInputObj)
		ePerAtomBulk = stdInputObj.workflow.output[0].bulkEPerAtom
		ePerAtomDefect = stdInputObj.workflow.output[0].defectEPerAtom
		outTable += [self.ePerAtomFmtStr.format(x) for x in [ePerAtomBulk,ePerAtomDefect]]  
		return outTable	

	def _getTableHeadings(self):
		return ["Basis Set", "Defect Energy (eV)"]

	def _getTableWithEPerAtomHeadings(self):
		return ["Basis Set", "Defect Energy (eV)", "E per atom (bulk, eV)", "E per atom (defect, eV)",]

	def __call__(self, stdInputObj):
		stdInputObj.workflow.run()
		assert len(stdInputObj.workflow.output)==1
		assert len(stdInputObj.label)==1
		output = types.SimpleNamespace(tableData=None, tableWithEPerAtomVals=None,
		                               tableHeaders=None, tableHeadersWithEPerAtom=None)

		output.tableData = self._getTableData(stdInputObj)
		output.tableHeaders = self._getTableHeadings()
		output.tableWithEPerAtomVals = self._getTableDataWithEPerAtom(stdInputObj)
		output.tableHeadersWithEPerAtom = self._getTableWithEPerAtomHeadings()
		return output



