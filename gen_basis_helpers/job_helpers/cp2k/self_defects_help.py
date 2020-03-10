
from .. import self_defects as selfDefectHelp
from ...cp2k import cp2k_creator as cp2kCreatorHelp


class SelfDefectStandardInputCreator(selfDefectHelp.CodeSpecificStandardInputCreatorTemplate):

	registeredKwargs = set(selfDefectHelp.CodeSpecificStandardInputCreatorTemplate.registeredKwargs)
	registeredKwargs.add("absGridCutoff")
	registeredKwargs.add("addedMOsBulk")
	registeredKwargs.add("addedMOsDefect")
	registeredKwargs.add("basisObjs")
	registeredKwargs.add("cp2kMethodStr")
	registeredKwargs.add("relGridCutoff")

	#Called form above
	def _createCalcObjCreatorBulk(self):
		outCreator = self._createCalcObjCreatorAll()
		outCreator.addedMOs = self.addedMOsBulk
		return outCreator

	#Called from above
	def _createCalcObjCreatorDefect(self):
		outCreator = self._createCalcObjCreatorAll()
		outCreator.addedMOs = self.addedMOsDefect
		return outCreator

	def _createCalcObjCreatorAll(self):
		kwargDict = {"basisObjs":self.basisObjs, "methodStr":self.cp2kMethodStr,
		             "absGridCutoff":self.absGridCutoff, "relGridCutoff":self.relGridCutoff}
		return cp2kCreatorHelp.CP2KCalcObjFactoryStandard(**kwargDict)
