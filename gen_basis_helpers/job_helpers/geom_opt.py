

import copy
from . import standard_template_obj as stdTemplate

from ..shared import calc_runners as calcRunners
from ..workflows import geom_opt_workflow as goptFlow

class CodeSpecificStandardInputCreatorTemplate(stdTemplate.StandardInputCreatorTemplateBase):

	registeredKwargs = set(stdTemplate.StandardInputCreatorTemplateBase.registeredKwargs)
	registeredKwargs.add("baseCreator")
	registeredKwargs.add("fileName")

	def _createFromSelf(self):
		workflow = self._getWorkflow()
		label = self.label
		outObj = calcRunners.StandardInputObj(workflow,label)
		return outObj

	def _getWorkflow(self):
		creator = self._getCreatorObj()
		calcObj = creator.create()
		outWorkflow = goptFlow.GeomOptWorkflow(calcObj)
		return outWorkflow
	
	def _getCreatorObj(self):
		outCreator = self._getBaseCreator()
		outCreator.runType = "geomOpt"
		outCreator.workFolder=  self.outFolder
		outCreator.fileName = self.fileName if self.fileName is not None else "geom_opt"
		return outCreator

	def _getBaseCreator(self):
		""" Creates a CalcMethodFactoryBase with settings that are shared (e.g. cutoff energy) between all calculations. This part generally needs to be overwritten by subclasses. Note that many of these paramters (e.g. k-points) will be overwritten during the create() method	
		"""
		if self.baseCreator is not None:
			return copy.deepcopy(self.baseCreator)
		else:
			raise ValueError("baseCreator attribute not set")

	@property
	def baseCreator(self):
		""" (CalcMethodFactoryBase) Object. Most of the settings in this object will be used for ALL calculations, though things like folder names will be overridden
		"""
		return self._baseCreator

	@baseCreator.setter
	def baseCreator(self,val):
		self._baseCreator = val

