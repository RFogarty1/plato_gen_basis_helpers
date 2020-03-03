
import os

from ..shared import creator_resetable_kwargs as baseCreator
from ..shared import calc_runners as calcRunners
from ..shared import label_objs as labelHelp
from ..workflows import elastic_workflows as elasticFlow

class HcpElasticStandardInputCreator(baseCreator.CreatorWithResetableKwargsTemplate):
	""" Class for creating standard input objects for calculating Hcp elastic constants

	"""

	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs)
	registeredKwargs.add("baseGeom")
	registeredKwargs.add("creator")
	registeredKwargs.add("eleKey")
	registeredKwargs.add("structKey")
	registeredKwargs.add("methodKey")
	registeredKwargs.add("strainValues")
	registeredKwargs.add("baseWorkFolder")
	registeredKwargs.add("extToWorkFolder")

	def _createFromSelf(self):
		workFlow = self._createWorkFlow()
		label = self._createLabel()
		return calcRunners.StandardInputObj( workFlow, label )

	def _createWorkFlow(self):
		factory = elasticFlow.HcpElasticWorkflowCreator(baseGeom=self.baseGeom, creator=self.creator,
		                                                        strainValues=self.strainValues,
		                                                        workFolder=self._outFolder)
		return factory.create()

	def _createLabel(self):
		return labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)

	def _setDefaultInitAttrs(self):
		self.extToWorkFolder = "elastic_hcp"
		self.structKey = "hcp"

	@property
	def _outFolder(self):
		extension = self.extToWorkFolder if self.extToWorkFolder is not None else ""
		return os.path.join(self.baseWorkFolder, extension)

