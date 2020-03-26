
import os
from ..shared import calc_runners as calcRunners
from ..shared import creator_resetable_kwargs as baseCreator
from ..shared import label_objs as labelHelp

#Main purpose is to factor out a couple of standard methods, such as those used to generate fodler paths and labels
class StandardInputCreatorTemplateBase(baseCreator.CreatorWithResetableKwargsTemplate):

	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs)
	registeredKwargs.add("baseWorkFolder")
	registeredKwargs.add("eleKey")
	registeredKwargs.add("methodKey")
	registeredKwargs.add("structKey")

	@property
	def label(self):
		return labelHelp.StandardLabel(eleKey=self.eleKey, methodKey=self.methodKey, structKey=self.structKey)

	@property
	def outFolder(self):
		return os.path.join(self.baseWorkFolder, self.eleKey, self.structKey, self.methodKey)

