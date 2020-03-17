

from ..shared import creator_resetable_kwargs as baseCreator


class OrbitalBasisSetParserTemplate(baseCreator.CreatorWithResetableKwargsTemplate):

	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs)


	def parse(self, **kwargs):
		""" Parse the required basis set; as determined by object attributes and **kwargs (which overwrite attributes temporarily). This is really an alias for create(); so see that docstring as well
		
		Args:
			kwargs: Must be one of self.registeredKwargs
				 
		Returns
			outBasis: (OrbitalBasisSetBase object) This has all the info needed to plot the radial parts of the basis sets
	 
		"""
		return self.create(**kwargs)


	def _createFromSelf(self):
		return self._parseFromSelf()


	#This needs to be overwritten by the code specific versions
	def _parseFromSelf(self):
		raise NotImplementedError("")
