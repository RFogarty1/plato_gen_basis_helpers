

import copy
from . import standard_template_obj as stdTemplate

from ..shared import calc_runners as calcRunners
from ..workflows import reaction_energies as reactFlow
from ..workflows import total_energies as totEnergyFlow

class CodeSpecificStandardInputCreatorTemplate(stdTemplate.StandardInputCreatorTemplateBase):

	registeredKwargs = set(stdTemplate.StandardInputCreatorTemplateBase.registeredKwargs)

	registeredKwargs.add("reactantGeoms")
	registeredKwargs.add("productGeoms")
	registeredKwargs.add("reactantStoichiometries")
	registeredKwargs.add("productStoichiometries")
	registeredKwargs.add("baseCreator")

	def _createFromSelf(self):
		workflow = self._getReactionEnergyWorkflow()
		label = self.label
		return calcRunners.StandardInputObj(workflow, label)

	def _getBaseCreator(self):
		""" Creates a CalcMethodFactoryBase with settings that are shared (e.g. cutoff energy) between all calculations. This part generally needs to be overwritten by subclasses. Note that many of these paramters (e.g. k-points) will be overwritten during the create() method	
		"""
		if baseCreator is not None:
			return copy.deepcopy(self.baseCreator)
		else:
			raise ValueError("baseCreator attribute not set")

	def _getReactionEnergyWorkflow(self):
		reactantFlow, productFlow = self._getReactantsTotalEnergyWorkflow(), self._getProductsTotalEnergyWorkflow()
		return reactFlow.ReactionEnergyWorkflow(reactantFlow, productFlow)


	def _getReactantsTotalEnergyWorkflow(self):
		calcObjs = self._getReactantCalcObjs()
		workflows = [totEnergyFlow.TotalEnergyWorkflow(x) for x in calcObjs]
		assert len(workflows)==len(self.reactantStoichiometries)
		return totEnergyFlow.TotalEnergyGroupWorkflow(workflows, self.reactantStoichiometries)

	def _getProductsTotalEnergyWorkflow(self):
		calcObjs = self._getProductCalcObjs()
		workflows = [totEnergyFlow.TotalEnergyWorkflow(x) for x in calcObjs]
		assert len(workflows)==len(self.productStoichiometries)
		return totEnergyFlow.TotalEnergyGroupWorkflow(workflows, self.productStoichiometries)

	def _getReactantCalcObjs(self):
		startCreators = self._getCreatorsForSetOfGeoms(self.reactantGeoms)
		for idx,creator in enumerate(startCreators):
			creator.fileName = "reactant_{}".format(idx)
		outObjs = [x.create() for x in startCreators]
		return outObjs

	def _getProductCalcObjs(self):
		startCreators = self._getCreatorsForSetOfGeoms(self.productGeoms)
		for idx,creator in enumerate(startCreators):
			creator.fileName = "product_{}".format(idx)
		outObjs = [x.create() for x in startCreators]
		return outObjs

	def _getCreatorsForSetOfGeoms(self, geoms):
		outCreators = list()
		for geom in geoms:
			currCreator = self._getBaseCreator()
			currCreator.geom = geom
			self._applySharedOptionsToCreator(currCreator)
			outCreators.append(currCreator)
		return outCreators

	def _applySharedOptionsToCreator(self, creator):
		creator.workFolder = self.outFolder
		creator.kPts = [1,1,1]


	#Properties for the sake of doc-stringing
	@property
	def reactantStoichiometries(self):
		""" iter of floats, each element is the stoichiometry of a reactant. """
		return self._reactantStoichiometries

	@reactantStoichiometries.setter
	def reactantStoichiometries(self, val):
		self._reactantStoichiometries = val

	@property
	def productStoichiometries(self):
		""" iter of floats, each element is the stoichiometry of a reactant. """
		return self._productStoichiometries

	@productStoichiometries.setter
	def productStoichiometries(self, val):
		self._productStoichiometries = val

	@property
	def baseCreator(self):
		""" (CalcMethodFactoryBase) Object. Most of the settings in this object will be used for ALL calculations, though things like folder names and kPts will be overridden
		"""
		return self._baseCreator

	@baseCreator.setter
	def baseCreator(self,val):
		self._baseCreator = val

