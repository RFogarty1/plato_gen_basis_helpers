
import os
from ..shared import creator_resetable_kwargs as baseCreator
from ..shared import label_objs as labelHelp
from ..shared import calc_runners as calcRunners
from ..workflows import total_energies as totEnergyFlow
from ..workflows import reaction_energies as reactFlow

class CodeSpecificStandardInputCreatorTemplate(baseCreator.CreatorWithResetableKwargsTemplate):

	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs)

	registeredKwargs.add("baseWorkFolder")
	registeredKwargs.add("eleKey")
	registeredKwargs.add("methodKey")
	registeredKwargs.add("structKey")
	registeredKwargs.add("bulkWithAdsorbedGeom")
	registeredKwargs.add("bulkWithoutAdsorbedGeom")
	registeredKwargs.add("gasPhaseReactantGeoms")
	registeredKwargs.add("gasPhaseProductGeoms")
	registeredKwargs.add("gasPhaseReactantStoichiometries") #Note the bulk always have stoichiometries of 1
	registeredKwargs.add("gasPhaseProductStoichiometries")
	registeredKwargs.add("kPtsBulk")

	def _createFromSelf(self):
		outLabel = labelHelp.StandardLabel(eleKey=self.eleKey, methodKey=self.methodKey, structKey=self.structKey)
		outWorkflow = self._getReactionWorkflow()
		outObj = calcRunners.StandardInputObj( outWorkflow, outLabel )
		return outObj

	#Code-specific; please overwrite.
	def _getBaseCreator(self):
		""" Creates a CalcMethodFactoryBase with settings that are shared (e.g. cutoff energy) between all calculations. This part needs to be overwritten by subclasses. Note that many of these paramters (e.g. k-points) will be overwritten during the create() method	
		"""
		raise NotImplementedError("")

	def _getReactionWorkflow(self):
		reactants, products = self._getReactantsWorkflow(), self._getProductsWorkflow()
		reactionFlow = reactFlow.ReactionEnergyWorkflow(reactants, products)
		return reactionFlow

	def _getReactantsWorkflow(self):
		reactantWorkflows = self._getReactantTotalEnergyWorkflows()
		reactantStoics = [1] + self.gasPhaseReactantStoichiometries
		assert len(reactantWorkflows)==len(reactantStoics)
		return totEnergyFlow.TotalEnergyGroupWorkflow( reactantWorkflows, reactantStoics )

	def _getProductsWorkflow(self):
		productWorkflows = self._getProductTotalEnergyWorkflows()
		productStoics = [1] + self.gasPhaseProductStoichiometries
		assert len(productWorkflows)==len(productStoics)
		return totEnergyFlow.TotalEnergyGroupWorkflow( productWorkflows, productStoics )

	def _getReactantTotalEnergyWorkflows(self):
		reactantCalcObjs = self._getReactantCalcObjs()
		outFlows = [totEnergyFlow.TotalEnergyWorkflow(x) for x in reactantCalcObjs]
		return outFlows

	def _getProductTotalEnergyWorkflows(self):
		productCalcObjs = self._getProductCalcObjs()
		outFlows = [totEnergyFlow.TotalEnergyWorkflow(x) for x in productCalcObjs]
		return outFlows

	def _getReactantCalcObjs(self):
		bulkCalcObj = self._getBulkWithoutAdsorbedGeomCreator().create()
		gasPhaseFactories = self._getGasPhaseReacantObjs()
		gasPhaseCalcObjs = [x.create() for x in gasPhaseFactories]
		return [bulkCalcObj] + gasPhaseCalcObjs

	def _getProductCalcObjs(self):
		bulkCalcObj = self._getBulkWithAdsorbedGeomCreator().create()
		gasPhaseFactories = self._getGasPhaseProductObjs()
		gasPhaseCalcObjs = [x.create() for x in gasPhaseFactories]
		return [bulkCalcObj] + gasPhaseCalcObjs

	def _getGasPhaseReacantObjs(self):
		numbObjs = len(self.gasPhaseReactantGeoms)
		outObjs = self._getNGasPhaseCreatorObjsWithSharedOptions( numbObjs )
		for idx,x in enumerate(outObjs):
			x.fileName = "gas_phase_reactant_{}".format(idx)
			x.geom = self.gasPhaseReactantGeoms[idx]
		return outObjs

	def _getGasPhaseProductObjs(self):
		numbObjs = len(self.gasPhaseProductGeoms)
		outObjs = self._getNGasPhaseCreatorObjsWithSharedOptions( numbObjs )
		for idx,x in enumerate(outObjs):
			x.fileName = "gas_phase_product_{}".format(idx)
			x.geom = self.gasPhaseProductGeoms[idx]
		return outObjs

	def _getNGasPhaseCreatorObjsWithSharedOptions(self, numbObjs):
		outObjs = [self._getBaseCreator() for x in range(numbObjs)]
		for x in outObjs:
			x.workFolder = self.outFolder
			x.kPts = [1,1,1]
		return outObjs

	def _getBulkWithAdsorbedGeomCreator(self):
		outObj = self._getBaseCreator()
		outObj.geom = self.bulkWithAdsorbedGeom
		outObj.fileName = "bulk_with_adsorbed"
		self._modifyBulkCreatorWithSharedOptions(outObj)
		return outObj

	def _getBulkWithoutAdsorbedGeomCreator(self):
		outObj = self._getBaseCreator()
		outObj.geom = self.bulkWithoutAdsorbedGeom
		outObj.fileName = "bulk_calculation"
		self._modifyBulkCreatorWithSharedOptions(outObj)
		return outObj

	def _modifyBulkCreatorWithSharedOptions(self, creator):
		creator.workFolder = self.outFolder
		creator.kPts = self.kPtsBulk


	@property
	def gasPhaseReactantStoichiometries(self):
		""" iter of floats, each element is the stoichiometry of gas phase reactants. Note that the bulk calculations have fixed stoichiometries of 1. Will return len 0 list if not set to something (or if set to None)
		
		"""
		if self._gasPhaseReactantStoichiometries is None:
			return list()
		else:
			return self._gasPhaseReactantStoichiometries

	@gasPhaseReactantStoichiometries.setter
	def gasPhaseReactantStoichiometries(self, val):
		self._gasPhaseReactantStoichiometries = val

	@property
	def gasPhaseProductStoichiometries(self):
		""" iter of floats, each element is the stoichiometry of gas phase products. Note that the bulk calculations have fixed stoichiometries of 1. Will return len 0 list if not set to something (or if set to None)
		
		"""
		if self._gasPhaseProductStoichiometries is None:
			return list()
		else:
			return self._gasPhaseProductStoichiometries

	@gasPhaseProductStoichiometries.setter
	def gasPhaseProductStoichiometries(self,val):
		self._gasPhaseProductStoichiometries = val

	#This is the folder we run all calculations in.
	@property
	def outFolder(self):
		return os.path.join(self.baseWorkFolder, self.eleKey, self.methodKey, self.structKey)

