

import copy
import types
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
	registeredKwargs.add("reactantCharges")
	registeredKwargs.add("productCharges")

	def _createFromSelf(self):
		workflow = self._getReactionEnergyWorkflow()
		label = self.label
		return calcRunners.StandardInputObj(workflow, label)

	def _getBaseCreator(self):
		""" Creates a CalcMethodFactoryBase with settings that are shared (e.g. cutoff energy) between all calculations. This part generally needs to be overwritten by subclasses. Note that many of these paramters (e.g. k-points) will be overwritten during the create() method	
		"""
		if self.baseCreator is not None:
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
			creator.charge = self.reactantCharges[idx]
		outObjs = [x.create() for x in startCreators]
		return outObjs

	def _getProductCalcObjs(self):
		startCreators = self._getCreatorsForSetOfGeoms(self.productGeoms)
		for idx,creator in enumerate(startCreators):
			creator.fileName = "product_{}".format(idx)
			creator.charge = self.productCharges[idx]
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



class MapGasPhaseReactionEnergyWorkflowToUsefulFormatStandard():
	"""Callable class used to transform ReactionEnergyWorkflow output into a better format for tabulating

	   The callable interface takes a non-composite StandardInput object as the sole argument
	"""

	def __init__(self, reactEnergyFmt="{:.3f}", reactantEnergyFmt="{:.3f}",
	             productEnergyFmt="{:.3f}"):
		""" Initializer
		
		Args:
			reactEnergyFmt: (str,optional) Format str used for the total reaction energy format
				 
		"""
		self.reactEnergyFmt = reactEnergyFmt
		self.reactantEnergyFmt = reactantEnergyFmt
		self.productEnergyFmt = productEnergyFmt


	def _getSimpleTableHeadings(self):
		return ["Method", "Reaction Energy (eV)"]

	def _getSimpleTableData(self, stdInputObj):
		methodLabel = stdInputObj.label[0].methodKey
		reactEnergy = stdInputObj.workflow.output[0].energy
		return [methodLabel, self.reactEnergyFmt.format(reactEnergy)]

	def _getTableWithEnergyBreakdownsData(self, stdInputObj):
		outData = self._getSimpleTableData(stdInputObj)
		reactantEnergies = stdInputObj.workflow.output[0].reactantEnergies
		productEnergies = stdInputObj.workflow.output[0].productEnergies
		outData += [self.reactantEnergyFmt.format(x) for x in reactantEnergies]
		outData += [self.productEnergyFmt.format(x) for x in productEnergies] 
		return outData

	def _getTableWithEnergyBreakdownsHeadings(self, stdInputObj):
		outHeadings = self._getSimpleTableHeadings()
		numbReactants = len( stdInputObj.workflow.output[0].reactantEnergies )
		numbProducts = len( stdInputObj.workflow.output[0].productEnergies )
		outHeadings += ["Reactant {} Energy (eV)".format(x+1) for x in range(numbReactants)]
		outHeadings += ["Product {} Energy (eV)".format(x+1) for x in range(numbProducts)]
		return outHeadings

	def __call__(self, stdInputObj):
		stdInputObj.workflow.run()
		assert len(stdInputObj.workflow.output)==1
		assert len(stdInputObj.label)==1
		output = types.SimpleNamespace(tableData=None, tableHeadings=None,
		                               tableWithBreakdownData=None, tableWithBreakdownHeadings=None)
		output.tableData = self._getSimpleTableData(stdInputObj)
		output.tableHeadings = self._getSimpleTableHeadings()
		output.tableWithBreakdownData = self._getTableWithEnergyBreakdownsData(stdInputObj)
		output.tableWithBreakdownHeadings = self._getTableWithEnergyBreakdownsHeadings(stdInputObj)

		return output

