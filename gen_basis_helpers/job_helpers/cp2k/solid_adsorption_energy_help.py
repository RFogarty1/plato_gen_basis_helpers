
from ...cp2k import cp2k_creator as cp2kCreatorModule
from .. import solid_adsorption_energy as solidAdsorb



class SolidAdsorptionEnergyStandardInputCreator(solidAdsorb.CodeSpecificStandardInputCreatorTemplate):

	registeredKwargs = set(solidAdsorb.CodeSpecificStandardInputCreatorTemplate.registeredKwargs)
	registeredKwargs.add("absGridCutoff")
	registeredKwargs.add("addedMOsBulk")
	registeredKwargs.add("basisObjs")
	registeredKwargs.add("cp2kMethodStr")
	registeredKwargs.add("relGridCutoff")
	registeredKwargs.add("printAOMullikenPop")

	#Takes care of MOST the options (since im forcing cutoffs to be the same for all)
	def _getBaseCreator(self):
		kwargDict = {"absGridCutoff":self.absGridCutoff, "basisObjs":self.basisObjs,
		             "methodStr":self.cp2kMethodStr, "relGridCutoff":self.relGridCutoff,
		             "printAOMullikenPop":self.printAOMullikenPop}
		return cp2kCreatorModule.CP2KCalcObjFactoryStandard(**kwargDict)


	#Needed for bulk specific options
	def _modifyBulkCreatorWithSharedOptions(self,creator):
		super()._modifyBulkCreatorWithSharedOptions(creator)
		if self.addedMOsBulk is not None:
			creator.addedMOs = self.addedMOsBulk



