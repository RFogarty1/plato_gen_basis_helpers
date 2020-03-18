

""" CP2K specific code for creating surface energies """

import os
from .. import surface_energies as surfEHelp
from ...cp2k import cp2k_creator as cp2kCreatorHelp


class SurfaceEnergyStandardInputCreator(surfEHelp.CodeSpecificStandardInputCreatorTemplate):

	#TODO: probably split into args which can and CANNOT be passed directly to the creator (or similar)
	registeredKwargs = set(surfEHelp.CodeSpecificStandardInputCreatorTemplate.registeredKwargs)
	registeredKwargs.add("absGridCutoff")
	registeredKwargs.add("addedMOs")
	registeredKwargs.add("basisObjs") #Iter of basis objects
	registeredKwargs.add("cp2kMethodStr")
	registeredKwargs.add("relGridCutoff")
	registeredKwargs.add("extraCreatorKwargDict")

	#This is the important function; i.e. the one that gets called by higher level code
	def _createCalcObjCreator(self):
		kwargDict = { "addedMOs":self.addedMOs, "basisObjs":self.basisObjs, "methodStr":self.cp2kMethodStr,
		              "absGridCutoff":self.absGridCutoff, "relGridCutoff":self.relGridCutoff,
		              "workFolder":self._folderPath }
		if self.extraCreatorKwargDict is not None:
			kwargDict.update( self.extraCreatorKwargDict )
		outCreator = cp2kCreatorHelp.CP2KCalcObjFactoryStandard(**kwargDict)
		return outCreator


	@property
	def _folderPath(self):
		return os.path.join(self.baseWorkFolder, self.eleKey, self.methodKey, self.structKey)

