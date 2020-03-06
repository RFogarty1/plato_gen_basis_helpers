
""" Very specific functions for castep """

from .. import elastic_constants as elasHelp
from ...castep import castep_creator as castepCreator

class CastepHcpElasticStandardInputCreator(elasHelp.CodeSpecificStandardInputCreatorTemplate):

	registeredKwargs = set(elasHelp.CodeSpecificStandardInputCreatorTemplate.registeredKwargs)
	registeredKwargs.add("convDatabase")
	registeredKwargs.add("pseudoPotDict")
	registeredKwargs.add("symmetryGenerate")

	def _getCreatorWithoutGeom(self):
		thisEleConvDatabase = getattr(self.convDatabase,self.eleKey.capitalize())

		#1) Get the grid cutoff
		cutoffEnergy = thisEleConvDatabase.getCutoffEnergyPrimCell(self.structKey)

		#2) Get the k-points to use
		kPts = thisEleConvDatabase.getKptsPrimCell(self.structKey)

		#3) Create the factory
		factory = castepCreator.CastepCalcObjFactoryStandard(methodStr=self.methodKey, kPts=kPts, cutoffEnergy=cutoffEnergy,
		                                                     pseudoPotDict=self.pseudoPotDict, symmetryGenerate=self.symmetryGenerate)

		return factory


	def _setDefaultInitAttrs(self):
		super()._setDefaultInitAttrs()
		self.symmetryGenerate=True


class BaseCastepSingleElementConvDatabase():


	def getKptsPrimCell(self, structKey):
		raise NotImplementedError("")


	def getCutoffEnergyPrimCell(self, structKey):
		raise NotImplementedError("")

