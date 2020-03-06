
""" Module is meant to be VERY specific to plato. """

from .. import elastic_constants as elasHelp
from ...plato import plato_creator as platoCreator

class PlatoHcpElasticStandardInputCreator(elasHelp.CodeSpecificStandardInputCreatorTemplate):

	registeredKwargs = set(elasHelp.CodeSpecificStandardInputCreatorTemplate.registeredKwargs)
	registeredKwargs.add("convDatabase")


	def _getCreatorWithoutGeom(self):
		convDatabase = getattr(self.convDatabase,self.eleKey.capitalize())
		
		#1) Get the integration grid spacing
		if self.methodKey.startswith("dft_"):
			gridVals = convDatabase.integGridVals.getPrimCellDftFFTGrid(self.structKey)
		#tb1 doesnt need a grid, but it doenst hurt it
		elif self.methodKey.startswith("dft2_") or self.methodKey.startswith("tb1_"):
			gridVals = convDatabase.integGridVals.getPrimCellDft2AngularGrid(self.structKey)
		else:
			raise ValueError("{} is not a recognized methodKey".format(self.methodKey))

		#2) Get the k-points to use
		kPts = convDatabase.kptGridVals.getKptsPrimCell(self.structKey)

		#3) Get the dataset to use
		dataSet = getattr(self.eleDatabase, self.eleKey.capitalize()).modelFiles.dft2PlatoPath

		factory = platoCreator.PlatoCalcObjFactoryStandard(gridVals=gridVals, kPts=kPts,
		                                                   methodStr=self.methodKey, dataSet=dataSet)

		return factory

	#Setting some properties purely for documentation reasons
	@property
	def convDatabase(self):
		""" RefConvergenceDatabase object """
		return self._convDatabase

	@convDatabase.setter
	def convDatabase(self,val):
		self._convDatabase = val

	

