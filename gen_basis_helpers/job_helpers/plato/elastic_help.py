
""" Module is meant to be VERY specific to plato. """

import os

from .. import elastic_constants as elasHelp
from ...plato import plato_creator as platoCreator
from ...shared import creator_resetable_kwargs as baseCreator


class PlatoHcpElasticStandardInputCreator(baseCreator.CreatorWithResetableKwargsTemplate):

	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs)
	registeredKwargs.add("eleKey")
	registeredKwargs.add("methodKey")
	registeredKwargs.add("structKey")
	registeredKwargs.add("convDatabase")
	registeredKwargs.add("eleDatabase")
	registeredKwargs.add("baseWorkFolder")
	registeredKwargs.add("strainValues")
	registeredKwargs.add("eType")

	def _createFromSelf(self):
		kwargDict = dict()
		kwargDict["baseWorkFolder"], kwargDict["extToWorkFolder"] = self._getBaseAndExtPaths()
		kwargDict["baseGeom"] = getattr(self.eleDatabase, self.eleKey.capitalize()).getPlaneWaveGeom(self.structKey)
		kwargDict["strainValues"] = self.strainValues
		kwargDict["eleKey"], kwargDict["methodKey"], kwargDict["structKey"] = [getattr(self,x) for x in ["eleKey","methodKey","structKey"]]
		kwargDict["creator"] = self._getCreatorWithoutGeom()
		kwargDict["eType"] = self.eType
		factory = elasHelp.HcpElasticStandardInputCreator(**kwargDict)
		return factory.create()

	def _getCreatorWithoutGeom(self):
		convDatabase = getattr(self.convDatabase,self.eleKey.capitalize())
		

		#1) Get the integration grid spacing
		if self.methodKey.startswith("dft_"):
			gridVals = convDatabase.integGridVals.getPrimCellDftGrid(self.structKey)
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

	def _getBaseAndExtPaths(self):
		extPath = os.path.join("elastic",self.eleKey,"hcp",self.methodKey)
		return self.baseWorkFolder, extPath


	#Setting some properties purely for documentation reasons
	@property
	def convDatabase(self):
		""" RefConvergenceDatabase object """
		return self._convDatabase

	@convDatabase.setter
	def convDatabase(self,val):
		self._convDatabase = val

	@property
	def eleDatabase(self):
		""" RefElementalDataBase object """
		return self._eleDatabase


	@eleDatabase.setter
	def eleDatabase(self,val):
		self._eleDatabase = val	
	

