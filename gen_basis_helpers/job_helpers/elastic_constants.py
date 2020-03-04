
import itertools as it
import os
import types

import numpy as np

from ..shared import creator_resetable_kwargs as baseCreator
from ..shared import calc_runners as calcRunners
from ..shared import label_objs as labelHelp
from ..workflows import elastic_workflows as elasticFlow

class HcpElasticStandardInputCreator(baseCreator.CreatorWithResetableKwargsTemplate):
	""" Class for creating standard input objects for calculating Hcp elastic constants

	"""

	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs)
	registeredKwargs.add("baseGeom")
	registeredKwargs.add("creator")
	registeredKwargs.add("eleKey")
	registeredKwargs.add("structKey")
	registeredKwargs.add("methodKey")
	registeredKwargs.add("strainValues")
	registeredKwargs.add("baseWorkFolder")
	registeredKwargs.add("extToWorkFolder")
	registeredKwargs.add("eType")

	def _createFromSelf(self):
		workFlow = self._createWorkFlow()
		label = self._createLabel()
		return calcRunners.StandardInputObj( workFlow, label )

	def _createWorkFlow(self):
		factory = elasticFlow.HcpElasticWorkflowCreator(baseGeom=self.baseGeom, creator=self.creator,
		                                                        strainValues=self.strainValues,
		                                                        workFolder=self._outFolder, eType=self.eType)
		return factory.create()

	def _createLabel(self):
		return labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)

	def _setDefaultInitAttrs(self):
		self.extToWorkFolder = "elastic_hcp"
		self.structKey = "hcp"

	@property
	def _outFolder(self):
		extension = self.extToWorkFolder if self.extToWorkFolder is not None else ""
		return os.path.join(self.baseWorkFolder, extension)



class MapElasticflowOutputToUsefulFormatStandard():
	"""Callable class used to transform ElasticFlow output into more useful format for a single StandardInputObject

	   The callable interface takes the StandardInput object as the sole argument

	"""

	def __init__(self, structStr, fitStrainVals=None):
		allowedStructStrs = ["hcp"]
		if structStr != "hcp":
			raise ValueError("{} is currently an invalid structure type for this mapping function; valid values are {}".format(structStr, allowedStructStrs))
		self.structStr = structStr

		self.fitStrainVals = fitStrainVals if fitStrainVals is not None else None #None just means figure out the default at runtime


	def _getTableData(self, stdInputObj):
		assert len(stdInputObj.label)==1
		elastics = ["{:.2f}".format(x) for x in stdInputObj.workflow.output[0].elasticConsts.values()]
		total = [stdInputObj.label[0].methodKey] + elastics
		return total

	def _getPlotData(self, stdInputObj):
		fitStrainVals = self.fitStrainVals if self.fitStrainVals is not None else self._getFitStrainValsFromInpObj(stdInputObj)
		#Step 1 = get the fit data
		allFitData = list()
		for x in stdInputObj.workflow.output[0].stressStrainData:
			yVals = x.fitFunct(fitStrainVals)
			fitData = [[x,y] for x,y in it.zip_longest(fitStrainVals,yVals)]
			allFitData.append(fitData)

		#Step 2 = get the actual data
		allActData = list()
		for x in stdInputObj.workflow.output[0].stressStrainData:
			allActData.append( x.actVals )

		#Step 3 = get the combined (actData,fitData) array
		allOutputData = list()
		for act,fit in it.zip_longest(allActData,allFitData):
			allOutputData.append( [act,fit] )

		return allOutputData


	def _getFitStrainValsFromInpObj(self,stdInputObj):
		allStressStrainData = stdInputObj.workflow.output[0].stressStrainData
		allStrainVals = list()
		for x in allStressStrainData:
			allStrainVals.extend( [a[0] for a in x.strainVsEnergy] )

#		allStrainVals = [x.strainVsEnergy[0] for x in allStressStrainData]
		minCoeff = min(allStrainVals)
		maxCoeff = max(allStrainVals)
		step = 0.001
		return np.arange(minCoeff, maxCoeff, step)

	def __call__(self, stdInputObj):
		stdInputObj.workflow.run()
		assert len(stdInputObj.workflow.output)==1
		outputObj = types.SimpleNamespace(plotData=None, tableData=None, allData=None,
		                                  tableHeadings=None, strainStrs=None)
		outputObj.tableHeadings = ["Method", "c11", "c12", "c13", "c33", "c44"]
		outputObj.tableData = self._getTableData(stdInputObj)
		outputObj.strainStrs = [x.toStr() for x in stdInputObj.workflow.output[0].strains]
		outputObj.allData = stdInputObj.workflow.output[0]
		outputObj.plotData = self._getPlotData(stdInputObj)
		return outputObj #Should return an output object
