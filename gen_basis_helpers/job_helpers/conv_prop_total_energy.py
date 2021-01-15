
#Code to help create calculations for converging total energy with respect to some property (for example, plane-wave cutoff or the k-points used)
import copy
import itertools as it
import types

from ..shared import calc_runners as calcRunners
from ..shared import label_objs as labelHelp
from ..shared import creator_resetable_kwargs as baseCreator
from ..workflows import convergers as convFlow

import os

class CodeSpecificStandardInputCreatorTemplate(baseCreator.CreatorWithResetableKwargsTemplate):

	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs)
	registeredKwargs.add("baseWorkFolder")
	registeredKwargs.add("eleKey")
	registeredKwargs.add("methodKey")
	registeredKwargs.add("structKey")
	registeredKwargs.add("geom")
	registeredKwargs.add("mapConvValToNumber") #Only needed if the convergence values arent single numbers
	registeredKwargs.add("baseCreatorObj")
	registeredKwargs.add("convVals")
#	registeredKwargs.add("modCreatorFunct")
	registeredKwargs.add("convKwarg")

	def _createFromSelf(self):
		self._checkConvValsUnique()
		calcObjs = self._getCalcObjs()
		mappedConvVals = [self._getNumberFromConvVal(x) for x in self.convVals]
		outWorkflow = convFlow.GridConvergenceEnergyWorkflow (calcObjs, mappedConvVals)
		labelObj = labelHelp.StandardLabel(eleKey=self.eleKey, methodKey=self.methodKey,
		                                   structKey=self.structKey)
		return calcRunners.StandardInputObj(outWorkflow, labelObj)


	@property
	def outFolder(self):
		return os.path.join(self.baseWorkFolder, self.eleKey, self.structKey, self.methodKey)

	@property
	def _outBaseFileNames(self):
		convFmt = "conv_val_{:.3f}"
		convNumbs = [self._getNumberFromConvVal(x) for x in self.convVals]
		return [convFmt.format(x).replace(".","pt") for x in convNumbs]

	def _checkConvValsUnique(self):
		outFileNames = self._outBaseFileNames
		numbNames, numbUnique = len(outFileNames), len(set(outFileNames))
		if numbNames != numbUnique:
			raise ValueError("ConvVals dont lead to unique fileNames, we get {} names for {} files. Conv vals are {}".format(numbUnique,numbNames, self.convVals))

	def _getNumberFromConvVal(self, convVal):
		if self.mapConvValToNumber is None:
			return convVal
		else:
			return self.mapConvValToNumber(convVal)

	def _getCalcObjs(self):
		creators = self._getCreatorsForWorkflow()
		return [x.create() for x in creators]

	def _getCreatorsForWorkflow(self):
		outCreators = self._getBasicCreatorObjForEachRequiredCalc()
		
		#Modify with shared options
		for creator in outCreators:
			self._modifyCreatorObjWithSharedOptions(creator)

		#Modify with new filenames (which depend on convVals)
		for creator, fileName in it.zip_longest(outCreators, self._outBaseFileNames):
			creator.fileName = fileName

		#Modify with the actual convergence values
		for creator, convVal in it.zip_longest(outCreators, self.convVals):
			setattr(creator,self.convKwarg,convVal) #Can easily add a modCreatorFunct later if more complexity needed ever

		return outCreators

	def _modifyCreatorObjWithSharedOptions(self, creatorObj):
		creatorObj.geom = self.geom
		creatorObj.workFolder = self.outFolder

	def _getBasicCreatorObjForEachRequiredCalc(self):
		return [copy.deepcopy(self.baseCreatorObj) for x in self.convVals]



	#Properties created for the sake of documentation
	@property
	def baseCreatorObj(self):
		""" (CalcMethodFactoryBase) Object which defines all the parameters for each calculation. At minimum, the geometry, fileName and workFolder attributes wont be used
		"""
		return self._baseCreatorObj

	@baseCreatorObj.setter
	def baseCreatorObj(self,val):
		self._baseCreatorObj = val

	@property
	def mapConvValToNumber(self):
		""" Function f(x) used to convert each convVal to a number (float). Needed in generating fileNames and creating the workflow. No need to set if each convVal is ALREADY a number (e.g. when converging plane-wave cutoffs). An example of when this is needed is when dealing with k-points	
		"""
		return self._mapConvValToNumber

	@mapConvValToNumber.setter
	def mapConvValToNumber(self, val):
		self._mapConvValToNumber = val





class MapConvergersWorkflowToUsefulFormatStandard():
	""" Callable class used to map the output from a total energy convergence workflow, assuming its part of a StandardInput object, to a more convenient format """

	def __init__(self, deltaE=True):
		""" Initializer
		
		Args:
			deltaE: (bool) If True then convert the output energies into relative values. The reference is taken as the last value by default( presuming that convVals is passed in order of least to most converged)
				 
		"""
		self.deltaE=deltaE

	def _getStdOutData(self, stdInpObj):
		xVals = [x[0] for x in stdInpObj.workflow.output.convResults]		
		yVals = [x[1] for x in stdInpObj.workflow.output.convResults]
		deltaYVals = [y-yVals[-1] for y in yVals] #Assuming that we want it relative to the LAST (most converged) valu
		if self.deltaE:
			return [[x,y] for x,y in it.zip_longest(xVals, deltaYVals)]
		else:
			return [[x,y] for x,y in it.zip_longest(xVals, yVals)]

	def __call__(self, stdInpObj):
		stdInpObj.workflow.run()
		outObj = types.SimpleNamespace(outData=None)
		outObj.outData = self._getStdOutData(stdInpObj)
		return outObj









