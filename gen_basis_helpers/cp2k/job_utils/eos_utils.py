
import os

from .. import basis_register as basReg
from .. import method_register as methReg
from .. import cp2k_creator as cp2kCreator
from ...shared import label_objs as labelHelp
from ...shared import creator_resetable_kwargs as baseCreator
from ...shared import calc_runners as calcRunners
from ...workflows import base_flow as baseFlow
from ...workflows import eos_workflow as eosFlow

class CP2KEosStandardObjCreator(baseCreator.CreatorWithResetableKwargsTemplate):
	""" Purpose of this is to create the StandardInput object for CP2K EoS calculations using minimal input """
	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs)
	registeredKwargs.add("workFolder")
	registeredKwargs.add("absGridCutoff")
	registeredKwargs.add("relGridCutoff")
	registeredKwargs.add("cp2kMethodStr")
	registeredKwargs.add("basisStrDict")
	registeredKwargs.add("basisAlias")
	registeredKwargs.add("structStrs")
	registeredKwargs.add("addedMOs")
	registeredKwargs.add("maxFunctEvals")
	registeredKwargs.add("eosStr")
	registeredKwargs.add("structStrParamMapper")
	registeredKwargs.add("eleStr")

	#This is essentially main(); the parent class handles temporarily setting kwargs to user values
	def _createFromSelf(self):
		calcObjCreator = self._getCreatorObject()
		fitFunct = self._getFitFunction()
		workflow = self._createWorkflowFromCreatorAndFitFunct(calcObjCreator,fitFunct)
		label = self._createLabel()
		return calcRunners.StandardInputObj(workflow, label)
	
	def _setDefaultInitAttrs(self):
		self.eosStr = "murnaghan"
		self.maxFunctEvals = 10000

	def _createLabel(self):
		return labelHelp.StandardLabel(eleKey=self.eleStr, methodKey=self.basisAlias, structKey="eos")


	def _createWorkflowFromCreatorAndFitFunct(self, creator, fitFunct):
		factory =  CP2KEosWorkflowCreator(calcObjCreator=creator, eosFitFunction=fitFunct,
		                              structStrParamMapper=self.structStrParamMapper,
		                              eleKey=self.eleStr, methodKey=self.basisAlias, structStrs=self.structStrs)
		return factory.create()

	def _getCreatorObject(self):
		#Some of our kwargs correspond directly to options (with different names) in the creator objectt
		ourKwargsToCreatorKwargMaps = {"workFolder":"folderPath",
		                               "cp2kMethodStr":"methodStr"}
		outKwargDict = dict()

		#First get all our keywords (and vals) which are ALSO in the factory kwargs
		sharedKwargs = self.registeredKwargs.intersection(cp2kCreator.CP2KCalcObjFactoryStandard.registeredKwargs)
		for sharedKwarg in sharedKwargs:
			outKwargDict[sharedKwarg] = getattr(self,sharedKwarg)

		#Then we need to add in any options that arent in our kwargs but need setting
		for kwarg in ourKwargsToCreatorKwargMaps.keys():
			val = getattr(self, kwarg)
			outKey = ourKwargsToCreatorKwargMaps[kwarg]
			outKwargDict[outKey] = val

		outKwargDict["basisObjs"] = self._getBasisObjects()

		return cp2kCreator.CP2KCalcObjFactoryStandard(**outKwargDict)

	def _getBasisObjects(self):
		outObjs = list()
		for eleKey,basKey in self.basisStrDict.items():
			outObjs.append( basReg.createCP2KBasisObjFromEleAndBasisStr(eleKey, basKey) )
		return outObjs

	def _getFitFunction(self):
		return eosFlow.StandardEosFitFunction(eosStr=self.eosStr, maxFev=self.maxFunctEvals)


class CP2KEosWorkflowCreator(baseCreator.CreatorWithResetableKwargsTemplate):
	""" Purpose of this is to create a workflow object for running EoS calculations with CP2K """
	registeredKwargs = set(baseCreator.CreatorWithResetableKwargsTemplate.registeredKwargs)
	registeredKwargs.add("calcObjCreator")
	registeredKwargs.add("structStrParamMapper")
	registeredKwargs.add("eosFitFunction")
	registeredKwargs.add("eleKey") #Needs to be set but value likely rarely relevant
	registeredKwargs.add("methodKey") #Needs to be set but value likely rarely relevant
	registeredKwargs.add("structStrs")

	def _createFromSelf(self):
		outObjs = list()
		for structStr in self.structStrs:
			outObjs.append( self._createSingleWorkflow(structStr) )
		return baseFlow.StandardLabelledWorkflowComposite( outObjs )

	def _createSingleWorkflow(self, structStr):
		calcObjs = self._createCalcObjsSingleWorkflow(structStr)
		label = labelHelp.StandardLabel(eleKey=self.eleKey, methodKey=self.methodKey,
		                                structKey=structStr)
		return eosFlow.EosWorkflow(calcObjs, self.eosFitFunction, label)

	def _createCalcObjsSingleWorkflow(self, structStr):
		geoms = self.structStrParamMapper.getGeomsForStructStr( structStr )
		startFolder = self.calcObjCreator.folderPath
		outFolder = os.path.join(startFolder, structStr)
		extraKwargs = self.structStrParamMapper.getKwargDictForStructStr( structStr )
		outObjs = list()
		for x in geoms:
			outObjs.append( self.calcObjCreator.create(geom=x, folderPath=outFolder, fileName=self._getFileNameFromGeom(x),
			                                           **extraKwargs) )
		return outObjs

	def _getFileNameFromGeom(self,geom):
		vol = geom.volume
		volStr = "vol_{:.3f}".format(vol).replace(".","pt")
		return volStr

	#Properties which are implemented for the sake of having some docstrings
	@property
	def calcObjCreator(self):
		return self._calcObjCreator

	@calcObjCreator.setter
	def calcObjCreator(self, val):
		""" CP2KCalcObjFactoryStandard object. Geom and fileName will be set as part of create. Anything else needs setting either in the input object or from parameters returned by structStrParamMapper.
		"""
		self._calcObjCreator = val

	@property
	def eosFitFunction(self):
		""" Function (usually callable class) representing an eos fit. The interface of the call function is:
	
	Args:
		volumes: (float iter) Volume per atom values for structures used to calculate energy-volume curves. Units = bohr^3
		energies: (float iter) Energy per atom values for structures used to calculated energy-volume curves. Units = eV
	Returns
		fitDict: A dictionary containing various parameters relating to the EoS fit. GPa used for bulk-modulus, bohr^3 per atom for volume

		"""
		return self._eosFitFunction

	@eosFitFunction.setter
	def eosFitFunction(self,val):
		self._eosFitFunction = val



class StructStrToParamMapper():
	"""Class is responsible for taking structure strs (e.g. hcp) and returning {key,val} for any specific properties we need setting for ALL relevant calcs + the set of geoms for standard EoS calcs

	"""

	def __init__(self, structDatabase, convDatabase):
		""" Initializer
		
		Args:
			structDatabase: (RefElementalDataBase obj) Contains a getStructsForEos(self,key) function
			convDatabase: (RefConvergenceDatabase) Contains information on k-points and integration grids to use for various structures; we only care about k-points for CP2K (CP2K grid convergences arent in the integ grid class anyway)
	
		"""
		self.structDatabase = structDatabase
		self.convDatabase = convDatabase


	def getGeomsForStructStr(self, structStr):
		""" Gets a set of geometries to use for calculating equation of state for structStr
		
		Args:
			structStr: (str) Alias for the set of EoS geometries we want; e.g. hcp
				
		Returns
			outGeoms: (iter of plato_pylib UnitCell objects) List of geometries to return	
		"""
		return self.structDatabase.getStructsForEos(structStr)

	def getKwargDictForStructStr(self, structStr):
		""" Gets a dictionary, which can be passed as kwargs to CP2K creator object, containing any options we need to modify for structStr (e.g. the number of k-points is the most likely thing)
		
		Args:
			structStr: (str) Alias for the set of EoS geometries we want; e.g. hcp
				
		Returns
			outDict: (dict) This can be passed to CP2k creator objects create function as kwargs to overwrite any defaults with values specific for these calculations
		
		"""
		kPts = self.convDatabase.kptGridVals.getKptsPrimCell(structStr)
		outDict = {"kPts":kPts}
		return outDict
