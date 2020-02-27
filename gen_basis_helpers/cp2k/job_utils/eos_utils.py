
from ...shared import label_objs as labelHelp
from ...shared import creator_resetable_kwargs as baseCreator
from ...workflows import base_flow as baseFlow
from ...workflows import eos_workflow as eosFlow

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
		extraKwargs = self.structStrParamMapper.getKwargDictForStructStr( structStr )
		outObjs = list()
		for x in geoms:
			outObjs.append( self.calcObjCreator.create(geom=x, fileName=self._getFileNameFromGeom(x),
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
		outDict = {"kpts":kPts}
		return outDict
