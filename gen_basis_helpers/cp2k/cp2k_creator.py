
""" Purpose of this module is to make it simpler to create CP2KCalcObj objects """

import os

from . import method_register as methRegister
from . import basis_register as basRegister
from . import cp2k_file_helpers as fileHelpers
from . import cp2k_calc_objs as calcObjs
from ..shared import data_plot_base as dPlotBase


#Should really inherit from CalcMethodFactoryBase, but slightly annoying to refactor it
class BaseCP2KCalcObjFactory():

	def create(self, **kwargs):
		""" Create a CP2KCalcObj object based on attributes of this object
		
		Args:
			**kwargs: These are attributes of the factory object. Passing them here is equivalent to setting them on this object, calling create and then resetting them to the original value. For example:
			outObj = obj.create(geom=someGeom)
			
			is the same as
			
			origGeom = obj.geom
			obj.geom = someGeom
			outObj = obj.create()
			obj.geom = origGeom

		Returns
			CP2KCalcObj

		Raises:
			KeyError: If an unregistered kwarg is passed
			ValueError: If a required option isnt set within the object itself (e.g. if geom isnt set, the implementation MIGHT decided to throw this error - note this is the Base class docstring)
		
		"""
		raise NotImplementedError("")


class CP2KCalcObjFactoryStandard(BaseCP2KCalcObjFactory):

	registeredKwargs = set()
	registeredKwargs.add("methodStr")
	registeredKwargs.add("requiredArgsToBeSet")
	registeredKwargs.add("kPts")
	registeredKwargs.add("addedMOs")
	registeredKwargs.add("geom")
	registeredKwargs.add("basisObjs")
	registeredKwargs.add("folderPath")
	registeredKwargs.add("fileName")
	registeredKwargs.add("absGridCutoff")
	registeredKwargs.add("relGridCutoff")

	def __init__(self,**kwargs):
		""" Initializer for CP2K calc-object factory
		
		Raises:
			KeyError: If an unregistered kwarg is passed
		"""
		#First initialise all arguments
		for key in self.registeredKwargs:
			setattr(self,key,None)

		#Decide which arguments NEED to be set in the final creator; If one isnt set this will throw once create() is called
		self._baseReqArgsToBeSet = ["methodStr","geom"] #These ALWAYS need setting; unless user decides to overwrite this hidden attrib
		self._additionalReqArgsToBeSet = set()

		#Now set all arguments, if the keyword is valid
		self._updateAttrsFromKwargs(**kwargs)


	def _updateAttrsFromKwargs(self, **kwargs):
		for key in kwargs:
			if key in self.registeredKwargs:
				setattr(self,key,kwargs[key])
			else:
				raise KeyError("{} is an invalid keyword.\n Available kwargs are {}".format(key , self.registeredKwargs))


	@property
	def requiredArgsToBeSet(self):
		""" These are attributes which NEED to be set (i.e. not None) when calling the self.create() function.

		Note: When setting this attribute you are actually only setting values ON TOP OF those defined in self._baseReqArgsToBeSet. i.e. Even if you set this property to None there will still be some required attributes for the create function
		
		"""
		return list( set(self._baseReqArgsToBeSet).union(self._additionalReqArgsToBeSet) )

	@requiredArgsToBeSet.setter
	def requiredArgsToBeSet(self, vals):
		if vals is None:
			self._additionalReqArgsToBeSet = set()
		else:
			self._additionalReqArgsToBeSet = set(vals)

	#Making these properties purely so they can have doc-strings
	@property
	def geom(self):
		""" Geometry of the system as represented by a plato_pylib UnitCell object
		"""
		return self._geom

	@geom.setter
	def geom(self,val):
		self._geom = val

	@property
	def basisObjs(self):
		""" An iter of CP2KBasisObjBase objects, each represents the basis set to use for one element """
		return self._basisObj

	@basisObjs.setter
	def basisObjs(self,val):
		self._basisObj = val


	def _ensureReqArgsAllSet(self):
		for key in self.requiredArgsToBeSet:
			if getattr(self,key) is None:
				raise ValueError("{} Needs to be set to create a CP2K Calc Object".format(key))

	def create(self, **kwargs):
		with dPlotBase.temporarilySetDataPlotterRegisteredAttrs(self,kwargs):
			self._ensureReqArgsAllSet()
			outObj = self._createOutputObj()
		return outObj

	def _createOutputObj(self):
		basicObj = methRegister.createCP2KObjFromMethodStr(self.methodStr)
		self._modPycp2kObj(basicObj)
		outputObj = calcObjs.CP2KCalcObj(basicObj, basePath=self._getPathToPassCalcObj())
		return outputObj

	def _modPycp2kObj(self,pycp2kObj):
		#Modify basis set info and geometry; these need a special function essentially
		fileHelpers.addGeomAndBasisInfoToSimpleCP2KObj(pycp2kObj, self.geom, self.basisObjs)

		#Modify any remaining properties we care about
		modDict = dict()
		if self.kPts is not None:
			modDict["kpts"] = self.kPts
		if self.addedMOs is not None:
			modDict["addedMOs".lower()] = self.addedMOs
		if self.absGridCutoff is not None:
			modDict["gridCutAbs".lower()] = self.absGridCutoff
		if self.relGridCutoff is not None:
			modDict["gridCutRel".lower()] = self.relGridCutoff

		fileHelpers.modCp2kObjBasedOnDict(pycp2kObj, modDict)

	def _getPathToPassCalcObj(self):

		if (self.folderPath is not None) and (self.fileName is not None):
			return os.path.join(self.folderPath, self.fileName)

		elif (self.folderPath is not None) or (self.fileName is not None):
			raise AssertionError("Both folderPath AND fileName need to be set; or neither need to be set" 
			                     "(i.e. both should be None) but folderPath={},fileName={}".format(self.folderPath,self.fileName))

		elif (self.folderPath is None) and (self.fileName is None):
			return None

		else:
			raise AssertionError("Not sure why, but something went wrong with figuring out a file path")

