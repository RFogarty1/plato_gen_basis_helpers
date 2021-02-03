
""" Purpose of this module is to make it simpler to create CP2KCalcObj objects """

import contextlib
import os

from . import method_register as methRegister
from . import basis_register as basRegister
from . import cp2k_file_helpers as fileHelpers
from . import cp2k_calc_objs as calcObjs
from . import cp2k_basis_obj as basisObjHelp
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
	registeredKwargs.add("workFolder")
	registeredKwargs.add("folderPath") #DEPRECATED KWARG. workFolder should be used instead to be concsistent with CalcMethodFactoryBase
	registeredKwargs.add("fileName")
	registeredKwargs.add("absGridCutoff")
	registeredKwargs.add("relGridCutoff")
	registeredKwargs.add("printAOMullikenPop")
	registeredKwargs.add("charge")
	registeredKwargs.add("runType")
	registeredKwargs.add("geomConstraints")
	registeredKwargs.add("extraModPyC2pkOpts")
	registeredKwargs.add("saveRestartFile")
	registeredKwargs.add("epsScf")
	registeredKwargs.add("fragmentsBSSE")
	registeredKwargs.add("xcFunctional")
	registeredKwargs.add("grimmeDisp")
	registeredKwargs.add("nonLocalDisp")
	registeredKwargs.add("surfaceDipoleCorr")
	registeredKwargs.add("mdOpts")
	registeredKwargs.add("walltime")
	registeredKwargs.add("extrapolationMethod")
	registeredKwargs.add("print_every_n_md_steps")
	registeredKwargs.add("print_every_n_scf_steps")
	registeredKwargs.add("restart_file_every_n_md_steps")
	registeredKwargs.add("prefDiagLib")
	registeredKwargs.add("epsDef")

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

		#Do a check for workFolder/folderPath. If one is None but the other is not then we need to NOT depend on arbitrarily argument ordering
		outKwargs = {k:v for k,v in kwargs.items()}
		outKwargs["workFolder"] = outKwargs.get("workFolder", None)
		outKwargs["folderPath"] = outKwargs.get("folderPath", None)
		if (outKwargs["workFolder"] is None) and (outKwargs["folderPath"] is not None):
			self.workFolder = outKwargs["folderPath"]
		elif (outKwargs["workFolder"] is not None) and (outKwargs["folderPath"] is None):
			self.workFolder = outKwargs["workFolder"]
		else:
			pass #If both are None, then obviously do nothing


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


	#These need properties so that the old "folderPath" kwarg and the new "workFolder" kwargs are linked
	@property
	def folderPath(self):
		""" Now just an alias for workFolder; which is the attr whcih SHOULD actually be set """
		return self.workFolder #The new kwarg

	@folderPath.setter
	def folderPath(self,val):
		self.workFolder = val

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
		if self.runType is None:
			md = False
		else:
			md = True if self.runType.lower()=="md" else False
		keepRestartFile = True if self.saveRestartFile is None else self.saveRestartFile #Usually not ever written, so passing False can cause issues (attempt to rm a non-existent file can throw an error)
		outputObj = calcObjs.CP2KCalcObj(basicObj, basePath=self._getPathToPassCalcObj(), saveRestartFile=keepRestartFile, md=md)
		return outputObj

	def _modPycp2kObj(self,pycp2kObj):
		#Modify basis set info and geometry; these need a special function essentially
		basisObjs = self._getBasisObjs()
		fileHelpers.addGeomAndBasisInfoToSimpleCP2KObj(pycp2kObj, self.geom, basisObjs)

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
		if self.printAOMullikenPop:
			modDict["printAOMullikenPop".lower()] = True
		if self.charge is not None:
			modDict["charge"] = self.charge
		if self.epsScf is not None:
			modDict["epsScf"] = self.epsScf
		if self.xcFunctional is not None:
			modDict["xcFunctional".lower()] = self.xcFunctional
		if self.grimmeDisp is not None:
			currDict = self.grimmeDisp.modPyCP2KDict
			modDict.update(currDict)
		if self.nonLocalDisp is not None:
			currDict = self.nonLocalDisp.modPyCP2KDict
			modDict.update(currDict)
		if self.surfaceDipoleCorr is not None:
			currDict = self.surfaceDipoleCorr.modPyCP2KDict
			modDict.update(currDict)
		if self.walltime is not None:
			modDict["walltime"] = self.walltime
		if self.extrapolationMethod is not None:
			modDict["qsExtrapolationMethod"] = self.extrapolationMethod
		if self.print_every_n_md_steps is not None:
			modDict["trajPrintEachMd"] = self.print_every_n_md_steps
		if self.print_every_n_scf_steps is not None:
			modDict["trajPrintEachScf"] =  self.print_every_n_scf_steps
		if self.restart_file_every_n_md_steps is not None:
			modDict["restartPrintEachMd"] = self.restart_file_every_n_md_steps
		if self.prefDiagLib is not None:
			modDict["prefDiagLib"] = self.prefDiagLib
		if self.epsDef is not None:
			modDict["epsDef"] = self.epsDef

		modDict["scfPrintRestart".lower()] = False

		runTypeModDict = self._getModDictBasedOnRunType()
		modDict.update(runTypeModDict)

		if self.extraModPyC2pkOpts is not None:
			modDict.update(self.extraModPyC2pkOpts)

		fileHelpers.modCp2kObjBasedOnDict(pycp2kObj, modDict)


	def _getBasisObjs(self):
		basisObjs = self.basisObjs
		if self.runType is not None:
			if self.runType.lower() == "bsse":
				basisObjs = basisObjHelp.getBasisObjsWithGhostVersionsIncluded(self.basisObjs)
		return basisObjs

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


	def _getModDictBasedOnRunType(self):
		outDict = dict()
		runStr = self.runType if self.runType is not None else "None" #Need to be able to apply .lower() to it

		if runStr.lower() == "geomOpt".lower():
			if self.geomConstraints is None:
				outDict["runType".lower()] = "cell_opt" #Later, may need to pass geo_opt (or similar) if constraints set a certain way
			else:
				geomOptDict = getCP2KModDictBasedOnGeomConstraints(self.geomConstraints)
				outDict.update(geomOptDict)

		if runStr.lower() == "bsse":
			outDict["runType".lower()] = "bsse"
			outDict["fragmentsBSSE".lower()] = self.fragmentsBSSE

		if runStr.lower() == "md":
			outDict = self.mdOpts.optDict
			outDict["runType".lower()] = "md"

		return outDict


def getCP2KModDictBasedOnGeomConstraints(geomConstraints):
	outDict = dict()
	if geomConstraints.constraintsPresent is False:
		outDict["runtype"] = "cell_opt"
		return outDict

	outDict = _getModDictBasedOnCellConstraints(geomConstraints.cellConstraints)

	if outDict == dict():
		raise ValueError("Cant handle geomConstraints object")
	return outDict

def _getModDictBasedOnCellConstraints(cellConstraints):
	outDict = dict()
	if all(cellConstraints.anglesToFix) and all(cellConstraints.lattParamsToFix):
		outDict["runtype"] = "geo_opt"
		return outDict

	if all(cellConstraints.anglesToFix):
		outDict["geo_constrain_cell_angles"] = [True,True,True] #Constrain all angles case
		outDict["runtype"] = "cell_opt"
	return outDict

