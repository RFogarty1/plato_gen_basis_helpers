

import copy
from . import standard_template_obj as stdTemplate

from ..cp2k import cp2k_creator_to_dict as creatorToDictHelp
from ..shared import calc_runners as calcRunners
from ..workflows import geom_opt_workflow as goptFlow

class CodeSpecificStandardInputCreatorTemplate(stdTemplate.StandardInputCreatorTemplateBase):

	registeredKwargs = set(stdTemplate.StandardInputCreatorTemplateBase.registeredKwargs)
	registeredKwargs.add("baseCreator")
	registeredKwargs.add("fileName")

	def _createFromSelf(self):
		workflow = self._getWorkflow()
		label = self.label
		outObj = calcRunners.StandardInputObj(workflow,label)
		return outObj

	def _getWorkflow(self):
		creator = self._getCreatorObj()
		calcObj = creator.create()
		outWorkflow = goptFlow.GeomOptWorkflow(calcObj)
		return outWorkflow
	
	def _getCreatorObj(self):
		outCreator = self._getBaseCreator()
		outCreator.runType = "geomOpt"
		outCreator.workFolder=  self.outFolder
		outCreator.fileName = self.fileName if self.fileName is not None else "geom_opt"
		return outCreator

	def _getBaseCreator(self):
		""" Creates a CalcMethodFactoryBase with settings that are shared (e.g. cutoff energy) between all calculations. This part generally needs to be overwritten by subclasses. Note that many of these paramters (e.g. k-points) will be overwritten during the create() method	
		"""
		if self.baseCreator is not None:
			return copy.deepcopy(self.baseCreator)
		else:
			raise ValueError("baseCreator attribute not set")

	@property
	def baseCreator(self):
		""" (CalcMethodFactoryBase) Object. Most of the settings in this object will be used for ALL calculations, though things like folder names will be overridden
		"""
		return self._baseCreator

	@baseCreator.setter
	def baseCreator(self,val):
		self._baseCreator = val




def getBasicOutDictFromGeomStdInpCreator(stdInpCreator, calcObjCreatorToDict=None):
	""" Get a dictionary for writing to a database from a stdInpCreator. This assumes jobs have actually been run, and all options were present in the creator factory (rather than some being passed at creation time)
	
	Args:
		stdInpCreator: (job_helpers/geomOpt creator object)
		creatorToDict: (Optional, function) Function which takes a CP2K object creator and returns a dictionary representation

	"""
	creatorToDict = calcObjCreatorToDict if calcObjCreatorToDict is not None else creatorToDictHelp.getSimpleCreatorObjToDictMapObj()
	outDict = creatorToDict(stdInpCreator.baseCreator)
	outDict["compound"] = stdInpCreator.eleKey
	outDict["method"] = stdInpCreator.methodKey
	outDict["structure"] = stdInpCreator.structKey
	currStdOut = stdInpCreator.create().createOutputObj()
	nCreated = len(currStdOut.data[0])
	assert nCreated==1
	outDict["out_geom"] = currStdOut.data[0][0].geom.toDict()
	outDict["energies"] = currStdOut.data[0][0].parsedFile.energies.toDict()
	return outDict









