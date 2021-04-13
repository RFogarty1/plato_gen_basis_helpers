
import copy
import itertools as it
from . import standard_template_obj as stdTemplate
from ..cp2k import cp2k_creator_to_dict as creatorToDictHelp
from ..shared import calc_runners as calcRunners
from ..workflows import parsed_file_workflow as parsedFileFlow
from ..workflows import base_flow as baseFlow

class ParsedFileObjsForMultiGeomsStandardInputCreator(stdTemplate.StandardInputCreatorTemplateBase):

	registeredKwargs = set(stdTemplate.StandardInputCreatorTemplateBase.registeredKwargs)

	registeredKwargs.add("geoms")
	registeredKwargs.add("baseCreator")
	registeredKwargs.add("catchParserErrors")
	registeredKwargs.add("baseFileName")
	registeredKwargs.add("fileNames")

	def _createFromSelf(self):
		workflow = self._getWorkflow()
		outObj = calcRunners.StandardInputObj(workflow, self.label)
		return outObj

	def _getWorkflow(self):
		calcObjs = self._getAllCalcObjs()
		catchParserErrors = False if self.catchParserErrors is None else True
		singleWorkflows = [parsedFileFlow.ParsedFileWorkflow(x,catchParserErrors=catchParserErrors) for x in calcObjs]
		compositeWorkflow = baseFlow.StandardLabelledWorkflowComposite(singleWorkflows)
		return compositeWorkflow

	def _getAllCalcObjs(self):
		return [x.create() for x in self._getAllCreators()]

	def _getAllCreators(self):
		outCreators = list()
		for geom, fileName in it.zip_longest(self.geoms, self._getFileNamesForAll()):
			currCreator = copy.deepcopy(self.baseCreator)
			currCreator.workFolder = self.outFolder
			currCreator.geom = geom
			currCreator.fileName = fileName
			outCreators.append(currCreator)
		return outCreators

	def _getFileNamesForAll(self):
		if self.fileNames is not None:
			return self.fileNames
		baseFileName = "inp_file" if self.baseFileName is None else self.baseFileName
		return [baseFileName + "_{}".format(idx) for idx,unused in enumerate(self.geoms,1)]

def getBasicOutDictsFromParsedFileMultiGeomsStdInpCreator(stdInpCreator, calcObjCreatorToDict=None):
	""" Get a list of basic dictionaries containing relevant options for this stdInpCreator; assuming jobs have already been run
	
	Args:
		stdInpCreator: (job_helpers/parsed_out_files)
		calcObjCreatorToDict: (Optional, function) Function which takes a CP2K object creator and returns a dictionary representation

	Returns
		outDicts: (iter of dicts) One dict per geometry in stdInpCreator
	
	"""
	creatorToDict = calcObjCreatorToDict if calcObjCreatorToDict is not None else creatorToDictHelp.getSimpleCreatorObjToDictMapObj()
	outDicts = list()	
	calcObjCreatorDict = creatorToDict(stdInpCreator.baseCreator)
	stdInpCreatorDict = {"compound":stdInpCreator.eleKey, "method":stdInpCreator.methodKey, "structure":stdInpCreator.structKey}
	currStdOut = stdInpCreator.create().createOutputObj()
	assert len(currStdOut.data)==1
	for currData in currStdOut.data[0]:
		currDict = dict()
		currDict["out_geom"] = currData.parsedFile.unitCell.toDict()
		currDict["energies"] = currData.parsedFile.energies.toDict()
		currDict.update(calcObjCreatorDict)
		currDict.update(stdInpCreatorDict)
		outDicts.append(currDict)
	return outDicts

