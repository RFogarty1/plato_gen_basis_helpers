
import copy
import itertools as it
from . import standard_template_obj as stdTemplate
from ..shared import calc_runners as calcRunners
from ..workflows import parsed_file_workflow as parsedFileFlow
from ..workflows import base_flow as baseFlow

class ParsedFileObjsForMultiGeomsStandardInputCreator(stdTemplate.StandardInputCreatorTemplateBase):

	registeredKwargs = set(stdTemplate.StandardInputCreatorTemplateBase.registeredKwargs)

	registeredKwargs.add("geoms")
	registeredKwargs.add("baseCreator")


	def _createFromSelf(self):
		workflow = self._getWorkflow()
		outObj = calcRunners.StandardInputObj(workflow, self.label)
		return outObj

	def _getWorkflow(self):
		calcObjs = self._getAllCalcObjs()
		singleWorkflows = [parsedFileFlow.ParsedFileWorkflow(x) for x in calcObjs]
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
		return ["inp_file_{}".format(idx) for idx,unused in enumerate(self.geoms,1)]




