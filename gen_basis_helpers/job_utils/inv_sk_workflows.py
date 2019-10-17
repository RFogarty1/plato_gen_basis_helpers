
import copy
import itertools as it
import os
import pathlib
from types import SimpleNamespace

import plato_fit_integrals.core.workflow_coordinator as WFlowCoordinate
import plato_pylib.plato.mod_plato_inp_files as modInp
import plato_pylib.utils.job_running_functs as jobRun
import plato_pylib.plato.parse_inv_sk as parseInvSk



class CompositeInvSkWorkFlow(WFlowCoordinate.WorkFlowBase):
	""" WorkFlow combining a group of InvSkWorkFlow objects """

	def __init__(self, invSkWorkFlows, label):
		""" Create an Inverse-SK Composite workFlow (same essential interface parts as Inverse-SK WorkFlow)
		
		Args:
			invSkWorkFlows: iter, each element is a single InvSkWorkFlow object
			label: str, A way of identifying workFlow composites
		"""

		self.workFlows = invSkWorkFlows
		self.label = label

	@property
	def preRunShellComms(self):
		outList = list()
		for x in self.workFlows:
			currShellComms = x.preRunShellComms
			if currShellComms is not None:
				outList.extend(currShellComms)
		return outList

	@property
	def allLabels(self):
		outList = [self.label]
		for x in self.branches:
			outList.extend(x.allLabels)
		return outList

	@property
	def workFolder(self):
		return None

	@property
	def branches(self):
		return self.workFlows


	@property
	def namespaceAttrs(self):
		allAttrs = list()
		for x in self.workFlows:
			allAttrs.extend(x.namespaceAttrs)
		uniqueAttrs = list(set(allAttrs))

		return sorted(uniqueAttrs) #Want ordering to be consistent across calls; use of set removes this


	def fetch(self,*args):
		""" Traverses the tree structure to get the object specified by *args (each value of arg represents the label for one level down)
		
		Args:
			Each arg represents the label for a Composite workFlow. E.g. fetch("Zr","perfect_crystals","hcp") would first look
		for the Zr label at the top level; once found it would look for the perfect_crystals label in the next level of
		composites and then the hcp label in those. fetch("Zr","perfect_crystals") would also work fine.
				
		Returns
			InvSkWorkFlow(Composite or otherwise)
		
		Raises:
			Errors
		"""
		nArgs = len(args)
		if nArgs == 0:
			return self

		for b in self.branches:
			if b.label == args[0]:
				if len(args)==1:
					return b
				else:
					return b.fetch(*args[1:])

		raise ValueError("Could not find branch with label {}".format(args[0]))

		return None



	#NOTE: Not tested because its hard to test.....
	def run(self):
		#Step 1 = actually run the stuff
		for x in self.workFlows:
			x.run()

		#Step 2 - get a copy of each ouput
		allOutputs = list()
		for x in self.workFlows:
			allOutputs.append( copy.deepcopy(x.output) )

		#Step 3 = get a merged output
		mergedOutput = SimpleNamespace()
		allEleCombos = self.namespaceAttrs
		for eleCombo in allEleCombos:
			currBaseInvObj = None
			if len(allOutputs)==1:
				currBaseInvObj = copy.deepcopy( getattr(allOutputs[0],eleCombo) )

			for x in allOutputs:
				try:
					currAppendObj = copy.deepcopy (getattr(x,eleCombo))
					if currBaseInvObj is None:
						currBaseInvObj =  currAppendObj
					else:
						currBaseInvObj.addInvSKParsedFileData(currAppendObj)
				except AttributeError: #e.g. if mergeing Zr_Zr and Mg_Mg results, the Zr_Zr will be missing from the Mg_Mg parts
					pass
			setattr(mergedOutput, eleCombo, currBaseInvObj)

		self.output = mergedOutput


class InvSkWorkFlow(WFlowCoordinate.WorkFlowBase):
	""" WorkFlow for running inverse-SK calculations """

	def __init__(self, workFolder, optDict, elements, structs, label, removeXtal=True, units="ryd"):
		""" Create an Inverse-SK workFlow
		
		Args:
			workFolder: directory to run calculations in
			optDict: dict, options for running calcs in plato_pylib format
			elements: elements we want the inverse-SK data for
			structs: list of plato_pylib UnitCell structures
			label: Str, a way of identifying workFlow components
			removeXtal: Bool, whether to remove crystal-field terms from the parsed results
			units: units for the output energies; ryd or ev
		"""
		self._workFolder = os.path.abspath(workFolder)
		self.optDict = {k.lower():v for k,v in dict(optDict).items()}
		self.optDict["inversesk"] = 1
		self.structs = structs
		self.elements = list(elements)
		self._platoComm = "dft2" #only code that does inverse-SK calculations
		self.label = label
		self.removeXtal = removeXtal
		self.units=units

		#Only need to write the input files upon initialisation
		pathlib.Path(self.workFolder).mkdir(parents=True,exist_ok=True)
		self._writeFiles()

	@property
	def preRunShellComms(self):
		#TODO: Need to remove any previous temporary directories in a safe way before starting
		inpPaths = [x + ".in" for x in self._baseFilePaths]
		return jobRun.invSkInputPathsToBashComms(inpPaths)

	@property
	def workFolder(self):
		return self._workFolder

	@property
	def allLabels(self):
		return [self.label]

	@property
	def namespaceAttrs(self):
		allElementCombos = [x for x in it.product(self.elements,self.elements)]
		prefix = "invSkResults_"
		outSuffixes = ["{}_{}".format(*x) for x in allElementCombos]
		outAttrs = [prefix+suffix for suffix in outSuffixes]
		return outAttrs

	@property
	def _inpFilePaths(self):
		return [x + ".in" for x in self._baseFilePaths]

	@property
	def _baseFilePaths(self):
		basePaths = [os.path.join(self._workFolder,x) for x in self._baseFileNames]
		return basePaths

	@property
	def _baseFileNames(self):
		numbFiles = len(self.structs)
		allFileNames = list()
		for idx in range(numbFiles):
			outName = "inv_sk_inp_{}".format(idx)
			allFileNames.append(outName)
		return allFileNames

	@property
	def branches(self):
		return None


	def _writeFiles(self):
		inpPaths = self._inpFilePaths
		for inpPath,struct in it.zip_longest(inpPaths,self.structs):
			geomSection = modInp.getPlatoGeomDictFromUnitCell(struct)
			strDict = modInp.getStrDictFromOptDict(self.optDict, self._platoComm)
			strDict.update(geomSection)
			modInp.writePlatoOutFileFromDict(inpPath,strDict)

	#NOTE: not covered by a unit-test, because it relies so much on external functions
	def run(self):
		tempDict = dict()
		for comboStr in self.namespaceAttrs:
			eleCombo = comboStr.split("_")[1:]
			allFiles = self.getOutFilePathsOneEleCombo(*eleCombo)
			allParsed = [parseInvSk.parseInvSK(x, units=self.units) for x in allFiles]
			for x in allParsed:
				if self.removeXtal:
					x.removeXtalFieldTerms()
			combined = _getCombinedSetOfParsedInvSkObjsForSameEleCombo(allParsed)
			tempDict[comboStr] = combined
		self.output = SimpleNamespace(**tempDict)

	def getOutFilePathsOneEleCombo(self,eleA,eleB):
		eleStr = "{}_{}".format(eleA.capitalize(),eleB.capitalize())
		return [x+"_"+eleStr+"_SK.csv" for x in self._baseFilePaths]


def _getCombinedSetOfParsedInvSkObjsForSameEleCombo(allParsedFiles:list):
	outObj = allParsedFiles[0]
	for x in allParsedFiles[1:]:
		outObj.addInvSKParsedFileData(x)
	return outObj


