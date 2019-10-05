
import numpy as np
import types

import plato_pylib.plato.parse_tbint_files as parseTbint

from . import data_plotter_matrix_eles as dPlotter
from ..shared.label_objs import StandardLabel
from ..job_utils import matrix_elements_helpers as matHelp
from ..shared import misc_utils as misc

def swapEosRunMethodToExtractDiagMatrixElements(workFlow):
	replaceRunMethodInWorkFlow(workFlow, volVsDiagMatrixElementsFromWorkFlow)


def replaceRunMethodInWorkFlow(workFlow, newRunMethod):
	workFlow.run = types.MethodType(newRunMethod, workFlow)


#If you dont want volPerAtom you can probably use functtools to wrap it
def volVsDiagMatrixElementsFromWorkFlow(inpWorkFlow, volPerAtom=True):
	inpFileList = inpWorkFlow._inpFilePaths
	outVals = matHelp.getVolVsOnSiteDiagTerms(inpFileList, volPerAtom=volPerAtom)

	try:
		startNamespace = getattr(inpWorkFlow,"output")
	except AttributeError:
		startNamespace = types.SimpleNamespace()
	startNamespace.volVsDiagOnSite=outVals
	inpWorkFlow.output = startNamespace

def createMatrixEleAnalysersFromWFlowCoord(wFlowCoord):
	analysers = list()
	for x in wFlowCoord._workFlows:
		analysers.append( createMatrixEleAnalyserFromWorkFlow(x)  )
	return analysers

def createMatrixEleAnalyserFromWorkFlow(workFlow):
	rawOutput = workFlow.output.volVsDiagOnSite
	forNp = list()
	for volData in rawOutput:
		currData = list([volData[0]])
		for x in volData[1]:
			currData.extend(x)
		forNp.append(currData)
	outData = np.array(forNp)

	return EleEosMatrixDiagElementsResult(outData, eleKey=workFlow.element,structKey=workFlow.structKey,methodKey=workFlow.methodKey)

class EleEosMatrixDiagElementsComposite():

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=True)
	def __init__(self, inpObjs):
		self.objs = list(inpObjs)
		self._ensureLabelsAllDifferent()

	def _ensureLabelsAllDifferent(self):
		if len(list(self.label)) != len(set(self.label)):
			raise ValueError("Duplicate labels detected")

	@property
	def label(self):
		outList = list()
		for x in self.objs:
			outList.extend(x.label)

		return outList

	def getPlotData(self):
		outList = list()
		for x in self.objs:
			outList.extend(x.getPlotData())
		return outList

class EleEosMatrixDiagElementsResult():
	""" Class for holding data related to diagonal matrix elements for structures used to calculate eos

	Attributes:
		label (1-ele list): Contains StandardLabel obj, used to differentiate objects

	"""
	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=False)
	def __init__(self, inpData, eleKey=None, structKey=None, methodKey=None):
		""" Initializer
		
		Args:
			inpData (nxm iter): np.array(inpData) called. First column should be volumes, all other columns energies
			                    for matrix elements
			eleKey(required keyword, str): Element label. Used with structKey/methodKey to uniquely identify/search for this object
			structKey(required keyword, str): structure label
			methodKey(required keyword, str): method label

		"""
		self._label = StandardLabel(eleKey=eleKey, structKey=structKey,methodKey=methodKey)
		self._data = np.array(inpData)


	@property
	def label(self):
		return [self._label]

	def getPlotData(self):
		return [np.array(self._data)]


def createShellMapperFromAdtFile(inpFile):
	angMomToLabel = {0:"s", 1:"p", 2:"d"}
	parsedFile = parseTbint.parseAdtFile(inpFile)
	shellToAngMomMap = parsedFile["shellToAngMom"]

	shellOrder = list()
	for key in sorted(shellToAngMomMap.keys()):
		currShellAngMom = shellToAngMomMap[key]
		currShellLabel = angMomToLabel[currShellAngMom]
		shellOrder.append( currShellLabel )

	outObj = ShellMapper.fromShellOrderings(shellOrder)

	return outObj

class ShellMapper():
	"""Class used to map shell angular mom+index to a list of orbital indices

	"""

	_angLabelToNumbOrbs = {"s":1, "p":3, "d":5}
	def __init__(self, sStart=None, pStart=None, dStart=None):
		""" Initialiser. Generally should not be called directly
		
		Args:
			sStart(iter): List of indices corresponding to the first orbital index for the
			              nth shell with s-type angular momentum (orbAngMom=0) 
			pStart(iter): As above but for p-shells (orbAngMom=1)
			dStart(iter): As above but for d-shells (orbAngMom=2)	
		"""
		self.sStart = list(sStart) if sStart is not None else None
		self.pStart = list(pStart) if pStart is not None else None
		self.dStart = list(dStart) if dStart is not None else None
		self._angLabelToShellStartList = {"s":self.sStart, "p":self.pStart, "d":self.dStart}

		self._eqAttrs = ["sStart","pStart","dStart"]


	@classmethod
	def fromShellOrderings(cls, shellOrder):
		""" Alternative intialiser
		
		Args:
			shellOrder(iter): List of strings denoting the ordering of shels in terms of angular momenta of orbitals. e.g. ["s","p","d"] 
				
		"""
		sStart, pStart, dStart = list(), list(), list()
		_angLabelToList = {"s":sStart, "p":pStart, "d":dStart}
		total = 0
		for key in shellOrder:
			currList = _angLabelToList[key]
			currList.append( total )
			total += cls._angLabelToNumbOrbs[key]

		return cls(sStart=sStart, pStart=pStart, dStart=dStart)


	@property
	def numbOrbs(self):
		total = 0
		total += sum([1 for x in self.sStart])
		total += sum([3 for x in self.pStart])
		total += sum([5 for x in self.dStart])
		return total


	def getDiagIndicesNthShellOfType(self, angMomLabel, nthShell=0):
		""" Returns the orbital-indices for a user-specified shell
		
		Args:
			angMomLabel (str): s,p or d label for the shell angular momentum
			nthShell (opt,int): Which shell (starting at zero) you want for that angular momentum. E.g. values of 0 or 1 make sense for a double-zeta basis
				
		Returns
			diagIndices(int iter): List of the orbital indices for requested shell.
		
		Raises:
			Errors
		"""
		startIdx = self._angLabelToShellStartList[angMomLabel][nthShell]
		endIdx = startIdx + self._angLabelToNumbOrbs[angMomLabel]
		return [x for x in range(startIdx,endIdx)]


	def __eq__(self, other):
		if not isinstance(other,ShellMapper):
			return False

		for attr in self._eqAttrs:
			if getattr(self,attr) != getattr(other,attr):
				return False

		return True

class MatEleEosDataPresenter():

	def __init__(self, analyser, shellMapper, dataPlotter=None):
		""" Initialiser
		
		Args:
			analyser (EleEosMatrixDiagElementsComposite): Object with getPlotData and composite search implemented in same way as EleEosMatrixDiagElementsComposite
			shellMapper (simpleNamespace): attrs are eleKeys, values are ShellMapper objects for each element
			dataPlotter: Object that controls how we plot the data. Generally should leave the default on and just modify its options
		"""
		self.analyser = analyser
		self.shellMapper = shellMapper
		if dataPlotter is None:
			self.dataPlotter = dPlotter.DataPlotterDiagMatrixEles.fromDefaultPlusKwargs()


	#WARNING: NOT UNIT-TESTED
	#TODO: We need to have titleStr sorted out for this to really work for more than 1 ele/struct at once
	def plotData(self, shellOrdering, eleKeys=None, structKeys=None, methodKeys=None, plotRelE=False, refMethod=None):
		outFigs = list()
		if eleKeys is None:
			eleKeys = set([x.eleKey for x in self.analyser.label])
		if structKeys is None:
			structKeys = set([x.structKey for x in self.analyser.label])
		if methodKeys is None:
			methodKeys = set([x.methodKey for x in self.analyser.label])


		for eKey in eleKeys:
			for sKey in structKeys:
				print("eKey={}".format(eKey))
				print("sKey={}".format(sKey))
				currData = self.getMappedPlotDataOneEleAndStruct(methodKeys, sKey, eKey, shellOrdering, subRefData=plotRelE, refMethod=refMethod)
				currPlot = self.dataPlotter.createPlot(currData)
				outFigs.append(currPlot)
		return outFigs


	def getMappedPlotDataOneEleAndStruct(self, methodKeys, structKey, eleKey, shellOrdering, subRefData=False, refMethod=None):
		startDataDict = self.getPlotDataWithoutMapping(methodKeys, structKey, eleKey)
		mapper = getattr(self.shellMapper,eleKey)
		mappedDataDict = dict()
		for key in methodKeys:
			mappedDataDict[key] = self._getShellMappedVersionOfArray( startDataDict[key], shellOrdering, mapper)

		if subRefData:
			refData = np.array( mappedDataDict[refMethod] )
			for key in methodKeys:
				mappedDataDict[key][:,1:] -= refData[:,1:]

		#Output order is the same as input order
		outList = list()
		for x in methodKeys:
			outList.append( mappedDataDict[x] )

		return outList

	def getPlotDataWithoutMapping(self, methodKeys:iter, structKey, eleKey):
		outDict = dict()
		for method in methodKeys:
			outDict[method] = self._getSingleMethodPlotDataWithoutMapping(method, structKey, eleKey)
		return outDict

	def _getSingleMethodPlotDataWithoutMapping(self, methodKey, structKey, eleKey):
		relevantObjs = self.analyser.getObjectsWithComponents([methodKey, structKey, eleKey],caseSensitive=False)
		assert len(relevantObjs)==1, "{} objects (instead of 1) found with components method={},struct={},eleKey={}".format(len(relevantObjs), methodKey,structKey,eleKey)
		outData = relevantObjs[0].getPlotData()


		assert len(outData) == 1, "number of output arrays needs to be 1, not {}".format( len(outData) )

		return outData[0]


	def _getShellMappedVersionOfArray(self, inpArray, shellOrdering, shellMapper):

		outArray = np.array( [inpArray[:,0]] ).T #This gets us the volumes at least
#		outArray = np.array( inpArray[:,0] ).T #This gets us the volumes at least

		for sLabel in shellOrdering:
			currCols = [1+x for x in shellMapper.getDiagIndicesNthShellOfType(sLabel)]
			arraySlice = inpArray[:,currCols]
			outArray = np.concatenate( [outArray, arraySlice], axis=1 )

		return outArray




