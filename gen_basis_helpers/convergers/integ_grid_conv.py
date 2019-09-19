

""" Set of objects/functions for helping converge integration grid spacing in plato """

import copy
import os

from . import convergers as convObjs

def createDft2AngularSpacingOptDictMapper():
	return VaryDft2AngularSpacingParamToOptDictMapper()

def createDft2RadialSpacingOptDictMapper():
	return VaryDft2RadialSpacingParamToOptDictMapper()

def createDftFFTSpacingOptDictMapper():
	return VaryFFTSpacingOptDictMapper()

def createKPointsMPGridOptDictMapper():
	return VaryKPtMPGridOptDictMapper()


class VaryDft2AngularSpacingParamToOptDictMapper(convObjs.VaryParamToOptDictMapper):

	def __init__(self):
		pass

	def modOptDictWithVariableParam(self, optDict, varParam):
		_ensureDft2MeshTypeIsAtomCentred(optDict)
		radialParam = optDict["integralmeshspacing"][0]
		optDict["integralmeshspacing"] = [radialParam, varParam, varParam]


class VaryFFTSpacingOptDictMapper(convObjs.VaryParamToOptDictMapper):

	def __init__(self):
		pass

	def modOptDictWithVariableParam(self, optDict, varParam):
		optDict["fftGridSpacing".lower()] = varParam


class VaryDft2RadialSpacingParamToOptDictMapper(convObjs.VaryParamToOptDictMapper):

	def __init__(self):
		pass

	def modOptDictWithVariableParam(self, optDict, varParam):
		_ensureDft2MeshTypeIsAtomCentred(optDict)
		outGrid = copy.deepcopy( optDict["integralmeshspacing"] )
		outGrid[0] = varParam
		optDict["integralmeshspacing"] = outGrid


class VaryKPtMPGridOptDictMapper(convObjs.VaryParamToOptDictMapper):

	def __init__(self):
		pass

	def modOptDictWithVariableParam(self, optDict, varParam):
		optDict["blochstates"] = varParam


def _ensureDft2MeshTypeIsAtomCentred(optDict):
	try:
		assert optDict["integralmeshtype"]=="atom", "integralmeshtype must be equal to atom" 
	except KeyError:
		raise AssertionError( "integralmeshtype key missing in optDict (Needs to be present and set to atom)" )


#Possibly overly-tied to the (non-enforced) interface of CreateStructEnergiesWorkFlow
#class AtomCentreGridJobRunner(convObjs.PropConvJobRunner):
#
#	def __init__(self, structSetLabel, workFlowFactory, variableParams, variableToGridInfoObj):
#		self.variableParams = variableParams
#		self.variableToGridInfoObj = variableToGridInfoObj
#		self.workFlowFactory = workFlow
##		self._defWorkFolder = os.path.join(self._baseFolder, endFolderName)
##		self._workFolder = None
##		self._baseFolder = os.path.abspath(baseFolder)
##		endFolderName = "convgrid_{}_{}".format(variableToGridInfoObj.varyType, structSetLabel)
#
#
#	@property
#	def label(self):
#		return None
#
#
#	@property
#	def fileNames(self):
#		baseName = "acg"
#		gridTypeName = self.variableToGridInfoObj.varyType
#		paramContrib = [int(x) for x in self.variableParams]
#		outNames = ["{}_{}_{}.in".format(baseName, gridTypeName, x) for x in paramContrib]
#		return outNames
#
#	@property
#	def runComms(self):
#		allComms = list()
#		for x in self._createWorkFlows():
#			allComms.extend(x.preRunShellComms)
#		return allComms
#
#
#	def _createWorkFlows():
#		outWorkFlows = list()
#		for vParam in varParams:
#			currFactory = copy.deepcopy(self.workFlowFactory)
#			currGridParams = self.variableToGridInfoObj.getGridFromVariable(vParam)
#			baseFactory.modOptsDict[self.variableToGridInfoObj.gridKwarg] = currGridParams
#			currFactory = os.path.join(self.workFolder, "val_{}".format(vParam))
#			currWorkFlow = currFactory()
#			outWorkFlows.append(currWorkFlow)
#		return outWorkFlows
#
#
#	def _parseDataAfterRunningJobs(self):
#		pass
##		for x in self._createWorkFlows():
#			
#	def createAnalyser(self, dataPlotter=None):
#		return  GeneralConvPropAnalyser(self.label, self._parseDataAfterRunningJobs(), dataPlotter=dataPlotter)
#
#
#	@property
#	def workFolder(self):
#		if self._workFolder is None:
#			return self._defWorkFolder
#		else:
#			return self._workFolder
#
#	@workFolder.setter
#	def workFolder(self, value):
#		self._workFolder = value
#
#
#
#class GeneralConvPropAnalyser(convObjs.PropConvAnalyserComposite):
#
#	def __init__(self, label, data, dataPlotter = None):
#		self.label = label 
#		self.data = data
#		self.dataPlotter = dataPlotter
#
#	@property
#	def data(self):
#		return data
#
#	@property
#	def label(self):
#		return label
#

