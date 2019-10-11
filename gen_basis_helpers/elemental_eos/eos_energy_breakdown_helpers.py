

import functools
import os
import types
import numpy as np

from collections import OrderedDict

import plato_pylib.plato.parse_plato_out_files as platoOut

from . import base_objs as baseObjs
from .data_plotter_energy_breakdowns import EosEnergyDataPlotter 
from .matrix_ele_helpers import replaceRunMethodInWorkFlow
from ..shared.label_objs import StandardLabel
from ..shared import misc_utils as misc



def swapEosRunMethodToExtactE0Energies(workFlow, rydToEv=False):
	if rydToEv:
		swapFunct = functools.partial(volVsE0EnergiesFromWorkFlow, rydToEv=True)
	else:
		swapFunct = volVsE0EnergiesFromWorkFlow
	replaceRunMethodInWorkFlow(workFlow, swapFunct)


def volVsE0EnergiesFromWorkFlow(inpWorkFlow, perAtom=True, rydToEv=False):
	inpPaths = inpWorkFlow._inpFilePaths
	outPaths = [os.path.splitext(x)[0] + ".out" for x in inpPaths]

	if rydToEv:
		parser = platoOut.parsePlatoOutFile_energiesInEv
	else:
		parser = platoOut.parsePlatoOutFile


	allVols, allE0 = list(), list()
	for currPath in outPaths:
		parsedFile = parser(currPath)
		vol, e0 = parsedFile["unitCell"].volume, parsedFile["energies"].e0Coh
		if perAtom:
			nAtoms = parsedFile["numbAtoms"]
			vol, e0 = vol/nAtoms, e0/nAtoms
		allVols.append(vol)
		allE0.append(e0)

	outVals = np.array( [(vol,e0) for vol,e0 in zip(allVols,allE0)] )

	try:
		startNamespace = getattr(inpWorkFlow,"output")
	except AttributeError:
		startNamespace = types.SimpleNamespace()

	startNamespace.volVsE0 = outVals
	inpWorkFlow.output = startNamespace


def createAnalysersFromWorkFlowCoord(wFlowCoord):
	analysers = list()
	for x in wFlowCoord._workFlows:
		currData = x.output.volVsE0
		analysers.append( EosEnergyAnalyser( currData, eleKey=x.element, structKey=x.structKey, methodKey=x.methodKey ) )
	return analysers


class EosEnergyAnalyser(baseObjs.BaseEosPropAnalyser):

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=False)
	def __init__(self, inpData, eleKey=None, structKey=None, methodKey=None):
		""" Initializer
		
		Args:
			inpData (nx2 iter): np.array(inpData) called. First column should be volumes, second column is energies
			eleKey(required keyword, str): Element label. Used with structKey/methodKey to uniquely identify/search for this object
			structKey(required keyword, str): structure label
			methodKey(required keyword, str): method label

		"""
		self._label = StandardLabel(eleKey=eleKey, structKey=structKey, methodKey=methodKey)
		self._data = np.array( inpData )


	@property
	def label(self):
		return [self._label]

	def getPlotData(self):
		return [np.array(self._data)]


class EosEnergyPresenter():

	def __init__(self, analyser, dataPlotter=None):
		self.analyser = analyser
		if dataPlotter is None:
			self.dataPlotter = EosEnergyDataPlotter.fromDefaultPlusKwargs()
		else:
			self.dataPlotter = dataPlotter


	def plotData(self, eleKeys=None, structKeys=None, methodKeys=None, plotRelE=None, refMethod=None):
		if eleKeys is None:
			eleKeys = set([x.eleKey for x in self.analyser.label])
		if structKeys is None:
			structKeys = set([x.structKey for x in self.analyser.label])
		if methodKeys is None:
			methodKeys = set([x.methodKey for x in self.analyser.label])

		outFigs = list()
		for eKey in eleKeys:
			currFig = self._createFigOneElement(eKey, structKeys, methodKeys, plotRelE=plotRelE, refMethod=refMethod)
			outFigs.append(currFig)

		return outFigs



	def _createFigOneElement(self, eleKey, structKeys, methodKeys, plotRelE=None, refMethod=None):
		figData = self._getPlotDataOneEle(eleKey, structKeys,methodKeys, plotRelE=plotRelE, refMethod=refMethod)

		outFig = self.dataPlotter.createPlot(figData)
		return outFig



	def _getPlotDataOneEle(self, eleKey, structKeys, methodKeys, plotRelE=None, refMethod=None):
		dataDict = OrderedDict()
		for mKey in methodKeys:
			dataDict[mKey] = self._getRawPlotDataOneMethod(eleKey, structKeys, mKey)

		if (plotRelE is not None) :
			refData = self._getRawPlotDataOneMethod(eleKey, structKeys, refMethod)

		if plotRelE:
			for mKey in methodKeys:
				for structIdx, unused in enumerate(structKeys):
					self._subDataBFromA( dataDict[mKey][structIdx], refData[structIdx] )

		outList = list()
		for mKey in methodKeys:
			outList.append( dataDict[mKey] )

		return outList


	def _getRawPlotDataOneMethod(self, eleKey, structKeys, methodKey):
		outList = list()
		for sKey in structKeys:
			relObj = self.analyser.getObjectsWithComponents([eleKey,sKey,methodKey])
			assert len(relObj) == 1, "{} objects found with {},{},{} keys".format(len(relObj),eleKey,sKey,methodKey)
			currData = relObj[0].getPlotData()[0]
			outList.append(currData)
		return outList


	def _subDataBFromA(self,dataA, dataB):
		assert dataA.shape == dataB.shape
		assert np.allclose( dataA[:,0], dataB[:,0] )
		dataA[:,1] -= dataB[:,1]



