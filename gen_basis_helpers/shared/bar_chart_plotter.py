

import itertools as it
import numpy as np
import matplotlib.pyplot as plt
from . import data_plot_base as dPlotBase
from . import misc_utils as misc

class LabelledBarChartPlotterStandard(dPlotBase.DataPlotterBase):
	registeredKwargs = set(dPlotBase.DataPlotterBase.registeredKwargs)
	registeredKwargs.discard("plotFunct")
	registeredKwargs.add("fixDataSeriesOrLabels")
	registeredKwargs.add("labels") 
	registeredKwargs.add("labelRotations")
	registeredKwargs.add("barWidths")
	registeredKwargs.add("barSeparationWithinLabels")
	registeredKwargs.add("barSeparationBetweenLabels")
	registeredKwargs.add("xStartPos")
	registeredKwargs.add("patterns")

	def createPlot(self, plotData=None, **kwargs):
		""" Takes data in plotData argument and creates a plot 
		
		Args:
			plotData: list of input data. Each list entry is either one data series or one labels data (made from one set of y-values). E.g. [ [1,2], [3,4] ]. See attribute fixDataSeriesOrLabels for more information
			kwargs: keyword-arguments in same format as when assigning to the object attributres (keys in self.registeredKwargs). These are assigned to the object for solely this function call
		Returns
			Handle to the overall figure (if axHandle not supplied) or None (if axHandle supplied)
		"""
		return super().createPlot(plotData=plotData, **kwargs)

	#Overidden base class method
	def _setDefaultInitAttrs(self):
		self.fixDataSeriesOrLabels = "labels"
		self.barWidths = 1
		self.barSeparationWithinLabels = 0
		self.barSeparationBetweenLabels = 1
		self.xStartPos = 1


	@property
	def labels(self):
		""" iter of strs. These are the labels that are placed along the x-axis
		"""
		return self._labels

	@labels.setter
	def labels(self, val):
		self._labels = val


	@property
	def fixDataSeriesOrLabels(self):
		""" Property determining the meaning of one input data series. This can be set to either "dataSeries" or "labels". 
		"""
		if self._fixDataSeriesOrLabels is None:
			return "labels"
		else:
			return self._fixDataSeriesOrLabels

	@fixDataSeriesOrLabels.setter
	def fixDataSeriesOrLabels(self,val):
		if val is None:
			self._fixDataSeriesOrLabels = None
			return None
		allowedVals = ["dataSeries".lower(), "labels".lower()]
		if val.lower() not in allowedVals:
			raise AttributeError("{} is not an allowed value for fixDataSeriesOrLabels; allowed values are {}".format(val, allowedVals))
		self._fixDataSeriesOrLabels = val.lower()


	#OVERIDING functions on the parent class
	def _getToPlotData(self, plotData=None):
		inpFormatData = self._getInputFormatToPlotData(plotData)
		if inpFormatData is None:
			return list()
		return self.mapInputFormatToPlotFormat(inpFormatData)

	def _getInputFormatToPlotData(self, plotData=None):
		if plotData is not None:
			toPlot = plotData
		else:
			if self.data is not None:
				toPlot = self.data
			else:
				toPlot = None
		return toPlot

	def mapInputFormatToPlotFormat(self, inpData):
		with misc.temporarilySetInstanceAttrs(self,{"data":inpData}):
			allXVals = self._getXVals()
			allYVals = self._getYValsInFixedLabelFormat()

		outSeries = [list() for x in range(len(allXVals[0]))]

		#Get the values
		for xVals,yVals in it.zip_longest(allXVals,allYVals): #1 label at a time?
			assert len(xVals)==len(yVals)
			for idxDataSeries, unused in enumerate(xVals):
				outSeries[idxDataSeries].append(   [xVals[idxDataSeries], yVals[idxDataSeries]] )

		return [np.array(x) for x in outSeries]



	def _getXVals(self):
		yValsAll = self._getYValsInFixedLabelFormat()
		outXVals = list()
		currXPos = self.xStartPos + (0.5*self.barWidths)
		for yVals in yValsAll:
			currXVals = list()
			#Get the x-values in the centre of each bar for THIS label
			for yVal in yVals:
				currXVals.append( currXPos )
				currXPos += self.barWidths+self.barSeparationWithinLabels
			outXVals.append(currXVals)
			#Now get to the centre of the bar for the NEXT data series
			currXPos += (self.barSeparationBetweenLabels - self.barSeparationWithinLabels)
		return outXVals


	def _getYValsInFixedLabelFormat(self):
		if self.fixDataSeriesOrLabels == "labels":
			return self.data
		elif self.fixDataSeriesOrLabels == "dataSeries".lower():
			return self._getYValsFromFixedDataFormat()
		else:
			raise ValueError("Hit a part of an elif loop i didnt think was possible")

	def _getYValsFromFixedDataFormat(self):
		numbLabels = len(self.data[0])
		assert all([len(x)==numbLabels for x in self.data])

		outYVals = [list() for x in range(numbLabels)]
		for yVals in self.data:
			for idx,yVal in enumerate(yVals):
				outYVals[idx].append(yVal)

		return outYVals


	def _plotSingleDataSeries(self, inpData, label, idx):
		if self.patterns is None:
			patternKwarg = None
		else:
			patternKwarg = misc.getCycleListToMaxLength(self.patterns,idx+1)[idx]

		plt.bar(inpData[:,0], inpData[:,1], self.barWidths, label=label, hatch=patternKwarg)


	#Overrides parent hook function
	def _modifyAxisPlot(self):
		plotData = self.data #gauranteed by parent class
		self._putLabelsOnPlot(plotData)

	def _putLabelsOnPlot(self, plotData):
		if plotData == list():
			return None

		#Get data in label-centrice format
		xTickPositions = list()
		labelCentricXVals = [list() for x in range(len(plotData[0][:,0]))]
		for dSeries in plotData:
			for idx,xVal in enumerate(dSeries[:,0]):
				labelCentricXVals[idx].append( [xVal, 0] ) #y-value is needed, but doesnt matter what it is 

		labelCentricXVals = [np.array(x) for x in labelCentricXVals]

		#Modify the tick positions and then the labels
		xTickPositions = _getXTickPositionsForLabels(labelCentricXVals)
		plt.xticks(xTickPositions, self.labels, rotation=self.labelRotations)	


def _getXTickPositionsForLabels(inpData):
	outVals = list()
	for inpVals in inpData:
		sumXVals = sum([x for x in inpVals[:,0] ])
		averageX = sumXVals / len(inpVals[:,0])
		outVals.append(averageX)
	return outVals
