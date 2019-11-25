
import itertools as it
import numpy as np
import matplotlib.pyplot as plt

from ..shared import data_plot_base as basePlotter
from ..shared import misc_utils as misc


class EosEnergyDataPlotter(basePlotter.DataPlotterStandard):

	def __init__(self, **kwargs):
		self.registeredKwargs.add("sortXBeforePlot")
		self.registeredKwargs.add("lineMarkers")
		self.registeredKwargs.add("lineMarkerSizes")
		self.registeredKwargs.add("lineMarkerFillStyles")
		super().__init__(**kwargs)

		self.methodProps = ["lineStyles", "lineMarkers", "lineMarkerSizes", "dataLabels","lineMarkerFillStyles"]
		self.dataSeriesProps = ["lineColors"]



	@classmethod
	def fromDefaultPlusKwargs(cls, **kwargs):
		inpKwargs = dict()
		inpKwargs["xlabel"] = "Volume per atom ($a_0^3$)"
		inpKwargs["sortXBeforePlot"] = True
		inpKwargs.update(kwargs)
		return cls(**inpKwargs)


	#NOTE: Almost a total duplicate of DataPlotterDiagMatrixEles
	def createPlot(self, plotData, **kwargs):
		""" Takes data in plotData argument and creates a plot 
		
		Args:
			plotData: list of input data. Each list entry is data for one method. i.e. plotData=[methodAData, methodBData]. In turn these are lists of numpy arrays,
			          with xData in col1 and ydata in col2
			kwargs: keyword-arguments in same format as when assigning to the object attributres (keys in self.registeredKwargs). These are assigned to the object for solely this function call
		Returns
			Handle to the overall figure
		"""
		toPlot = self._getDataFormattedForSuperPlotter(plotData)

		with misc.fragile(basePlotter.temporarilySetDataPlotterRegisteredAttrs(self,kwargs)):

			for prop in self.methodProps:
				currVal = getattr(self,prop)
				if currVal is not None:
					setattr(self, prop, self._getMethodBasedArgListInCorrectFormat(plotData,currVal))

			for prop in self.dataSeriesProps:
				currVal = getattr(self,prop)
				if currVal is not None:
					setattr(self, prop, self._getDataSeriesBasedArgListInCorrectFormat(plotData, currVal))


			outFig = super().createPlot(toPlot)
			self.changeLineProp(outFig, "lineMarkers", lambda inpLine,value:inpLine.set_marker(value))
			self.changeLineProp(outFig, "lineMarkerSizes", lambda inpLine,value:inpLine.set_markersize(value), inclLegend=False)
			self.changeLineProp(outFig, "lineMarkerFillStyles", lambda inpLine,value:inpLine.set_fillstyle(value), inclLegend=False)


			if self.legend:
				currAx = outFig.get_axes()[0]
				currAx.legend()
				self.changeLegendEntriesToMethodOnly(outFig)


		return outFig



	def _getDataFormattedForSuperPlotter(self,plotData):
		outData = list()
		for methData in plotData:
			outData.extend(methData)

		if self.sortXBeforePlot:
			for idx,currData in enumerate(outData):
				sortedData = currData[ currData[:,0].argsort() ]
				outData[idx] = sortedData

		return outData


	def _getDataSeriesBasedArgListInCorrectFormat(self, plotData, propInInputFormat):
		outList = list()
		numbDataSeries = self._getNumbDataSeries(plotData)
		inpProps = misc.getCycleListToMaxLength(propInInputFormat, numbDataSeries)

		for mIdx,methData in enumerate(plotData):
			for cIdx, dSerieData in enumerate(methData):
				outList.append( inpProps[cIdx] )
		return outList


	def _getMethodBasedArgListInCorrectFormat(self, plotData, propInInputFormat):
		outList = list()
		numbDataSeries = self._getNumbDataSeries(plotData)
		inpProps = misc.getCycleListToMaxLength(propInInputFormat, numbDataSeries)

		for mIdx, methData in enumerate(plotData):
			for cIdx, dSerieData in enumerate(methData):
				outList.append( inpProps[mIdx] )
		return outList

	def _getNumbDataSeries(self, plotData):
		count = 0
		for methData in plotData:
			count += len(methData)
		return count


	#DUPLICATED FROM DataPlotterDiagMatrixEles
	def changeLegendEntriesToMethodOnly(self, outFig):
		currLegend = outFig.get_axes()[0].get_legend()
		usefulLines, usefulText = list(), list()
		currAx = outFig.get_axes()[0]

		#We grab the first instance of a new method
		for handle,textObj in it.zip_longest(currAx.get_lines(), currLegend.texts):
			if textObj.get_text() not in usefulText:
				usefulText.append( textObj.get_text() )
				usefulLines.append( handle )

		plt.legend(usefulLines,usefulText)

