
import itertools as it

import numpy as np

import matplotlib.pyplot as plt

from ..shared import data_plot_base as basePlotter
from ..shared import misc_utils as misc


class DataPlotterDiagMatrixEles(basePlotter.DataPlotterStandard):

	#Add extra Kwargs here
	def __init__(self, **kwargs):
		self.registeredKwargs.add("lineStyles") #TODO: Suspect this is redundant - check
		self.registeredKwargs.add("sortXBeforePlot")
		super().__init__(**kwargs)


	@classmethod
	def fromDefaultPlusKwargs(cls, **kwargs):
		inpKwargs = dict()
		inpKwargs["xlabel"] = "Volume / bohr^{3}"
		inpKwargs["sortXBeforePlot"] = True
		inpKwargs.update(kwargs)
		return cls(**inpKwargs)

	def createPlot(self, plotData, **kwargs):
		""" Takes data in plotData argument and creates a plot 
		
		Args:
			plotData: list of input data. Each list entry is data for one method. i.e. plotData=[methodAData, methodBData]. In turn these are single numpy arrays
                      with xData in column1, and yData in all other columns
			kwargs: keyword-arguments in same format as when assigning to the object attributres (keys in self.registeredKwargs). These are assigned to the object for solely this function call
		Returns
			Handle to the overall figure
		"""

		toPlot = self._getDataFormattedForSuperPlotter(plotData)

		#We need to temporarily change the format of self.lineStyles in order for the super() method to interpret/plot
		#them correctly
		with misc.fragile(basePlotter.temporarilySetDataPlotterRegisteredAttrs(self,kwargs)):
			if self.lineStyles is not None:
				self.lineStyles = self._getLineStylesFormattedForSuperPlotter(plotData, self.lineStyles)
			if self.lineColors is not None:
				self.lineColors = self._getLineColorsFormattedForSuperPlotter(plotData, self.lineColors)
			if self.dataLabels is not None:
				self.dataLabels = self._getDataLabelsFormattedForSuperPlotter(plotData, self.dataLabels)

			outFig = super().createPlot(toPlot) #Important not to pass any Kwargs, the context manager has already translated them into attributes
			if self.legend:
				self.changeLegendEntriesToMethodOnly(outFig)

		return outFig



	def _getDataFormattedForSuperPlotter(self, plotData):
		outData = list()
		for methodData in plotData:
			for cIdx in range(1,methodData.shape[1]):
				currData = np.array( [methodData[:,0], methodData[:,cIdx]] ).T
				if self.sortXBeforePlot:
					currData = currData[ currData[:,0].argsort() ] #Sorting by 1st column(x) values
				outData.append( currData )
		return outData

	def _getLineStylesFormattedForSuperPlotter(self, plotData, lineStyles):
		return self._getMethodBasedArgListInCorrectFormat(plotData, lineStyles)

	def _getLineColorsFormattedForSuperPlotter(self, plotData, lineColors):
		return self._getDataSeriesBasedArgListInCorrectFormat(plotData, lineColors)

	def _getDataLabelsFormattedForSuperPlotter(self, plotData, dataLabels):
		return self._getMethodBasedArgListInCorrectFormat(plotData, dataLabels)

	def changeLegendEntriesToMethodOnly(self, outFig):
		currLegend = outFig.get_axes()[0].get_legend()
		usefulLines, usefulText = list(), list()

		#We grab the first instance of a new method
		for handle,textObj in it.zip_longest(currLegend.legendHandles, currLegend.texts):
			if textObj.get_text() not in usefulText:
				usefulText.append( textObj.get_text() )
				usefulLines.append( handle )

		plt.legend(usefulLines,usefulText)
			


	def _getMethodBasedArgListInCorrectFormat(self, plotData, propInInputFormat):
		outData = list()
		for idx,methodData in enumerate(plotData):
			for cIdx in range(1, methodData.shape[1]):
				outData.append( propInInputFormat[idx] )
		return outData

	def _getDataSeriesBasedArgListInCorrectFormat(self, plotData, propInInputFormat):
		outData = list()
		for mIdx, methodData in enumerate(plotData):
			for cIdx in range(1,methodData.shape[1]):
				outData.append( propInInputFormat[cIdx-1] )
		return outData



