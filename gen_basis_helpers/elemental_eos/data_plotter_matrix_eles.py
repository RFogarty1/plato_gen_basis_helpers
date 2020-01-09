
import itertools as it

import numpy as np

import matplotlib.pyplot as plt

from ..shared import data_plot_base as basePlotter
from ..shared import misc_utils as misc


class DataPlotterDiagMatrixEles(basePlotter.DataPlotterStandard):

	#Add extra Kwargs here
	def __init__(self, **kwargs):
		self.registeredKwargs.add("sortXBeforePlot")
		super().__init__(**kwargs)


	@classmethod
	def fromDefaultPlusKwargs(cls, **kwargs):
		inpKwargs = dict()
		inpKwargs["xlabel"] = "Volume ($a_0^3$ per atom)"
		inpKwargs["ylabel"] = "$\Delta E$ On-site matrix elements (eV)"
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

		#We need to temporarily change the format of self.lineStyles and similar
		with misc.fragile(basePlotter.temporarilySetDataPlotterRegisteredAttrs(self,kwargs)):

			methodProps = ["lineStyles","lineMarkers", "lineMarkerSizes", "dataLabels"]
			dataSeriesProps = ["lineColors"]

			for prop in methodProps:
				currVal = getattr(self,prop)
				if currVal is not None:
					setattr(self, prop, self._getMethodBasedArgListInCorrectFormat(plotData,currVal))

			for prop in dataSeriesProps:
				currVal = getattr(self,prop)
				if currVal is not None:
					setattr(self, prop, self._getDataSeriesBasedArgListInCorrectFormat(plotData, currVal))

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


	def _getMethodBasedArgListInCorrectFormat(self, plotData, propInInputFormat):
		outData = list()
		inpData = misc.getCycleListToMaxLength(propInInputFormat, len(plotData))
		for idx,methodData in enumerate(plotData):
			for cIdx in range(1, methodData.shape[1]):
				outData.append( inpData[idx] )
		return outData

	def _getDataSeriesBasedArgListInCorrectFormat(self, plotData, propInInputFormat):
		outData = list()

		maxNumbDataSeries = max( [x.shape[1] for x in plotData] )
		inpData = misc.getCycleListToMaxLength(propInInputFormat, maxNumbDataSeries)

		for mIdx, methodData in enumerate(plotData):
			for cIdx in range(1,methodData.shape[1]):
				outData.append( inpData[cIdx-1] )
		return outData


	def changeLegendEntriesToMethodOnly(self, outFig):
		currLegend = outFig.get_axes()[0].get_legend()
		usefulLines, usefulText = list(), list()

		#We grab the first instance of a new method
		for handle,textObj in it.zip_longest(currLegend.legendHandles, currLegend.texts):
			if textObj.get_text() not in usefulText:
				usefulText.append( textObj.get_text() )
				usefulLines.append( handle )

		plt.legend(usefulLines,usefulText)
			



