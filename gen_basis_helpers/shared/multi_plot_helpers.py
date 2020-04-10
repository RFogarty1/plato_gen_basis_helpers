

""" Module is meant for containing simple convenience functions for making simple plots with multiple sub-figures """

from . import multi_plot_grids as multiGrids
from . import multi_plotters_base as mPlotHelp

def createSimpleMultiPlot(plotFactories, nCols=None, figSize=None, multiPlotterKwargDict=None):
	""" Create a matplotlib figure containing multiple, independent plots
	
	Args:
		plotFactories: (iter of DataPlotterBase objects) Each represents a single plot. 
		nCols: (optional, int) Default is two (except if only one plotFactory is passed, then 1 column)
		figSize: (optional, tuple) Passed to RectangularPlotGrid where its used to determinr final figure size. (width,height)
		multiPlotterKwargDict: (optional, dict) Dictionary containig kwarg:value for passing keyword arguments to the MultiPlotterStandard object (see documentation for THAT object for available options)

	Returns
		 outFig: Handle to matplotlib figure object
 
	"""
	if nCols is None:
		nCols = 1 if len(plotFactories)==1 else 2
	if multiPlotterKwargDict is None:
		multiPlotterKwargDict = dict()


	nPlots = len(plotFactories)
	gridConfig = multiGrids.RectangularPlotGrid(nCols, nPlots, figSize=figSize)
	multiPlotFactory = mPlotHelp.MultiPlotterStandard(plotFactories, gridConfig, **multiPlotterKwargDict)
	return multiPlotFactory.create()




