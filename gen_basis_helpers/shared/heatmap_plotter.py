
import copy
import numpy as np

import matplotlib.pyplot as plt

from . import data_plot_base as dPlotHelp


class HeatMapPlotterSimple(dPlotHelp.DataPlotterBase):

	registeredKwargs = set(dPlotHelp.DataPlotterBase.registeredKwargs)
	registeredKwargs.add("addColorbar")


	#Overwriting docstring
	def createPlot(self, plotData=None, **kwargs):
		""" Takes data in plotData argument and creates a plot 
		
		Args:
			plotData: [ binEdgesArray, binDataArray ] Awkward format because usually data plotters take an iter for each data series
			binEdgesArray:  ( nxmx2x2 numpy array) Each [n][m] entry contains the bin edges [[minX,maxX],[minY,maxY]]. Higher dimension cases would contain more len-2 arrays; but we're limited to mapping 2-dim data here. Note n and m are the number of bins in each dimension
			binDataArray: (nxm array) Each entry contains the z-value (e.g. counts) for the bins defined in binEdgesArray[n][m]

		Returns
			Handle to the overall figure
	
		NOTES:
			a) the output from the function plt.pcolormesh is stored in self.plotOutput; this could be useful if wanting to add a colorbar later. However, this needs to be deleted if you want to copy the object later (e.g. after plotting with it once)
 
		"""
		return super().createPlot(plotData=plotData, **kwargs)


#			toPlot = self._getToPlotData(plotData)
#			outFig = self._getOutfigHandleAndSetCurrentAxis()

	def _getToPlotData(self, plotData=None):
		if plotData is not None:
			outX, outY, outZ = mapBinEdgesAndDataToXYZStandard(plotData[0], plotData[1])
			toPlot = [ [outX,outY,outZ] ] #Need to loop over this; hence has to be an array
		else:
			if self.data is not None:
				outX, outY, outZ = mapBinEdgesAndDataToXYZStandard(self.data[0], self.data[1])
				toPlot = [ [outX,outY,outZ] ] #Need to loop over this; hence has to be an array
			else:
				toPlot = list()

		return toPlot


	def _plotSingleDataSeries(self, inpData, label, idx):
		self.plotOutput = plt.pcolormesh(*inpData)


	def _modifyAxisPlot(self):
		if self.addColorbar:
			currAx = plt.gca() #We need both the axis and the figure to gaurantee adding to the right place
			currFig = currAx.get_figure()
			currFig.colorbar(self.plotOutput,ax=currAx)


def mapBinEdgesAndDataToXYZStandard(binEdgesArray, binDataArray):
	""" Standard function for mapping bin edges and data to the format matplotlib needs for making a "pcolormesh" plot (basically a heat map)
	
	Args:
		binEdgesArray: ( nxmx2x2 numpy array) Each [n][m] entry contains the bin edges [[minX,maxX],[minY,maxY]]. Higher dimension cases would contain more len-2 arrays; but we're limited to mapping 2-dim data here. Note n and m are the number of bins in each dimension
		binDataArray: (nxm array) Each entry contains the z-value (e.g. counts) for the bins defined in binEdgesArray[n][m]
 
	Returns
		X,Y,Z: Returned as a len-3 list. These should be passed sequentially into matplotlib.pyplot.pcolormesh to get a heatmap. X and Y are both (n+1)x(m+1) arrays and contain edges of the relevant rectangles we use in our plot (indices refer to nth rectangle across/up). Z contains all the values in an nxm matrix. 
 
	Raises:
		 ValueError: If the shape of binEdges/binData arrays are incorrect
	"""
	#TODO: Check that this is actually a 2-d array (need +1 since number of edges is +1 more than number of bins)
	dimX, dimY = len(binEdgesArray)+1, len(binEdgesArray[0])+1

	#We get a ValueError later anyway (at time of writing) but better to raise it here
	if len(binEdgesArray.shape)-2 != 2:
		nDims = len(binEdgesArray.shape)-2
		raise ValueError("Input array seems to be {} dimensions; we need 2 dimensions".format(nDims))

	outX = np.zeros( (dimX,dimY) )
	outY = np.zeros( (dimX,dimY) )
	outZ = copy.deepcopy(binDataArray)

	for idxX in range(dimX):
		for idxY in range(dimY):
			#Edge case where we take the upper edge of the previous bin (instead of lower edge of current)
			if (idxX==dimX-1) and (idxY==dimY-1):
				outX[idxX][idxY] = binEdgesArray[idxX-1][idxY-1][0][1] #Is this correct????
				outY[idxX][idxY] = binEdgesArray[idxX-1][idxY-1][1][1] 
			elif idxX==dimX-1:
				outX[idxX][idxY] = binEdgesArray[idxX-1][idxY][0][1]
				outY[idxX][idxY] = binEdgesArray[idxX-1][idxY][1][0]
			elif idxY==dimY-1:
				outX[idxX][idxY] = binEdgesArray[idxX][idxY-1][0][0]
				outY[idxX][idxY] = binEdgesArray[idxX][idxY-1][1][1]			
			else:
				outX[idxX][idxY] = binEdgesArray[idxX][idxY][0][0]
				outY[idxX][idxY] = binEdgesArray[idxX][idxY][1][0]

	return outX, outY, outZ



