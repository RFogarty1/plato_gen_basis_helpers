

import matplotlib.pyplot as plt
from . import multi_plotters_base as baseObjs


#See https://matplotlib.org/tutorials/intermediate/gridspec.html for description of underling mpl frunctions
class RectangularPlotGrid(baseObjs.MultiPlotGridBase):
	""" Class representing a rectangular grid of plots
	"""
	def __init__(self, nCols, nPlots):
		""" Initializer
		
		Args:
			nCols: (int) Number of columns per row
			nPlots: (int) Total number of plots required
		"""
		self.nCols = nCols
		self.nPlots = nPlots


	def create(self):
		outFig = plt.figure(constrained_layout=True)
		gridSpec = outFig.add_gridspec(*self.dims)
		outAxes = self._addSubPlotsToFigure(outFig,gridSpec)
		return (outFig,outAxes)

	#TODO: Later need to sort this out such that some plots can span multi-rows/columns
	def _addSubPlotsToFigure(self, figHandle, gridSpec):
		allAxes = list()
		nRows, nCols = self.dims
		for rIdx in range(nRows):
			for cIdx in range(nCols):
				currAxis = figHandle.add_subplot(gridSpec[rIdx,cIdx])
				allAxes.append(currAxis)
		return allAxes

	@property
	def dims(self):
		""" (nRows,nColums)
		"""
		nRows = 1 #1 row is the default/minimum
		outNPlots = self.nCols
		while outNPlots < self.nPlots:
			outNPlots += self.nCols
			nRows += 1
		return (nRows,self.nCols)
