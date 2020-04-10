
""" Purpose of this is to allow two independent sets of axes/data to be plotted on a single graph """


class MultiAxisPlotterStandard():
	""" Class for combining two individual plots onto a single graph with 2 indendepent x and y axes

	"""
	def __init__(self, plotFactories, axHandle=None):
		""" Initialiser
		
		Args:
			plotFactories: (len 2 iter of DataPlotterBase objects)
		"""
		self.plotFactories = plotFactories

	def create(self):
		""" Create the single plot with two independent sets of axes
		
		Returns
			outPlot: matplotlib figure if no axis was provided. If axHandle was provided then return value will be None, but the plot will have been added to the axis
	 
		"""
		nFactories = len(self.plotFactories)
		if nFactories!=2:
			raise ValueError("Multi axis plotter can only handle 2 independent axes; but {} given".format(nFactories))

		axHandle, outFig = self._getAxisHandleAndOutFigHandle()
		self.plotFactories[0].createPlot(axHandle=axHandle)
		secondAxis = self._getTheSecondIndependentAxis(axHandle)
		self.plotFactories[1].createPlot(axHandle=secondAxis)

		return outFig

	def _getAxisHandleAndOutFigHandle(self):
		if self.plotFactories[0].axHandle is None:
			outFig = self.plotFactories[0].createPlot()
			outAx = outFig.get_axes()[0]
		else:
			outFig = None
			outAx = self.plotFactories[0].axHandle	
		return outAx, outFig

	def _getTheSecondIndependentAxis(self, inpAxis):
		return inpAxis.twinx().twiny()




