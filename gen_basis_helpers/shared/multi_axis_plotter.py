

""" Purpose of this is to allow two independent sets of axes/data to be plotted on a single graph """


import matplotlib.pyplot as plt

class MultiAxisPlotterStandard():
	""" Class for combining two individual plots onto a single graph with 2 indendepent x and y axes

	"""
	def __init__(self, plotFactories)
		""" Initialiser
		
		Args:
			plotFactories: (len 2 iter of DataPlotterBase objects)
		"""
		self.plotFactories = plotFactories

	#This is pretty much purely so this works with multi-plot code
	def createPlot(self, axHandle=None):
		startAxHandle = self.plotFactories[0].axHandle
		if axHandle is not None:
			self.plotFactories[0].axHandle = axHandle
		output = self.create()
		self.plotFactories[0].axHandle = startAxHandle
		return output

	def create(self):
		""" Create the single plot with two independent sets of axes
		
		Returns
			outPlot: matplotlib figure if no axis was provided. If axHandle was provided then return value will be None, but the plot will have been added to the axis
	 
		"""
		nFactories = len(self.plotFactories)
		if nFactories!=2:
			raise ValueError("Multi axis plotter can only handle 2 independent axes; but {} given".format(nFactories))

		axHandle, outFig = self._getAxisHandleAndOutFigHandle()
		if outFig is None:
			self.plotFactories[0].createPlot(axHandle=axHandle)
		secondAxis = self._getSecondXIndependentAxis(axHandle)
		thirdAxis = self._getSecondYIndependentAxis(secondAxis)
		self.plotFactories[1].createPlot(axHandle=thirdAxis)
		self._ensureYLabelStillPresentOnSecondAxis(secondAxis)
		return outFig

	def _getAxisHandleAndOutFigHandle(self):
		if self.plotFactories[0].axHandle is None:
			outFig = self.plotFactories[0].createPlot()
			outAx = outFig.get_axes()[0]
		else:
			outFig = None
			outAx = self.plotFactories[0].axHandle	
		return outAx, outFig

	def _getSecondXIndependentAxis(self, inpAxis):
		return inpAxis.twinx()

	def _getSecondYIndependentAxis(self,inpAxis):
		return inpAxis.twiny()

	def _ensureYLabelStillPresentOnSecondAxis(self,inpAxis):
		if self.plotFactories[1].ylabel is None:
			return None

		if inpAxis.get_ylabel() == "":
			inpAxis.set_ylabel(self.plotFactories[1].ylabel)




