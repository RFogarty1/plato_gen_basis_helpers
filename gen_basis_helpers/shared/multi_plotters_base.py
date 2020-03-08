
""" Classes used to create figures which contain multiple plots together """


class MultiPlotterBase():
	"""Class use to create figures containing multiple plots

	"""

	def create(self):
		""" Create the multi-plot figure and return the handle to it
		"""
		raise NotImplementedError("")


class MultiPlotGridBase():
	""" Class encapsulates how we arrange single-plots in a figure with multiple plots.
	"""

	def create(self):
		""" Create the basic figure, with empty axes
		
		Returns
			(figHandle, axisHandles)
			figHandle: The handle to the matplotlib figure created
			axisHandles: (iter of axis handles) Returning this allows the plotter factories to fill these with data
	 
		"""
		raise NotImplementedError("")

class MultiPlotterStandard(MultiPlotterBase):
	
	def __init__(self, plotterFactories, gridCreator):
		""" Initializer
		
		Args:
			plotterFactories: (iter of DataPlotterBase objects), each has a create method which takes an axis asd a kwarg, and creates the relevant plot
			gridCreator: (MultiPlotGridBase object) has a create method which returns (figHandle, axisHandles). These axisHandles are then populated by the plotterFactories objects
				 
		"""
		self.plotterFactories = plotterFactories
		self.gridCreator = gridCreator

	def create(self):
		outFig, allAxes = self.gridCreator.create()
		for plotFactory, axis in zip(self.plotterFactories, allAxes):
			plotFactory.createPlot( axHandle=axis )
		return outFig

