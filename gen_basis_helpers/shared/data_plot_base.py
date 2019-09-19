



class DataPlotterBase():

	def createPlot(self, plotData):
		""" Takes data in plotData argument and creates a plot 
		
		Args:
			plotData: Data to plot, the exact format can be decided by implementation classes [this is the base class docstring]
				
		Returns
			Handle to the overall figure

		
		"""
		raise NotImplementedError

