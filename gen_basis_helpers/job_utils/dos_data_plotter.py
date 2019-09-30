

import itertools as it

import numpy as np
from ..shared import data_plot_base as basePlotter
from ..shared import misc_utils as misc

class DataPlotterDos(basePlotter.DataPlotterBase):

	def __init__(self, **kwargs):
		self.registeredKwargs.add("lineStyles")
		super().__init__(**kwargs)

	@classmethod
	def fromDefaultPlusKwargs(cls, **kwargs):
		inpKwargs = dict()
		inpKwargs["plotFunct"] = "plot"
		inpKwargs["xlabel"] = "Energy / eV"
		inpKwargs["ylabel"] = "Density of States"
		inpKwargs["showTitle"] = True
		inpKwargs["titleStr"] = "Default Title"
		inpKwargs.update(kwargs)
		return cls(**inpKwargs)


	def createPlot(self, plotData, **kwargs):
		outFig = super().createPlot(plotData, **kwargs)

		#Modify line styles if needed
		with misc.fragile(basePlotter.temporarilySetDataPlotterRegisteredAttrs(self,kwargs)):
			if self.lineStyles is None:
				raise misc.fragile.Break
			dataLines = outFig.get_axes()[0].get_lines()
			lineStyles = it.cycle(self.lineStyles)
			for handle,style in zip(dataLines,lineStyles):
				handle.set_linestyle(style)
				


		return outFig

