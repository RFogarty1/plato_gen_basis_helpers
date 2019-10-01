
import itertools as it

import numpy as np
import matplotlib.pyplot as plt

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
			for idx,(handle,style) in enumerate( zip(dataLines,lineStyles) ):
				handle.set_linestyle(style)
				#Need to modify the style in the legend too (which holds a copy of the line)
				if self.legend:
					legLineHandles = outFig.get_axes()[0].get_legend().get_lines()
					legLineHandles[idx].set_linestyle(style)


				
		#Which means i need to possibly redo the legend
		if self.legend:
			plt.legend()

		return outFig

