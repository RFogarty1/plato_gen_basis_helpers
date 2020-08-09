

import numpy as np
import matplotlib.pyplot as plt

from ..shared import data_plot_base as basePlotter

class DataPlotterDos(basePlotter.DataPlotterStandard):

	def __init__(self, **kwargs):
		self.registeredKwargs.add("lineStyles")
		super().__init__(**kwargs)

	@classmethod
	def fromDefaultPlusKwargs(cls, **kwargs):
		inpKwargs = dict()
		inpKwargs["plotFunct"] = "plot"
		inpKwargs["xlabel"] = "Energy (eV)"
		inpKwargs["ylabel"] = "Density of States"
		inpKwargs["showTitle"] = True
		inpKwargs["titleStr"] = "Default Title"
		inpKwargs.update(kwargs)
		return cls(**inpKwargs)


