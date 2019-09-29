
import numpy as np
import matplotlib.pyplot as plt
from ..shared import data_plot_base as basePlotter




class DataPlotterConvergers(basePlotter.DataPlotterBase):

	def __init__(self, **kwargs):
		for key in self.registeredKwargs:
			if key in kwargs:
				setattr(self,key,kwargs[key])
			else:
				setattr(self,key,None)


	@classmethod
	def fromDefaultPlusKwargs(cls, **kwargs):
		inpKwargs = dict()
		inpKwargs["xlabel"] = "Convergence parameter"
		inpKwargs["ylabel"] = "Property Value / units"
		inpKwargs["showTitle"] = True
		inpKwargs["titleStr"] = "Convergence Plot"
		inpKwargs["plotFunct"] = "scatter"
		inpKwargs.update(kwargs) #user-defined kwargs overwrite defualt ones where relevant
		return cls(**inpKwargs)




