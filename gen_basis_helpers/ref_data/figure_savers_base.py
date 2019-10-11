
import os
import pathlib

from ..shared.data_plot_base import temporarilySetDataPlotterRegisteredAttrs
from ..shared import misc_utils as misc

class FigureSaver():
	"""Contains all info needed to save a figure. A different base object is expected per paper
	"""
	registeredKwargs = set()
	registeredKwargs.add("fmt")
#	registeredKwargs.add("fullPath")
	registeredKwargs.add("baseFolder")
	registeredKwargs.add("fileName")

	def saveFigure(self, figHandle, **kwargs):
		""" Saves the figure with options specified by instance attributes and kwargs passed
		
		Args:
			figHandle: pyplot handle to the figure (or anything with the same .savefig interface) 
			kwargs: see instance.registeredKwargs. These are the same names as object attrs, and **kwargs are taken in preference
			to object attributes
		
		Raises:
			Errors
		"""
		raise NotImplementedError("")


class StandardFigureSaver(FigureSaver):

	def __init__(self, **kwargs):
		""" Initialiser for StandardFigureSaver; DO NOT call this directly, use fromDefaultPlusKwargs
		
		Args:
			kwargs: See self.registeredKwargs
		"""
		for key in self.registeredKwargs:
			setattr(self,key,None)

		for key in kwargs:
			if key in self.registeredKwargs:
				setattr(self,key,kwargs[key])
			else:
				raise KeyError("{} is an invalid keyword.\n Available kwargs are {}".format(key , self.registeredKwargs))


	@classmethod
	def fromDefaultPlusKwargs(cls, **kwargs):
		#Step 1 = check we have at least the minimal kwargs
		reqKwargs = ["baseFolder"]
		for kwarg in reqKwargs:
			assert kwarg in kwargs, "{} is a required argument".format(kwarg)
		#Step 2 = set any defaults
		outKwargs = dict()
		outKwargs["fileName"] = "test_fig"
		outKwargs["fmt"] = "eps"

		#Step 3 = update defaults with user options, then create the object
		outKwargs.update(kwargs)
		outKwargs["baseFolder"] = os.path.abspath(outKwargs["baseFolder"])
		return cls(**outKwargs)


	def saveFigure(self, figHandle, **kwargs):
		with temporarilySetDataPlotterRegisteredAttrs(self,kwargs):
			outPath = os.path.abspath( os.path.join(self.baseFolder, self.fileName) )
			outFolder = os.path.split(outPath)[0]
			pathlib.Path(outFolder).mkdir(exist_ok=True,parents=True)
			if os.path.splitext(outPath)[1]=="":
				outPath += ".{}".format( self.fmt )
			figHandle.savefig( outPath, format=self.fmt, bbox_inches='tight' )







