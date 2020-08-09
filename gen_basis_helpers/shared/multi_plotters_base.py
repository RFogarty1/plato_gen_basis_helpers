
import string
import itertools as it
from . import misc_utils as misc

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
	
	def __init__(self, plotterFactories, gridCreator, annotateGraphs=False, annotateStrs=None, annotatePositionsAbs=None,
	             annotatePositionsRel=None, annotateFontSize=None):
		""" Initializer
		
		Args:
			plotterFactories: (iter of DataPlotterBase objects), each has a create method which takes an axis asd a kwarg, and creates the relevant plot
			gridCreator: (MultiPlotGridBase object) has a create method which returns (figHandle, axisHandles). These axisHandles are then populated by the plotterFactories objects

		Optional Args:
			annotateGraphs: (bool) If True annotations will be added to the graphs to label them (e.g. as a, b, c or similar)
			annotateStrs: (iter of strs) Each element is the annotation label for one graph. If not set, then a sensible default is used
			annotatePositionsAbs: (iter of 2-element iters OR a single 2-element iter) Determines the positions of annotateStrs on 
			                   each plot. These are passed as the xy element to matplotlibs annotate function
			annotatePositionsRel: (iter of 2 element iters OR a single 2-element iter) Same as above except that the positions are defined in terms relative to the x/y limits. e.g
			                      [0.5,0.5] would mean putting it directly in the middle of the plot area
			annotateFontSize: (iter of floats OR a single float) Each element is the font size for each annotation
				 
		"""
		#We create this dict only to check input vals are consistent
		kwargDict = {"annotateGraphs":annotateGraphs, "annotateStrs":annotateStrs, "annotatePositionsAbs":annotatePositionsAbs, "annotatePositionsRel":annotatePositionsRel}
		self._checkKwargDictConsistent(kwargDict)

		self.plotterFactories = plotterFactories
		self.gridCreator = gridCreator
		self.annotateGraphs = annotateGraphs
		self.annotateStrs =  annotateStrs if annotateStrs is not None else [x for x in string.ascii_lowercase]
		self.annotatePositionsAbs = annotatePositionsAbs 
		self.annotatePositionsRel = annotatePositionsRel 
		self.annotateFontSize = annotateFontSize

		if (self.annotatePositionsAbs is None) and (self.annotatePositionsRel is None):
			self.annotatePositionsRel = [0.05,0.9]


	def _checkKwargDictConsistent(self, kwargDict):
		annotatePosAbs = kwargDict.get("annotatePositionsAbs",None)
		annotatePosRel = kwargDict.get("annotatePositionsRel",None)

		if ( (annotatePosAbs is not None) and (annotatePosRel is not None) ):
			raise ValueError( ("Either annotatePositionsAbs OR annotatePositionsRel can be set; not both at once."
			                   "Values are annotatePositionsAbs={}, annotatePositionsRel={}".format(annotatePosAbs,annotatePosRel) ) )

	def create(self, **kwargs):
		""" Create the multi-plot figure and return the handle to it. 

		Kwargs are object attributes; if any are set then those values are used to create the plot but the state of this object
		will NOT change 
		"""

		self._checkKwargDictConsistent(kwargs)
		with misc.temporarilySetInstanceAttrs(self,kwargs):
			outObj = self._createFromSelf()
		return outObj

	def _createFromSelf(self):
		outFig, allAxes = self.gridCreator.create()
		for plotFactory, axis in zip(self.plotterFactories, allAxes):
			plotFactory.createPlot( axHandle=axis )

		if self.annotateGraphs:
			self._addAnnotations(allAxes)
		return outFig

	@property
	def annotatePositionsAbs(self):
		return self._annotatePositionsAbs

	@annotatePositionsAbs.setter
	def annotatePositionsAbs(self,val):
		if val is None:
			self._annotatePositionsAbs = None
			return None
		else:
			self._annotatePositionsAbs = val

		self._annotatePositionsRel = None

	@property
	def annotatePositionsRel(self):
		return self._annotatePositionsRel

	@annotatePositionsRel.setter
	def annotatePositionsRel(self,val):
		if val is None:
			self._annotatePositionsRel = None
			return None
		else:
			self._annotatePositionsRel = val
		self._annotatePositionsAbs = None


	def _getIterFromAnnotatePos(self, val):
		try:
			_ = val[0][0] #Will work if this is a list of 2-element iterables
		except TypeError:
			outVal = it.cycle([val])	#Each element needs to be a 2-element iter	
		else:
			outVal = it.cycle(val)
		return outVal

	def _getCycledIterableFromIterOrSingleVal(self, val):
		try:
			_ = val[0] #Will work if we have an iterable
		except TypeError:
			outVal = it.cycle([val])
		else:
			outVal = it.cycle(val)
		return outVal

	def _addAnnotations(self, allAxes):
		useAbsPos = True if self.annotatePositionsAbs is not None else False
		if useAbsPos:
			annotatePos = self.annotatePositionsAbs
		else:
			annotatePos = self.annotatePositionsRel

		annotatePos = self._getIterFromAnnotatePos(annotatePos)
		annotateStrs = it.cycle(self.annotateStrs)
		fontSize = [x for unused,x in zip(allAxes,self._getCycledIterableFromIterOrSingleVal(self.annotateFontSize)) ]

		for idx, (ax,inpStr,inpPos) in enumerate( zip(allAxes, annotateStrs, annotatePos) ):
			annotateKwargs = {"fontsize":fontSize[idx]}
			if useAbsPos:
				ax.annotate(inpStr,inpPos, **annotateKwargs)
			else:
				currXlim, currYlim = ax.get_xlim(), ax.get_ylim()
				xRange, yRange = currXlim[1]-currXlim[0], currYlim[1]-currYlim[0]
				absPos = [ currXlim[0] + (xRange*inpPos[0]), currYlim[0] + (yRange*inpPos[1]) ]
				ax.annotate(inpStr,absPos, **annotateKwargs)
