

from ..shared import misc_utils as misc


class BaseEosPropAnalyser():

	@property
	def label(self):
		""" Return a list of objects implementing BaseLabel interface 
		
		"""
		raise NotImplementedError("")

	def getPlotData(self):
		""" Returns List of data. Exact format is up to the implementing class (this is base class doc-string)
		
		"""
		raise NotImplementedError("")

	def getObjectsWithComponents(self, components, caseSensitive=True):
		""" Returns analyser objects which match the components
		
		Args:
			components (str list): Strings which will be matched against self.label.components (or inherited cls). 
			caseSensitive(bool): Whether we need to match the case of each component or not 
		Returns
			objList (iter): List of objects. Empty if object components dont match, else length 1 with reference to self. The iter is required
			                since the same interface needs to be present for composite/non-composite objects
	 
		"""
		raise NotImplementedError("")


class StandardEosPropAnalyserComposite(BaseEosPropAnalyser):

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=True)
	def __init__(self, inpObjs):
		self.objs = list(inpObjs)
		self._ensureLabelsAllDifferent()

	def _ensureLabelsAllDifferent(self):
		if len(list(self.label)) != len(set(self.label)):
			raise ValueError("Duplicate labels detected")

	@property
	def label(self):
		outList = list()
		for x in self.objs:
			outList.extend(x.label)

		return outList

	def getPlotData(self):
		outList = list()
		for x in self.objs:
			outList.extend(x.getPlotData())
		return outList

