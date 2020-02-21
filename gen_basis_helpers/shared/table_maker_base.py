
import itertools as it
import tabulate

from .data_plot_base import temporarilySetDataPlotterRegisteredAttrs


class TableMakerBase():
	registeredKwargs = set()
	registeredKwargs.add("fmt")

	def __init__(self, **kwargs):
		for key in self.registeredKwargs:
			setattr(self,key,None)

		for key in kwargs:
			if key in self.registeredKwargs:
				setattr(self,key,kwargs[key])
			else:
				raise KeyError("{} is an invalid keyword.\n Available kwargs are {}".format(key , self.registeredKwargs))

	def createTable(self, tableData, **kwargs):
		""" Takes data in tableData and creates a table that can be passed to Tabulate (which will create it in whatever format you want)
		
		Args:
			tableData: list of input data. Each list entry is one data series (made from two columns). e.g. [ (np.array(xDataA,yDataA)), (np.array(xDataB,yDataB)) ]
			kwargs: keyword-arguments in same format as when assigning to the object attributres (keys in self.registeredKwargs). These are assigned to the object for solely this function call
	
		Returns
			outTable: Some kind of representation of a table, details depend on the fmt argument (which generally has the same options as the tabulate library)
		
		"""
		raise NotImplementedError("")


class TableMakerStandard(TableMakerBase):
	""" Class containing a createTable(self, tableData, **kwargs) used to create a specific type of table from data of a specific format

	Attributes:
		registeredKwargs: This contains a set of all keyword arguments associated with the class. These can be passed to the constructor
		                  or createTable. If passed to createTable they will only be set temporarily (the class state will be the same before/after the call).
		                  NOTE: Some of these are defined in the class, but some are only added upon initialisation. Hence you should query an INSTANCE rather 
		                  than the class

	Extra notes:
		mapInputToRowsFunct: This is possibly the most important of the registeredKwargs, it turns the input data (by default a list of x vs y) into a row format
		                     that tabulate can handle. Setting to "blankMapFunct" (defined in this module) means your passing data in the same format that tabulate accepts (i.e. its passed directly to tabulate)
	"""

	def __init__(self, **kwargs):
		self.registeredKwargs.add("headers")
		self.registeredKwargs.add("mapInputToRowsFunct")
		super().__init__(**kwargs)
		if self.mapInputToRowsFunct is None:
			self.mapInputToRowsFunct = InputToTabulateMapFunctionStandard()

	def createTable(self, tableData, **kwargs):
		with temporarilySetDataPlotterRegisteredAttrs(self, kwargs):
			outTable = self._createTable(tableData)
		return outTable


	def _createTable(self, tableData):
		reformedData = self.mapInputToRowsFunct(tableData)
		tabulateKwargDict = {"tablefmt":self.fmt,
		                     "headers":self.headers}
		tabulateKwargDict = {k:v for k,v in tabulateKwargDict.items() if v is not None}
		return tabulate.tabulate(reformedData, **tabulateKwargDict)



class InputToTabulateMapFunctionStandard():
	""" Callable class that determines the standard ways to turn input data to a way that works for tabulate. """

	def __init__(self):
		""" Initializer for callable-class which handles the mapping of input data to tabulate format. Empty for now (may add hooks later). 
			
			Note: The callable function takes data of the form [ [(x1,y1),(x2,y2)], [(x1,y3),(x2,y4)] ]. I.e. a list of iters which share the x-values	
		"""
		pass

	def _translateData(self, tableData):

		self._assertListLengthsEqual( tableData )
		self._assertAllXValsAreTheSame( tableData )

		allXVals = [x[0] for x in tableData[0]]

		#Could check the xvals here
		outData = list()
		for rowIdx,xVal in enumerate(allXVals):
			currRow = [xVal]
			for col in tableData:
				currRow.append( col[rowIdx][1] ) #the 1 is simply the y-index
			outData.append(currRow)

		return outData

	def _assertListLengthsEqual(self, tableData):
		listLens = [len(x) for x in tableData]
		assert all([x==listLens[0] for x in listLens]), "All input lists must be of the same length"

	def _assertAllXValsAreTheSame(self, tableData):
		allXVals = list()
		for col in tableData:
			xVals = [x[0] for x in col]
			allXVals.append(xVals)

		assert self._areAllListsTheSameWithinTolerance( allXVals ), "X-values must all be the same when creating a table"


	def _areAllListsTheSameWithinTolerance( self, allLists, tolerance=1e-8 ):
		refVals = allLists[0] #Use the first list as a reference
		for currList in allLists[1:]:
			diffList = [abs(x-ref) for x,ref in it.zip_longest(currList,refVals) ]
			if not all([x<tolerance for x in diffList]):
				return False 
		return True


	def __call__(self, tableData):
		return self._translateData(tableData)


def blankMapFunct(inpData):
	""" Simple function that returns the input argument as its output; "lambda x: x" would do the same """
	return inpData


