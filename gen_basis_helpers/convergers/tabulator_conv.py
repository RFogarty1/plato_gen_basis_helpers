
import tabulate
from . import data_plot_conv

class TableMakerBase():

	def createTable(self, tableData):
		""" Takes data in tableData and returns a table
		
		Args:
			tableData: Data to tabulate, the exact format can be decided by implementation classes [this is the base class docstring]
				
		Returns
			handle to the table (recommend the same struct tabulate library returns)
		
		"""
		raise NotImplementedError


#TODO: Most the logic of this class is identical to DataPlotterConvergers, whcih suggests i need a refactor
@data_plot_conv.addSetOfKwargDescriptorsToClass
class ConvDataTabulator():

	registeredKwargs = set()
	registeredKwargs.add("fmt")
	registeredKwargs.add("titleStr")
	registeredKwargs.add("tableHeadings")
	registeredKwargs.add("numbDecimalAll")
	def __init__(self, **kwargs):
		for key in self.registeredKwargs:
			if key in kwargs:
				setattr(self,key,kwargs[key])
			else:
				setattr(self,key,None)



	def _updateAttrsFromKwargs(self, **kwargs):
		for key in self.registeredKwargs:
			if key in kwargs.keys():
				setattr(self,key,kwargs[key])

	def _getDictForAllRegisteredAttrKwargs(self):
		outDict = {k:getattr(self,k) for k in self.registeredKwargs}
		return outDict


	def createTable(self, tableData, **kwargs):
		""" Creates table (tabulate format) from tableData
		
		Args:
			tableData: nx2 numpy array 
			kwargs: keyword-arguments in same format as when assigning to the object attributres (keys in self.registeredKwargs). These are assigned to the object for solely this function call
				
		Returns
			table: Format returned by tabulate
		
		Raises:
			Errors
		"""

		#Backup the state from before the function call. TODO: Need to use a context manager to do this more safely
		startKwargs = self._getDictForAllRegisteredAttrKwargs()
		self._updateAttrsFromKwargs(**kwargs)

	
		if self.numbDecimalAll is None:
			strFormat = "{}"
		else:
			strFormat = "{:." + str(self.numbDecimalAll) + "f}"


		#Get the title/heading in the table
		tabList = list() #Each entry is 1 row
		if self.titleStr is not None:
			tabList.append( [self.titleStr] )

		if self.tableHeadings is not None:
			tabList.append( self.tableHeadings )

		#Get the actual data into the table
		for rIdx in range(tableData.shape[0]):
			currRow = [strFormat.format(x) for x in tableData[rIdx,:]]
			tabList.append(currRow)

		if self.fmt is None:
			outTable = tabulate.tabulate( tabList )
		else:
			outTable = tabulate.tabulate( tabList, tablefmt=self.fmt )


		#Reset state to before the function call
		self._updateAttrsFromKwargs(**startKwargs)

		return outTable

