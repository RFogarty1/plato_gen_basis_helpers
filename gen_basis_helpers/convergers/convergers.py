
import itertools as it
import numpy as np
import plato_pylib.utils.job_running_functs as jobRun

from . import data_plot_conv as dPlot
from . import tabulator_conv as tabCreate

class DuplicateWorkFoldersError(ValueError):
	pass

#TODO: Havent tested for duplication upon initialisation (currently no error will be thrown)
class PropConvJobRunComposite():
	""" Composite Object used to run jobs for converging calculations w.r.t a parameter

	Attributes (including @properties):
		shouldWeRunCalcs (bool): True=run jobs, False=Dont. Even if set to true, it wont overwrite lower level branches (though if False that WILL overwrite them) 
		runComms: list of run commands for bash (which are used to run the jobs)

	"""

	def __init__(self, inpObjs, shouldWeRunCalcs=True):
		""" 
		Args:
			inpObjs: iter of PropConvJobRunner objects
		    shouldWeRunCalcs: (optional) If true (default) then jobs will be run, else they wont be
		"""
		self._branchObjs = list(inpObjs)
		self._shouldWeRunCalcs = shouldWeRunCalcs

		self._ensureNoDuplicateWorkFolders()


	def _ensureNoDuplicateWorkFolders(self):
		if ( len(self.workFolders) != len(set(self.workFolders)) ):
			raise DuplicateWorkFoldersError("Each PropConvJobRunnner needs to have a unique workFolder")


	@property
	def workFolders(self):
		allWorkFolders = list()
		for x in self._branchObjs:
			allWorkFolders.extend(x.workFolders)
		outWorkFolders = [x for x in allWorkFolders if x is not None]
		return outWorkFolders

	@property
	def shouldWeRunCalcs(self):
		return self._shouldWeRunCalcs

	@shouldWeRunCalcs.setter
	def shouldWeRunCalcs(self, value):
		assert(value is bool)
		self._shouldWeRunCalcs = True

	def runJobs(self,nCores=None):
		if nCores is None:
			nCores = 1
		jobRun.executeRunCommsParralel(self.runComms, nCores)


	@property
	def runComms(self):
		outComms = list()
		for x in self._branchObjs:
			if x.shouldWeRunCalcs is True:
				outComms.extend( x.runComms )
		return outComms

	def createAnalyser(self):
		outAnalysers = list()
		for branch in self._branchObjs:
			outAnalysers.append( branch.createAnalyser() )
		return PropConvAnalyserComposite(outAnalysers)


class PropConvAnalyserComposite():
	"""Object used to interact with data from a set of convergence calculations

	Attributes:
		data (list): Each entry contains an array with [convVals:propVals] 

	"""

	def __init__(self, inpObjs):
		self._branchObjs = list(inpObjs)
		self._useDataRelativeToLargestConvValue = False

	@property
	def data(self):
		outData = list()
		for x in self._branchObjs:
			outData.extend(x.data)
		return outData


	@property
	def useDataRelativeToLargestConvValue(self):
		raise NotImplementedError("Only the setter is implemented")

	@useDataRelativeToLargestConvValue.setter
	def useDataRelativeToLargestConvValue(self,value):
		for x in self._branchObjs:
			x.useDataRelativeToLargestConvValue = value


	def plotData(self, **kwargs):
		outData = list()
		for x in self._branchObjs:
			outData.extend( x.plotData(**kwargs) )
		return outData



	def tabulateData(self,**kwargs):
		outData = list()
		for x in self._branchObjs:
			outData.extend( x.tabulateData(**kwargs) )
		return outData





class PropConvJobRunner():
	""" Class used to run jobs for property convergence calculations. Analysis can then be done using the object created by createAnalyser()

	Attributes:
		label (str): Label describing the set of convergence calculations
		shouldWeRunCalcs (bool): Determines whether jobs should actually be run or not
		runComms (iter): Returns iter of run commands for bash (can be empty list)
		createAnalyser: Returns an analyser object if calcs have been run previously
	"""

	@property
	def label(self):
		raise NotImplementedError("")
		return None

	@property
	def shouldWeRunCalcs(self):
		raise NotImplementedError("")
		return None

	@property
	def runComms(self):
		raise NotImplementedError("")
		return None

	@property
	def workFolders(self):
		raise NotImplementedError("")
		return None

	def createAnalyser(self):
		raise NotImplementedError("")
		return None

	def runJobs(self):
		jobRun.executeRunCommsParralel(self.runComms) #cant parralelise the single job



class PropConvAnalyser():
	"""Object used for analysis/display of results from property-convergence calculations

	Attributes:
		label (str): Description of `attr1`.
		data (list): First entry contains array of convergence values vs calculated property values 

	"""

	@property
	def label(self):
		raise NotImplementedError("")
		return None

	@property
	def data(self):
		raise NotImplementedError("")
		return None

	@property
	def useDataRelativeToLargestConvValue(self):
		return self._useDataRelativeToLargestConvValue

	@useDataRelativeToLargestConvValue.setter
	def useDataRelativeToLargestConvValue(self,value):
		self._useDataRelativeToLargestConvValue = value


	def plotData(self):
		""" Creates a plot for convergence data and returns a list with handle to the graph in first index
		
		Returns
			outFig: list, first entry has figure handle. List is used to make creating composite objects easier
		
		"""
		raise NotImplementedError("")
		return None

	def tabulateData(self):
		raise NotImplementedError("")
		return None







class VaryParamToOptDictMapper():
	"""Class to modify an option dict (for a plato file, though could be extended) based on the value of varParam.
	e.g. varParam may represent a grid spacing, the logic of modifying the optDict is then encoded into the class
	[Base class docstring]

	"""

	def modOptDictWithVariableParam(self, optDict, varParam):
		""" Modifies optDict (IN PLACE) based on the value of varParam [Base class docstring]
		
		Args:
			optDict(dict): Dictionary representing options for an input file
			varParam: Parameter (e.g. maybe grid spacing) which you want modified in optDict
	
		Returns
			Nothing, modifies optDict in place
		
		"""
		raise NotImplementedError("")
		return None



class PropConvJobRunnerStandard(PropConvJobRunner):

	def __init__(self, workFlows, varyParams, label, shouldWeRunCalcs=True):
		""" Description of function
		
		Args:
			workFlows: (list) of WorkFlow objects (see WorkFlowBase in fit_plato_ints code)
			varyParams: (list) of the value taken by variable parameters. These numbers arent used in the calculations, just in the creation of an analyser object.
			label: (str) Labels the prop-conv runner, used so user can differentiate between them 
			shouldWeRunCalcs: (optional Bool, default=True) determines whether the plato calculations should be run (often set to False to save time, using previously run calculations)
				
		Raises:
			DuplicateWorkFoldersError: If two workFlows share a single folder
		"""
		self.workFlows = list(workFlows)
		self.varyParams = varyParams
		self._label = label
		self._shouldWeRunCalcs = shouldWeRunCalcs
		self._ensureNoDuplicateWorkFolders()

	def _ensureNoDuplicateWorkFolders(self):
		if len(self.workFolders) != len(set(self.workFolders)):
			raise DuplicateWorkFoldersError("Each workFlow needs a unique workFolder")

	@property
	def label(self):
		return self._label

	@property
	def shouldWeRunCalcs(self):
		return self._shouldWeRunCalcs

	@shouldWeRunCalcs.setter
	def shouldWeRunCalcs(self, value):
		assert isinstance(value,bool), "shouldWeRunCalcs needs to be a boolean"
		self._shouldWeRunCalcs = value

	@property
	def runComms(self):
		allComms = list()
		for x in self.workFlows:
			currComms = x.preRunShellComms
			if currComms is not None:
				allComms.extend(currComms)
		return allComms

	@property
	def workFolders(self):
		outFolders = list()
		for x in self.workFlows:
			outFolders.append( x.workFolder )
		outFolders = [x for x in outFolders if x is not None]
		return outFolders


	def runJobs(self):
		if self.shouldWeRunCalcs:
			jobRun.executeRunCommsParralel(self.runComms) #cant parralelise the single job

	def createAnalyser(self, outputField=None):
		#Do all the parsing/analysis required to generate the output. Possibily cache results if it turns out to be slow
		for x in self.workFlows:
			x.run()

		if outputField is None:
			keyVals = vars(self.workFlows[0].output).keys() 
			if ( len(keyVals)==1 ):
				outputField = [k for k in keyVals][0]
			else:
				raise ValueError("outputField value ambigous (please set explicitly); could be any of {}".format(keyVals))

		xVals = self.varyParams
		yVals = [getattr(x.output, outputField) for x in self.workFlows]

		convVals = [(x,y) for x,y in it.zip_longest(xVals,yVals)]

		return PropConvAnalyserStandard(convVals, self.label)



class PropConvAnalyserStandard(PropConvAnalyser):

	def __init__(self, convData, label, dataPlotter=None, tabulator=None):
		self._data = [np.array(convData)]
		self._label = label
		self._useDataRelativeToLargestConvValue = False

		if tabulator is None:
			self.tabulator = tabCreate.ConvDataTabulator(fmt='html', titleStr=self.label)
		else:
			self.tabulator = tabulator

		if dataPlotter is None:
			self.dataPlotter = dPlot.DataPlotterConvergers.fromDefaultPlusKwargs()
		else:
			self.dataPlotter = dataPlotter


	@property
	def label(self):
		return self._label

	@property
	def useDataRelativeToLargestConvValue(self):
		return self._useDataRelativeToLargestConvValue

	@useDataRelativeToLargestConvValue.setter
	def useDataRelativeToLargestConvValue(self,value):
		self._useDataRelativeToLargestConvValue = value
		

	@property
	def data(self):
		if self.useDataRelativeToLargestConvValue:
			outData = list()
			for x in self._data:
				currData = np.array(x)
				maxIndex = np.argmax(x[:,0])
				refVal = currData[maxIndex,1]
				currData[:,1] -= refVal
				outData.append(currData)
			return outData
		else:
			return self._data

	def plotData(self, **kwargs):
		inpKwargs = {"titleStr":self.label}
		inpKwargs.update(kwargs) #Overwrite any inpKwargs with use options
		return [self.dataPlotter.createPlot(self.data, **inpKwargs)]

	def tabulateData(self, **kwargs):
		inpKwargs = {"titleStr":self.label}
		inpKwargs.update(kwargs)
		return [self.tabulator.createTable(self.data[0], **inpKwargs)] #0 is pretty dumb really, we shouldnt actually be storing this as a list since there can only ever be 1 entry (unless i want to combine analysers in a weird way)

