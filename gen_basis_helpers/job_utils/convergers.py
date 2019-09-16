
import plato_pylib.utils.job_running_functs as jobRun

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

	@property
	def shouldWeRunCalcs(self):
		return self._shouldWeRunCalcs

	@shouldWeRunCalcs.setter
	def shouldWeRunCalcs(self, value):
		assert(value is bool)
		self._shouldWeRunCalcs = True

	def runJobs(nCores=None):
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
		pass


class PropConvAnalyserComposite():
	"""Object used to interact with data from a set of convergence calculations

	Attributes:
		data (list): Each entry contains an array with [convVals:propVals] 

	"""


	def __init__(self, inpObjs):
		self._branchObjs = list(inpObjs)

	@property
	def data(self):
		outData = list()
		for x in self._branchObjs:
			outData.extend(x.data)
		return outData

	def plotData(self):
		pass

	def tabulateData(self):
		pass





class PropConvJobRunner():
	"""The summary line for a class docstring should fit on one line.

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

	def createAnalyser(self):
		raise NotImplementedError("")
		return None

	def runJobs(self):
		jobRun.executeRunCommsParralel(self.runComms) #cant parralelise the single job



#class PropConvAnalyser():

