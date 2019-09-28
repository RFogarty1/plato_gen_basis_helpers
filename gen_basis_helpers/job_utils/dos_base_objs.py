




class DosAnalyserBase():
	"""Class for plotting/analysing DoS data

	Attributes:
		label (iter): List of DosLabel objects, used to identify/sort individual calculation
		data (nx2 iter or np array): density of states data (first col=energy, second=density of states)
		refData(nx2 iter or np array or None): Reference data (first col=energy, second=density of states)
		dataPlotter(DataPlotterBase object): Object used to handle how the data is plotted 

	"""


	def getObjectsWithComponents(self, components, caseSensitive=True):
		""" Returns analyser objects which match the components
		
		Args:
			components (str list): Strings which will be matched against DosLabel components. Any analyser with 
		all of these components present will be returned 
			caseSensitive(bool): Whether we need to match the case of each component or not 
		Returns
			objList (iter): List of Analyser objects. Function to combine into a composite is recommended, but not
		part of the required interface
	 
		"""
		raise NotImplementedError("")


	def plotData(self ,**kwargs):
		""" Create a plot of the density of states data
		
		Args:
			kwargs: Totally implementation dependent[this is the base class docstring]. No Kwargs should still return a reasonable plot though
				 
		Returns
			figHandle (iter): Handles to created figures (composite objects may create >1, hence we use an iter) [base-class docstring]
	 
		"""
		raise NotImplementedError("")


	def attachRefData(dosData, components, ignoreComponents=False, errorIfNoMatches=True, caseSenstiveComponents=True):
		""" Add reference data to relevant (or all) objects
		
		Args:
			dosData(nx2 iter or np array): First column is energy values, second is density of states value
			components(str list): Identifiers for the type of calculation, they are matched against components in labels 
			errorIfNoMatches(bool): True means will throw an error if input components dont match components of any objects
			caseSenstiveComponents(bool): Whether we need to match the case of each ignore Components or not 
 
		Raises:
			ValueError: If components dont match any object and errorIfNoMatches set to True
		"""
		raise NotImplementedError("")


class DosRunnerBase():
	"""Class for running DoS calculations.

	Attributes:
		dosComms: Bash commands for creating density of states data
		singlePointEnergyComms: Bash commands for running the single-point energy calculations required

	"""

	@property
	def dosComms(self):
		""" iter, each entry is a string for the bash command to create density of states data
		
		"""
		raise NotImplementedError()


	@property
	def singlePointEnergyComms(self):
		""" iter, each entry is a string for the bash command to run the total energy calculation (which is needed to generate Dos Data
		
		"""
		raise NotImplementedError()

	def runSinglePointEnergyCalcs(self, nCores=1):
		""" Runs the single-point energy calculations required to get density of states data
		
		Args:
			nCores: (int) Number of cores to run calculations over (default=1)
				
		Returns
			Nothing
		
		"""
		raise NotImplementedError()


	def runDosGeneratingCalcs(self, nCores=1):
		""" Calculates (and saves somewhere accesible) the density of states data 
		
		Args:
			nCores: (int, optional) Number of cores to run calculations over (default=1)
				
		Returns
			Nothing
		
		"""
		raise NotImplementedError()

	def createAnalyser(self):
		""" Creates a DosAnalyser (possibly composite) object - see DosAnalyserBase or similar
		
		Returns
			DosAnalyser object
		
		"""
		raise NotImplementedError()

