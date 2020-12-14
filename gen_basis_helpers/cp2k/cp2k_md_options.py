

class MolDynamicsOptsCP2KBase():
	""" Object representing molecular dynamics specific options (e.g. ensemble to use) for CP2K"""

	@property
	def optDict(self):
		raise NotImplementedError("")




class MolDynamicsOptsCP2KStandard():

	def __init__(self, timeStep=None, nSteps=None, ensemble=None, temperature=None):
		""" Initializer
		
		Args:
			timeStep: (float)
			nSteps: (int) Number of steps to run
			ensemble: (str) e.g. NVE, NPT. Passed directly to ensemble keyword in Motion/MD section
			temperature: (float) In Kelvin
				 
		Returns
			What Function Returns
	 
		Raises:
			Errors
		"""
		self.timeStep = timeStep
		self.nSteps = nSteps
		self.ensemble = ensemble
		self.temperature = temperature

		#Mapping from numbers to strings
		self.tempFmt = "{:.2f}"
		self.timeStepFmt = "{:.2f}"

	@property
	def optDict(self):
		outDict = {"mdEnsemble": self.ensemble, "mdSteps": self.nSteps,
		           "mdTimeStep": self.timeStepFmt.format(self.timeStep) if self.timeStep is not None else None,
		           "mdTemperature": self.tempFmt.format(self.temperature) if self.temperature is not None else None}

		outDict = {k:v for k,v in outDict.items() if v is not None}

		return outDict
