

class MolDynamicsOptsCP2KBase():
	""" Object representing molecular dynamics specific options (e.g. ensemble to use) for CP2K"""

	@property
	def optDict(self):
		raise NotImplementedError("")




class MolDynamicsOptsCP2KStandard():

	def __init__(self, timeStep=None, nSteps=None, ensemble=None, temperature=None, thermostatType=None):
		""" Initializer
		
		Args:
			timeStep: (float)
			nSteps: (int) Number of steps to run
			ensemble: (str) e.g. NVE, NPT. Passed directly to ensemble keyword in Motion/MD section
			temperature: (float) In Kelvin
			thermostatType: (Str) e.g. Nose
 
		Returns
			What Function Returns
	 
		Raises:
			Errors
		"""
		self.timeStep = timeStep
		self.nSteps = nSteps
		self.ensemble = ensemble
		self.temperature = temperature
		self.thermostatType = thermostatType

		#Mapping from numbers to strings
		self.tempFmt = "{:.2f}"
		self.timeStepFmt = "{:.2f}"

	@property
	def optDict(self):
		outDict = {"mdEnsemble": self.ensemble, "mdSteps": self.nSteps,
		           "mdTimeStep": self.timeStepFmt.format(self.timeStep) if self.timeStep is not None else None,
		           "mdTemperature": self.tempFmt.format(self.temperature) if self.temperature is not None else None,
	               "mdThermostatType": self.thermostatType if self.thermostatType is not None else None}

		outDict = {k:v for k,v in outDict.items() if v is not None}

		return outDict


class MetaDynamicsOptsCP2KStandard():

	def __init__(self, metaVars, ntHills=None, doHills=True, hillHeight=None, printColvarCommonIter=3, heightNumbFmt="{:.4f}"):
		""" Initializer
		
		Args:
			metaVars: (iter of MetaVarStandard objects)
			ntHills: (int) create a hill at most every N steps
			doHills: (Bool) True means create hills during the simulation; basically should always be true (maybe false when reading hills in)
			hillHeight: (float) Hill height in Ha (default CP2K units for this)
			heightNumbFmt: (str) format string for converting hillHeight into a string 

		"""
		self.metaVars = metaVars
		self.ntHills = ntHills
		self.doHills = doHills
		self.hillHeight = hillHeight
		self.printColvarCommonIter = printColvarCommonIter
		self.heightNumbFmt = heightNumbFmt

	@property
	def optDict(self):
		outDict = {"metaVars":self.metaVars, "metaDyn_doHills":self.doHills, 
		           "metaDyn_hillHeight": self.heightNumbFmt.format(self.hillHeight) if self.hillHeight is not None else None,
		           "metaDyn_printColvarCommonIterLevels": self.printColvarCommonIter,
		           "metaDyn_ntHills":self.ntHills}

		outDict = {k:v for k,v in outDict.items() if v is not None}

		return outDict

