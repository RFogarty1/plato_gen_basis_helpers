

class MolDynamicsOptsCP2KBase():
	""" Object representing molecular dynamics specific options (e.g. ensemble to use) for CP2K"""

	@property
	def optDict(self):
		raise NotImplementedError("")


class ThermostatOptsBase():
	""" Object representing options for the thermostat. See addToPyCp2kObj for interface """

	def addToPyCp2kObj(self, pyCp2kObj):
		""" Adds the thermostat info to the pycp2k obj
		
		Args:
			pyCp2kObj: Backend object used to generate input files
				 
		Returns
			Nothing but adds the thermostat info to the pycp2k obj
	 
		Raises:
			Errors
		"""
		raise NotImplementedError("")



class MolDynamicsOptsCP2KStandard():

	def __init__(self, timeStep=None, nSteps=None, ensemble=None, temperature=None,
	             thermostatType=None, printKindTemp=None):
		""" Initializer
		
		Args:
			timeStep: (float)
			nSteps: (int) Number of steps to run
			ensemble: (str) e.g. NVE, NPT. Passed directly to ensemble keyword in Motion/MD section
			temperature: (float) In Kelvin
			thermostatType: (Str) e.g. Nose
			printKindTemp: (Bool) If True print temperatures of atomic kinds (dumps to a file)
 
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
		self.printKindTemp = printKindTemp

		#Mapping from numbers to strings
		self.tempFmt = "{:.2f}"
		self.timeStepFmt = "{:.2f}"

	@property
	def optDict(self):
		outDict = {"mdEnsemble": self.ensemble, "mdSteps": self.nSteps,
		           "mdTimeStep": self.timeStepFmt.format(self.timeStep) if self.timeStep is not None else None,
		           "mdTemperature": self.tempFmt.format(self.temperature) if self.temperature is not None else None,
	               "mdThermostatType": self.thermostatType if self.thermostatType is not None else None,
		           "mdPrintKindTemps": self.printKindTemp if self.printKindTemp is not None else None}

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


class NoseThermostatOpts(ThermostatOptsBase):

	def __init__(self, length=None, mts=None, yoshida=None, timeCon=None, timeConFmt="{:.1f}"):
		""" Initializer
		
		Args:
		Most map directly to a cp2k keyword. Leaving as None means dont specifiy in the input file 
		(and therefore use the CP2K defaults)
			timeConFmt: (str) Format string for the 
				 
		"""
		self.length = length
		self.mts = mts
		self.yoshida = yoshida
		self.timeCon = timeCon
		self.timeConFmt = timeConFmt

	def addToPyCp2kObj(self, pyCp2kObj):
		thermoSect = pyCp2kObj.CP2K_INPUT.MOTION.MD.THERMOSTAT
		thermoSect.Type = "nose"
		thermoSect.NOSE.Length = self.length
		thermoSect.NOSE.Mts = self.mts
		thermoSect.NOSE.Timecon = self.timeConFmt.format(self.timeCon)


class LangevinThermostatOpts(ThermostatOptsBase):

	def __init__(self, gamma=None, noisyGamma=None):
		""" Initializer
		
		Args:
			Gamma: Langevin friction parameter
			noisyGamma: Intrinsic friction parameter; this is used to correct for energy dissipation effects in second generation CP-like dynamics runs
		
		"""
		self.gamma = gamma
		self.noisyGamma = noisyGamma

	def addToPyCp2kObj(self, pyCp2kObj):
		langevinSection = pyCp2kObj.CP2K_INPUT.MOTION.MD.LANGEVIN
		langevinSection.Gamma = self.gamma
		langevinSection.Noisygamma = self.noisyGamma



class ThermalRegion():
	""" Object representing a thermal region in MD. Used to run various types of Langevin (e.g. different gamma values) or mixed Langevin/NVE simulations """

	def __init__(self, atomList=None, doLangevin=None, noisyGamma=None, temperature=None, baseZeroAtomList=True):
		""" Initializer
		
		Args:
			atomList: (iter of ints) List of atom integers to be put in this region. Whether we use base-zero or base-one numbering is determined by baseZeroAtomList
			doLangevin: (Bool)
			noisyGamma: (float, Optional) Value of noisy gamma to overwrite the global value
			temperature: (float, Optional) Temperature for this region; NOTE: If left blank 0K is assumed even if the region is set to Langevin (rather than NVE)
			baseZeroAtomList: (Bool) If True our atomList is assumed to be base-zero, if false its assumed base 1. CP2K uses base one (i.e. the first atom in CP2K has the index 1)

		"""
		self.atomList = atomList
		self.doLangevin = doLangevin
		self.noisyGamma = noisyGamma
		self.baseZeroAtomList = baseZeroAtomList
		self.temperature = temperature

	@property
	def outAtomList(self):
		if self.baseZeroAtomList:
			return [x+1 for x in self.atomList]
		else:
			return self.atomList

	def addToPyCp2kObj(self, pyCp2kObj):
		pyCp2kObj.CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION_add()
		currRegion = pyCp2kObj.CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION_list[-1]
		currRegion.Do_langevin = self.doLangevin
		currRegion.List = self.outAtomList
		currRegion.Noisy_gamma_region = self.noisyGamma
		currRegion.Temperature = self.temperature

