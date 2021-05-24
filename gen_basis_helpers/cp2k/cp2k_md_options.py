
import itertools as it

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

	def __init__(self, metaVars, ntHills=None, doHills=True, hillHeight=None, printColvarCommonIter=3,
	             heightNumbFmt="{:.4f}", printHills=True, printHillsCommonIter=3, spawnHillsOpts=None):
		""" Initializer
		
		Args:
			metaVars: (iter of MetaVarStandard objects)
			ntHills: (int) create a hill at most every N steps
			doHills: (Bool) True means create hills during the simulation; basically should always be true (maybe false when reading hills in)
			hillHeight: (float) Hill height in Ha (default CP2K units for this)
			heightNumbFmt: (str) format string for converting hillHeight into a string 
			printHills: (Bool) True means cp2k prints the hills created during the simulation
			printHillsCommonIter: (int) Determines how many values are written in a single output file. LEAVE AS DEFAULT
			spawnHillsOpts: (MetadynamicsSpawnHillsOptions) Object specifying options related to spawninig hills initially

		"""
		self.metaVars = metaVars
		self.ntHills = ntHills
		self.doHills = doHills
		self.hillHeight = hillHeight
		self.printColvarCommonIter = printColvarCommonIter
		self.heightNumbFmt = heightNumbFmt
		self.printHills = printHills
		self.printHillsCommonIter = printHillsCommonIter
		self.spawnHillsOpts = spawnHillsOpts

	@property
	def optDict(self):
		outDict = {"metaVars":self.metaVars, "metaDyn_doHills":self.doHills, 
		           "metaDyn_hillHeight": self.heightNumbFmt.format(self.hillHeight) if self.hillHeight is not None else None,
		           "metaDyn_printColvarCommonIterLevels": self.printColvarCommonIter,
		           "metaDyn_ntHills":self.ntHills, "metaDyn_printHills":self.printHills,
		           "metaDyn_printHillsCommonIterLevels":self.printHillsCommonIter}

		outDict = {k:v for k,v in outDict.items() if v is not None}
		if self.spawnHillsOpts is not None:
			outDict.update(self.spawnHillsOpts.optDict)

		return outDict

class MetadynamicsSpawnHillsOptions():
	""" Class used to hold info for spawning hills, and to return it in a useful format to modify a pycp2k obj """

	def __init__(self, hillDicts):
		""" Initializer
		
		Args:
			hillDicts: (iter of dicts)

		hillDict:
			Keys are "scale", "pos", "height". Values for scale and pos need to be an iter (even if length 1) while height is just a float. Units are whatever cp2k uses in restart files (i think positions are generally in bohr when its a length, while heights are generally in Hartree (I THINK). Regardless, if these were read from a metadynLog file use the same units there
				 
		"""
		self.hillDicts = hillDicts
		self._eqTol = 1e-6

	@classmethod
	def fromIters(cls, scales=None, heights=None, positions=None):
		""" Alternative initializer
		
		Args: (ALL need to be set, despite the fact their keyword args)
			scales: (iter of float-iters) Contains width-parameter of hills
			heights: (iter of floats) Contains height-paramter of hills
			positions: (iter of float-iters) Contains the positions of each hill 
				 
		"""
		assert all([x is not None for x in [scales,heights,positions]])
		assert all([len(x)==len(scales) for x in [scales,heights,positions]])
		outDicts = list()
		for scale, height, pos in it.zip_longest(scales, heights, positions):
			currDict = {"scale":scale, "pos":pos, "height":height}
			outDicts.append(currDict)
		return cls(outDicts)

	@classmethod
	def fromMetadynHillInfo(cls, metaHillInfo, allowDiffHeights=False, heightTol=1e-3):
		""" Alternative initializer
		
		Args:
			metaHillInfo: (MetadynHillsInfo instance) This contains all info on spawned hill, including their time, and is generally how they are stored/manipulated by me at time of writing
			allowDiffHeights: (Bool) Only False implemented for now; meaning for each hill each collective variable must have the same height parameter
			heightTol: (float) For a single hill, the height for each collective variable must be within this to not raise an error when allowDiffHeights=False
 
		"""
		scales, positions = metaHillInfo.scales, metaHillInfo.positions
		#TODO:  Check heights are consistent and raise if not
		heights = [x[0] for x in metaHillInfo.heights]
		return cls.fromIters(scales=scales, heights=heights, positions=positions)


	@property
	def optDict(self):
		heights, scales, positions = list(), list(), list()
		for hillDict in self.hillDicts:
			heights.append( hillDict["height"] )
			scales.append( hillDict["scale"] )
			positions.append( hillDict["pos"] )

		outDict = {"metaDyn_spawnHillsHeight":heights, "metaDyn_spawnHillsPos":positions, "metaDyn_spawnHillsScale":scales,
		           "metaDyn_nHillsStartVal": len(heights)}
		return outDict

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)
		if len(self.hillDicts) != len(other.hillDicts):
			return False

		floatIterKeys = ["height"]
		iterOfFloatIterKeys = ["scale","pos"]
		
		for hillA,hillB in zip(self.hillDicts, other.hillDicts):
			for attr in iterOfFloatIterKeys:
				if len(hillA[attr]) != len(hillB[attr]):
					return False
				for valA, valB in zip( hillA[attr], hillB[attr] ):
					if abs(valA-valB) > eqTol:
						return False

			for attr in floatIterKeys:
				valA, valB = hillA[attr], hillB[attr]
				if abs(valA-valB) > eqTol:
					return False

		return True

class ThermostatOptsMultiRegions(ThermostatOptsBase):
	""" Extension of ThermostatOptsBase that allows multiple regions to be defined """
	def __init__(self, basicThermoOpts, regions, regionKwarg="defined"):
		""" Initializer
		
		Args:
			basicThermoOpts: (ThermostatOptsBase object) e.g. NoseThermostatOpts. Contains all options except those related to the regions to define
			regions: (iter of ThermostatRegionInfo) Contains objects representing thermostat regions

		"""
		self.basicThermoOpts = basicThermoOpts
		self.regions = regions
		self.regionKwarg = regionKwarg

	def addToPyCp2kObj(self, pyCp2kObj):
		self.basicThermoOpts.addToPyCp2kObj(pyCp2kObj)
		for reg in self.regions:
			reg.addToPyCp2kObj(pyCp2kObj)
		pyCp2kObj.CP2K_INPUT.MOTION.MD.THERMOSTAT.Region = self.regionKwarg


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


class AdaptiveLangevinThermostatOpts(ThermostatOptsBase):

	def __init__(self, timeConLangevin=None, timeConNose=None):
		""" Initializer
		
		Args:
			timeConLangevin: (float) Time constant (in fs) for the Langevin part of the thermostat (inverse friction coefficient)
			timeConNose: (float) Time constant (in fs) for the Nose part of the thermostat
 
		"""
		self.timeConLangevin = timeConLangevin
		self.timeConNose = timeConNose

	def addToPyCp2kObj(self, pyCp2kObj):
		thermostatSection = pyCp2kObj.CP2K_INPUT.MOTION.MD.THERMOSTAT
		thermostatSection.Type = "AD_LANGEVIN"
		thermostatSection.AD_LANGEVIN.Timecon_langevin = self.timeConLangevin
		thermostatSection.AD_LANGEVIN.Timecon_nh = self.timeConNose



class _OutAtomListMixin():

	@property
	def outAtomList(self):
		if self.baseZeroAtomList:
			return [x+1 for x in self.atomList]
		else:
			return self.atomList

class ThermalRegion(_OutAtomListMixin):
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

	def addToPyCp2kObj(self, pyCp2kObj):
		pyCp2kObj.CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION_add()
		currRegion = pyCp2kObj.CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION_list[-1]
		currRegion.Do_langevin = self.doLangevin
		currRegion.List = self.outAtomList
		currRegion.Noisy_gamma_region = self.noisyGamma
		currRegion.Temperature = self.temperature


class ThermostatRegionInfo(_OutAtomListMixin):
	""" Class representing a defined thermostat region under CP2K_INPUT/MOTION/MD/THERMOSTAT"""

	def __init__(self, atomList=None, baseZeroAtomList=True):
		""" Initializer
		
		Args:
			atomList: (iter of ints)
			baseZeroAtomList: (Bool) If True our atomList is assumed to be base-zero, if false its assumed base 1. CP2K uses base one (i.e. the first atom in CP2K has the index 1)
 
		"""
		self.atomList = atomList
		self.baseZeroAtomList = baseZeroAtomList

	def addToPyCp2kObj(self, pyCp2kObj):
		thermostatSection = pyCp2kObj.CP2K_INPUT.MOTION.MD.THERMOSTAT
		thermostatSection.DEFINE_REGION_add()
		thermostatSection.DEFINE_REGION_list[-1].List = self.outAtomList
