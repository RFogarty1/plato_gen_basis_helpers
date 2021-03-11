
from . import shared_io as sharedIoHelp

class NudgedBandPathStandard():

	def __init__(self, steps):
		""" Initializer
		
		Args:
			steps: (iter of NudgedBandStepStandard) Each represents a single step in the nudged band pathway
				 
		"""
		self.steps = list(steps)

	@classmethod
	def fromDict(cls, inpDict):
		return cls([ NudgedBandStepStandard.fromDict(x) for x in inpDict["steps"] ])

	def toDict(self):
		return {"steps": [x.toDict() for x in self.steps]}

	def __eq__(self, other):
		if len(self.steps) != len(other.steps):
			return False

		for stepA,stepB in zip(self.steps, other.steps):
			if stepA!=stepB:
				return False

		return True

class NudgedBandStepStandard():

	def __init__(self, geom=None, energies=None, dist=None):
		""" Iniitalizer
		
		Args:
			geom: (plato_pylib UnitCell object)
			energies: (plato_pylib energies object)
			dist: (float) Distance along the pathway from the START (so zero for the start geom)
 
		Returns
			What Function Returns
	 
		Raises:
			Errors
		"""
		self._eqTol = 1e-5
		self.geom = geom
		self.energies = energies
		self.dist = dist

	@classmethod
	def fromDict(cls, inpDict):
		return cls(**inpDict)

	def toDict(self):
		outAttrs = ["geom", "energies", "dist"]
		outDict = dict()
		for attr in outAttrs:
			outDict[attr] = getattr(self, attr)
		return outDict

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)

		numbAttrs = ["dist"]
		directCmpAttrs = ["geom", "energies"]

		#Compare numerical attributes (make sure their within error tolerance)
		for attr in numbAttrs:
			valA, valB = getattr(self, attr), getattr(other, attr)
			if (valA is None) and (valB is not None):
				return False
			elif (valA is not None) and (valB is None):
				return False
			elif (valA is None) and (valB is None):
				pass
			elif abs(valA-valB)>eqTol:
				return False

		#Compare attributes which only need a direct a==b test
		for attr in directCmpAttrs:
			valA, valB = getattr(self, attr), getattr(other,attr)
			if valA != valB:
				return False

		return True



def dumpNudgedBandPathwayToFile(nebPathway, outFile):
	""" Dump a nudged elastic band pathway to an output file in standard format (which is json)
	
	Args:
		nebPathway: (NudgedBandPathStandard)
		outFile: (str)
			 
 
	"""
	sharedIoHelp.dumpObjWithToDictToJson(nebPathway, outFile)


def readNudgedBandPathwayFromJsonFileStandard(inpFile):
	""" Read in a nudge band pathway from a dump in the standard format
	
	Args:
		inpFile: (str) Path containing the dumped NEB pathway
			 
	Returns
		outPath: (NudgedBandPathStandard) Object containing the NEB pathway
 
	"""
	return sharedIoHelp.readObjWithFromDictFromJsonFile(NudgedBandPathStandard,inpFile)

