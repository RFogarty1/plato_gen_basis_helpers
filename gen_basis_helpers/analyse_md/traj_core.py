
import plato_pylib.shared.ucell_class as uCellHelp

class TrajectoryBase():
	"""Object representing the total trajectory of a simulation. Main job is to yield a generator/iterable of TrajStepBase objects

	"""

	def __iter__(self):
		raise NotImplementedError("")


class TrajectoryInMemory(TrajectoryBase):

	def __init__(self, trajSteps):
		""" Initializer
		
		Args:
			trajSteps: (iter of TrajStepBase) Each object represents one step of an MD run
				 
		"""
		self.trajSteps = trajSteps

	def __iter__(self):
		return iter(self.trajSteps)

	@classmethod
	def fromDict(cls, inpDict):
		trajSteps = [TrajStepBase.fromDict(x) for x in inpDict["trajSteps"]]
		return cls(trajSteps)

	def toDict(self):
		outDict = {"trajSteps":[x.toDict() for x in self.trajSteps]}
		return outDict

	def __eq__(self, other):
		if self.trajSteps!=other.trajSteps:
			return False
		return True


class TrajStepBase():
	""" Represents a single step in a trajectory. Used to define the minimal expected attributes

	Attributes:
		unitCell: (plato_pylib UnitCell) Contains geometry at this step (central cell atoms only)
		step: (int) Contains the step number
		time: (float) Contains the time at this step number
	"""
	def __init__(self, unitCell=None, step=None, time=None):
		self._eqTol = 1e-5
		self.normalCmpAttrs = ["unitCell", "step"]
		self.numericalCmpAttrs = ["time"]
		self.unitCell = unitCell
		self.step = step
		self.time = time

	@classmethod
	def fromDict(cls, inpDict):
		useDict = {k:v for k,v in inpDict.items() if k!="unitCell"}
		if inpDict.get("unitCell",None) is not None:
			useDict["unitCell"] = uCellHelp.UnitCell.fromDict(inpDict["unitCell"])
		return cls(**useDict)

	def toDict(self):
		outDict = {"step":self.step, "time":self.time, "unitCell":self.unitCell.toDict()}
		return outDict

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)

		for attr in self.normalCmpAttrs:
			if getattr(self,attr)!=getattr(other,attr):
				return False

		for attr in self.numericalCmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if (valA is None) and (valB is None):
				pass
			elif (valA is None) or (valB is None):
				return False
			elif abs(valA-valB)>eqTol:
				return False

		return True

