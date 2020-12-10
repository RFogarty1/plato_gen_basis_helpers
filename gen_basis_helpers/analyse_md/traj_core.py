
import json
import os
import pathlib
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

	def applyFunctToEachTrajStep(self, funct):
		""" Applies a function to each trajectory step. Useful for doing things like converting length units
		
		Args:
			funct: f(TrajStepBase) Modifies the object in place
				 
		"""
		for tStep in self.trajSteps:
			funct(tStep)

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


def dumpTrajObjToFile(trajObj, outFile):
	""" Dump TrajectoryBase to file. Format involves writing each step as a dict (i.e. JSON notation).
	
	Args:
		trajObj: (TrajectoryBase object) Contains all steps in the trajectory
		outFile: (str) Path to the output file
	 
	"""
	outDir = os.path.split(outFile)[0]
	pathlib.Path(outDir).mkdir(parents=True, exist_ok=True)
	with open(outFile,"wt") as f:
		for step in trajObj:
			f.write(json.dumps(step.toDict()))
			f.write("\n")

def readTrajObjFromFileToTrajectoryInMemory(inpFile):
	""" Reads a trajectory file (format defined by dumpTrajObjToFile) and returns a TrajectoryInMemory object
	
	Args:
		inpFile: (str) Path to file containing the trajectory
			 
	Returns
		 trajObj: (TrajectoryInMemory) Holds the full trajectory in memory. Simple/flexible but may not be ideal for large simulations

	"""
	outTrajObjs = list()
	with open(inpFile,"rt") as f:
		for line in f:
			currDict = json.loads(line)
			currObj = TrajStepBase.fromDict(currDict)
			outTrajObjs.append(currObj)

	return TrajectoryInMemory(outTrajObjs)

def readLastTrajStepFromFile(inpFile):
	""" Reads only the final trajectory step from a trajectory file (format defined by dumpTrajObjToFile). Faster than reading the whole thing and taking the last step
	
	Args:
		inpFile: (str) Path to the *.traj file
			 
	Returns
		outStep: (TrajStepBase) Contains information (including geometry) for the final step in an MD trajectory
 
	"""
	with open(inpFile,"rt") as f:
		line = getFinalNLinesFromFileObj(f)[-1]
		currDict = json.loads(line)
		outObj = TrajStepBase.fromDict(currDict)
	return outObj


#Taken from:
#https://stackoverflow.com/questions/136168/get-last-n-lines-of-a-file-similar-to-tail
#Original function was called tail
def getFinalNLinesFromFileObj(f, lines=1, _buffer=4098):
    """Tail a file and get X lines from the end"""
    # place holder for the lines found
    lines_found = []

    # block counter will be multiplied by buffer
    # to get the block size from the end
    block_counter = -1

    # loop until we find X lines
    while len(lines_found) < lines:
        try:
            f.seek(block_counter * _buffer, os.SEEK_END)
        except IOError:  # either file is too small, or too many lines requested
            f.seek(0)
            lines_found = f.readlines()
            break

        lines_found = f.readlines()

        # we found enough lines, get out
        # Removed this line because it was redundant the while will catch
        # it, I left it for history
        # if len(lines_found) > lines:
        #    break

        # decrement the block counter to get the
        # next X bytes
        block_counter -= 1

    return lines_found[-lines:]


