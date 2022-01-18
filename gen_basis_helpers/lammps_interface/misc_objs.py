
import collections

class TypeMaps():
	""" Simple class for holding mappings from elements to indices. See initializer for attributes

	"""
	def __init__(self, eleToTypeIdx=None, bondToTypeIdx=None, angleToTypeIdx=None):
		""" Initializer
		
		Args:
			eleToTypeIdx: (dict) Keys are elements (e.g. "Mg") while values are indices (e.g. 1, 2) used to refer to them in LAMMPS files
			bondToTypeIdx: (dict) Keys are pairs of elements (e.g. ("O","H")) while values are "bond type" indices used in LAMMPS files
			angleToTypeIdx: (dict) Keys are len-3 tuples of elemets (e.g ("H","O","H")) while values are "angle type" indices used in LAMMPS files
 
		"""
		self.eleToTypeIdx = eleToTypeIdx
		self.bondToTypeIdx = bondToTypeIdx
		self.angleToTypeIdx = angleToTypeIdx



class Ensemble():
	"""Represents an ensemble, e.g. NVT, NVE, NPT"""

	@property
	def fixStr(self):
		raise NotImplementedError("")


class NVEEnsemble(Ensemble):

	def __init__(self):
		pass

	@property
	def fixStr(self):
		return "all nve"

class NVTEnsemble(Ensemble):

	def __init__(self, startTemp, thermostat="Nose-Hoover", endTemp=None, dampTime=None, numbFmt="{:.1f}", seed=1):
		""" Initializer
		
		Args:
			startTemp: (Float) Target temperature at the start of the run
			endTemp: (Float,Optional) Target temperature at the end of the run. Default is startTemp (meaning run the whole simulation at the same temperature)
			thermostat: (Str) The type of thermostat to use
			dampTime: (Float) Parameter represents ~the time expected to go from starting temperature (determined by initial velocities) to the target value. 100x timestep is generally recommended
			numbFmt: (Str) Reperesents the formatting to use for writing numerical values
			seed: (positive int) The seed for any random number generation
 
		"""
		self.startTemp = startTemp
		self.thermostat = thermostat
		self.dampTime = dampTime
		self.endTemp = endTemp if endTemp is not None else startTemp
		self.numbFmt = numbFmt
		self.seed = seed

	@property
	def fixStr(self):
		self._checkAttrsValid()
		thermoStr = self._getThermostatStr()
		if self.thermostat.lower()=="Nose-Hoover".lower():
			outFmt = "all {} temp " + " ".join([self.numbFmt for x in range(3)])
		elif self.thermostat.lower()=="Langevin".lower():
			outFmt = "all {} " + " ".join([self.numbFmt for x in range(3)]) + " {}".format(self.seed)
		else:
			raise NotImplementedError("")

		return outFmt.format(thermoStr, self.startTemp, self.endTemp, self.dampTime)

	def _checkAttrsValid(self):
		if self.dampTime is None:
			raise ValueError("{} is an invalid value for dampTime".format(self.dampTime))

	def _getThermostatStr(self):
		if self.thermostat.lower()=="Nose-Hoover".lower():
			return "nvt"
		elif self.thermostat.lower()=="Langevin".lower():
			return "langevin"
		else:
			raise ValueError("{} is an invalid value for thermostat".format(self.thermostat))



class NPTEnsembleNoseHooverStandard(Ensemble):

	def __init__(self, startTemp, startPressure, pressureDims=None, endTemp=None,  endPressure=None, dampTimeTemp=None, dampTimePressure=None,
	             numbFmtTime="{:.1f}", numbFmtPressure="{:.1f}", numbFmtTemp="{:.1f}"):
		self.startTemp = startTemp
		self.startPressure = startPressure
		self.pressureDims = "iso" if pressureDims is None else pressureDims
		self.endTemp = startTemp if endTemp is None else endTemp
		self.endPressure = startPressure if endPressure is None else endPressure
		self.dampTimeTemp = dampTimeTemp
		self.dampTimePressure = dampTimePressure
		self.numbFmtTime = numbFmtTime
		self.numbFmtPressure = numbFmtPressure
		self.numbFmtTemp = numbFmtTemp

	def _checkAttrsValid(self):
		if self.dampTimePressure is None:
			raise ValueError("{} is an invalid value for dampTimePressure".format(self.dampTimePressure))
		if self.dampTimeTemp is None:
			raise ValueError("{} is an invalid value for dampTimeTemp".format(self.dampTimeTemp))

	@property
	def fixStr(self):
		self._checkAttrsValid()
		startTempStr, endTempStr = [self.numbFmtTemp.format(x) for x in [self.startTemp, self.endTemp]]
		startPressStr, endPressStr = [self.numbFmtPressure.format(x) for x in [self.startPressure,self.endPressure]]
		dampTimePress, dampTimeTemp = [self.numbFmtTime.format(x) for x in [self.dampTimePressure,self.dampTimeTemp]]
		outStrFmt = "all npt temp {sTemp:} {eTemp:} {dampTemp:} {pressDims:} {sPress:} {ePress:} {dampPress:}"
		outKwargs = {"sTemp":startTempStr, "eTemp":endTempStr, "dampTemp":dampTimeTemp, "pressDims":self.pressureDims,
		             "sPress":startPressStr, "ePress":endPressStr, "dampPress":dampTimePress }
		return outStrFmt.format(**outKwargs)


class VelocityCreateCommObj():

	def __init__(self, temp, seed=200, dist="uniform", group="all", numbFmt="{:.1f}"):
		""" Initializer
		
		Args:
			temp: (float) Temperature that the initial velocities correspond to
			seed: (int) A seed for the random number generator; explicitly setting allows more reproducible settings
			dist: (str) Distribution to use.
			group: (str) The group of atoms to apply this command to
			numbFmt: (str) Formatting string for floats (i.e. temperature value)
				 
		"""
		self.temp = temp
		self.seed = seed
		self.dist = dist
		self.group = group
		self.numbFmt = numbFmt

	@property
	def commandStr(self):
		outFmt = "{} create " + self.numbFmt + " {} dist {}"
		outStr = outFmt.format(self.group, self.temp, self.seed, self.dist)
		return outStr



class DumpCommObjStandard():

	def __init__(self, everyNSteps, groupId="all",dumpType="atom", fileExt="lammpstrj", scale=False):
		""" Initializer
		
		Args:
			everyNSteps: (int) Write to the dump file for every N steps
			groupId: (str) Id for the group of atoms
			dumpType: (str) Built in type for the dump file
			fileExt: (str) The extension to use for the dump file
			scale: (Bool) Whether to use scaled co-ordinates or not
				 
		"""
		self.everyNSteps = everyNSteps
		self.groupId = groupId
		self.dumpType = dumpType
		self.fileExt = fileExt
		self.scale = scale

	@property
	def commandDict(self):
		outDict = collections.OrderedDict()
		scaleStr = "no" if self.scale is False else "yes"

		dumpComm = "myDump {} {} {} {}".format(self.groupId, self.dumpType, self.everyNSteps, "dump."+self.fileExt)

		outDict["dump"] = dumpComm
		outDict["dump_modify"] = "myDump scale {}".format(scaleStr)

		return outDict


class ReflectiveWallFace():
	""" Class representing a reflective wall. This is placed on one face and causes particles to bounce off it
	"""
	def __init__(self, face, groupId="all"):
		""" Initializer
		
		Args:
			face: (str) Either xlo, xhi, ylo, yhi, zlo, zhi.
			groupId: (str) The group command passed to face. Presumably its possible to make this wall only work on certain sets of atoms

		"""
		self.face = face
		self.groupId = groupId


	@property
	def fixStr(self):
		outStr = "{} wall/reflect {} EDGE".format(self.groupId, self.face)
		return outStr


class CombChargeNeutralisationOpts():
	""" Class representing charge neutralisation options for COMB potentials

	"""

	def __init__(self, groupId="all", nEvery=100, precision=1e-3):
		""" Initializer
		
		Args:
			groupId: (str) Identifier for the relevant group; default of all means charge-neutralisation is applied to all atoms
			nEvery: (int) Frequency to carry out charge neutralisation (nEvery=100 means do it every 100 steps) 
			precision: (float) Parameter controlling numerical accuracy to carry out charge neutralisation to
 
		"""
		self.groupId = groupId
		self.nEvery = nEvery
		self.precision = precision

	@property
	def fixStr(self):
		outStr = "{} qeq/comb {} {:.5f}".format(self.groupId, self.nEvery, self.precision)
		return outStr


class RescaleVelocitiesSimple():
	""" Class representing options for periodically rescaling velocities IFF temperature goes out of selected range; this is mainly to deal with early steps where T can spike """

	def __init__(self, groupId="all", nEvery=1, tempStart=300, tempEnd=None, maxDeviation=200, fraction=1):
		""" Initializer
		
		Args:
			groupId: (str) Identifier for the relevant group; default of all means we treat all atoms as the group
			nEvery: (int) Check to rescale every N steps
			tempStart: (float) Temperature at the start of the run
			tempEnd: (float) Temperature at the end of the run; Default is to use the same as the start temperature
			maxDeviation: (float) We rescale velocities if the temperature is this far from where the average will be
			fraction: (float) Rescale velocities to get this close to expected average temperature
 
		"""
		self.groupId = groupId
		self.nEvery = nEvery
		self.tempStart = tempStart
		self.tempEnd = tempStart if tempEnd is None else tempEnd
		self.maxDeviation = maxDeviation
		self.fraction = fraction

	@property
	def fixStr(self):
		outFmt = "{} temp/rescale {:d} {:.2f} {:.2f} {:.2f} {:.2f}"
		currArgs = [ self.groupId, self.nEvery, self.tempStart, self.tempEnd, self.maxDeviation, self.fraction ]
		return outFmt.format(*currArgs)



class AtomGroup():
	""" Represents a Group of atoms. These are used in LAMMPS input files to apply certain commands ONLY to certain groups

	"""

	@property
	def groupStr(self):
		raise NotImplementedError("")



class AtomGroupByType(AtomGroup):

	def __init__(self, groupId, typeIndices):
		""" Initializer
		
		Args:
			groupId: (str) name given to identify this group; essentially an attribute name within the lammps script file
			typeIndices: (iter of ints) Indices corresponding to the atom groups you want to include in this overall group (e.g. it may be "2" and "3" to group all water together, if O/H have ids of 2/3) 
				 
		"""
		self.groupId = groupId
		self.typeIndices = typeIndices

	@property
	def groupStr(self):
		return "{} type ".format(self.groupId) + " ".join([str(x) for x in self.typeIndices])
	





