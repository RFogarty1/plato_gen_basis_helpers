
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


class NVTEnsemble(Ensemble):

	def __init__(self, startTemp, thermostat="Nose-Hoover", endTemp=None, dampTime=None, numbFmt="{:.1f}"):
		""" Initializer
		
		Args:
			startTemp: (Float) Target temperature at the start of the run
			endTemp: (Float,Optional) Target temperature at the end of the run. Default is startTemp (meaning run the whole simulation at the same temperature)
			thermostat: (Str) The type of thermostat to use
			dampTime: (Float) Parameter represents ~the time expected to go from starting temperature (determined by initial velocities) to the target value. 100x timestep is generally recommended
			numbFmt: (Str) Reperesents the formatting to use for writing numerical values
 
		"""
		self.startTemp = startTemp
		self.thermostat = thermostat
		self.dampTime = dampTime
		self.endTemp = endTemp if endTemp is not None else startTemp
		self.numbFmt = numbFmt

	@property
	def fixStr(self):
		self._checkAttrsValid()
		thermoStr = self._getThermostatStr()
		outFmt = "all {} temp " + " ".join([self.numbFmt for x in range(3)])
		return outFmt.format(thermoStr, self.startTemp, self.endTemp, self.dampTime)

	def _checkAttrsValid(self):
		if self.dampTime is None:
			raise ValueError("{} is an invalid value for dampTime".format(self.dampTime))

	def _getThermostatStr(self):
		if self.thermostat.lower()=="Nose-Hoover".lower():
			return "nvt"
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

