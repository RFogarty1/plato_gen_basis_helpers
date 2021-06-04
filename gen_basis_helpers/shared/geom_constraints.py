import itertools as it


#TODO: Certain combinations of atomic/cell constraints wont be possible (e.g. free cell optimisation with fixed CARTESIAN co-ordinates); may need something to check that
class GeomConstraints():
	""" Class containing information on which geometrical parameters require constraining to initial values

	"""
	def __init__(self, atomicPositionConstraints, cellConstraints):
		self.atomicPositionConstraints = atomicPositionConstraints
		self.cellConstraints = cellConstraints

	@property
	def constraintsPresent(self):
		return self.atomicPositionConstraints.constraintsPresent or self.cellConstraints.constraintsPresent

	@classmethod
	def initWithNoConstraints(cls):
		atomicPosConstraints = AtomicPositionConstraints.initWithNoConstraints()
		cellConstraints = CellConstraints.initWithNoConstraints()
		return cls(atomicPosConstraints, cellConstraints)

	def __eq__(self, other):

		#Generally will be the faster comparison
		if self.cellConstraints != other.cellConstraints:
			return False

		if self.atomicPositionConstraints != other.atomicPositionConstraints:
			return False

		return True

	@classmethod
	def fromDict(cls, inpDict):
		atomicPosConstraints = AtomicPositionConstraints.fromDict( inpDict["atomicPositionConstraints"] )
		cellConstraints = CellConstraints.fromDict( inpDict["cellConstraints"] )
		return cls(atomicPosConstraints, cellConstraints)

	def toDict(self):
		return {"atomicPositionConstraints": self.atomicPositionConstraints.toDict(),
		        "cellConstraints": self.cellConstraints.toDict()}


class CellConstraints():
	"""Class containing information on which unit cell paramters require constraining to initial values

	Attributes:
		anglesToFix: (len-3 bool iter) Each element is whether to fix a lattice angle in order [alpha, beta, gamma]. True means fix to initial value
		lattParamsToFix: (len-3 bool iter) Each element is whether to fix a lattice parameter in order [a,b,c]. True means fix to initial value

	"""

	def __init__(self, anglesToFix, lattParamsToFix):
		""" Initialiser
		
		Args:
			anglesToFix: (len-3 bool iter) Each element is whether to fix a lattice angle in order [alpha, beta, gamma]. True means fix to initial value
			lattParamsToFix: (len-3 bool iter) Each element is whether to fix a lattice parameter in order [a,b,c]. True means fix to initial value
				 
		"""

		self.listedAttrs = ["anglesToFix", "lattParamsToFix"]
		self.anglesToFix = list(anglesToFix)
		self.lattParamsToFix = list(lattParamsToFix)
		self._checkParamsCorrectLength()

	def _checkParamsCorrectLength(self):
		if len(self.anglesToFix) != 3:
			raise AttributeError(" Len(3) bool iter required for anglesToFix; {} is not allowed".format(self.anglesToFix) )
		if len(self.lattParamsToFix) != 3:
			raise AttributeError(" Len(3) bool iter required for lattParamsToFix; {} is not allowed".format(self.lattParamsToFix) )


	@property
	def constraintsPresent(self):
		""" If constraints on cell are present, then this returns True, else False.	
		"""
		if any([x is True for x in self.anglesToFix]):
			return True

		if any([x is True for x in self.lattParamsToFix]):
			return True

		return False

	@classmethod
	def initWithNoConstraints(cls):
		anglesToFix = [False,False,False]
		lattParamsToFix = [False,False,False]
		return cls(anglesToFix, lattParamsToFix)

	@classmethod
	def fromDict(cls, inpDict):
		anglesToFix, lattParamsToFix = inpDict["anglesToFix"], inpDict["lattParamsToFix"]
		return cls(anglesToFix,lattParamsToFix)

	def toDict(self):
		outDict = dict()
		for attr in self.listedAttrs:
			outDict[attr] = getattr(self,attr)
		return outDict

	def __eq__(self,other):
		boolAttrs = ["anglesToFix","lattParamsToFix"]
		for attr in boolAttrs:
			if getattr(self,attr) != getattr(other,attr):
				return False
		return True



class AtomicPositionConstraints():
	""" This will eventually be used to hold information on which atomic positions are fixed to initial value (including fixing things like bond lengths); currently just a stub

	"""
	def __init__(self, atomicCartConstraints=None):
		self.atomicCartConstraints = list() if atomicCartConstraints is None else list(atomicCartConstraints)

	@property
	def constraintsPresent(self):
		if any([x.constraintsPresent for x in self.atomicCartConstraints]):
			return True
		return False

	@classmethod
	def initWithNoConstraints(cls):
		return cls()

	#Note: Slightly restrictive. If the same constraints are in different orders, then objects compare unequal
	def __eq__(self, other):
		#Filter out blank constraints
		atConstrA = [x for x in self.atomicCartConstraints if x.constraintsPresent]
		atConstrB = [x for x in other.atomicCartConstraints if x.constraintsPresent]

		if len(atConstrA) != len(atConstrB):
			return False

		for constrA,constrB in it.zip_longest(atConstrA,atConstrB):
			if constrA != constrB:
				return False

		return True

	@classmethod
	def fromDict(cls, inpDict):
		atomicCartConstraints = [AtomicCartesianConstraint.fromDict(x) for x in inpDict["atomicCartConstraints"]]
		return cls(atomicCartConstraints=atomicCartConstraints)

	def toDict(self):
		outDicts = [x.toDict() for x in self.atomicCartConstraints]
		return {"atomicCartConstraints":outDicts}

class AtomicCartesianConstraint():
	""" Represents constraints to apply to cartesian co-ordinates of a single atom

	"""

	def __init__(self, atomIdx, fixX=False, fixY=False, fixZ=False):
		""" Initializer
		
		Args:
			atomIdx (int): Index of the atom to apply constraints to, indexing starts at zero
				 
		"""
		self.atomIdx = atomIdx
		self.fixX = fixX
		self.fixY = fixY
		self.fixZ = fixZ

	def __eq__(self,other):
		directCmpAttrs = ["atomIdx","fixX","fixY","fixZ"]

		for attr in directCmpAttrs:
			if getattr(self,attr) != getattr(other,attr):
				return False

		return True

	@property
	def constraintsPresent(self):
		constraintsAttrs = ["fixX","fixY","fixZ"]
		for x in constraintsAttrs:
			if getattr(self,x) is True:
				return True
		return False

	@classmethod
	def fromDict(cls, inpDict):
		atomIdx = inpDict["atomIdx"]
		fixX, fixY, fixZ = [ inpDict[attr] for attr in ["fixX", "fixY", "fixZ"] ]
		kwargs = {"fixX":fixX, "fixY":fixY, "fixZ":fixZ}
		return cls(atomIdx, **kwargs)

	def toDict(self):
		attrs = ["atomIdx", "fixX", "fixY", "fixZ"]
		outDict = {attr: getattr(self,attr) for attr in attrs}
		return outDict




