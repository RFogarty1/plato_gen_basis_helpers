

#TODO: Add alternative initialiser for NO CONSTRAINTS to all these objects; which provides a natural way of specifying free-optimisations
class GeomConstraints():
	""" Class containing information on which geometrical parameters require constraining to initial values

	"""
	def __init__(self, atomicPostionConstraints, cellConstraints):
		self.atomicPostionConstraints = atomicPostionConstraints
		self.cellConstraints = cellConstraints

	@property
	def constraintsPresent(self):
		return self.atomicPostionConstraints.constraintsPresent or self.cellConstraints.constraintsPresent

	@classmethod
	def initWithNoConstraints(cls):
		atomicPosConstraints = AtomicPositionConstraints.initWithNoConstraints()
		cellConstraints = CellConstraints.initWithNoConstraints()
		return cls(atomicPosConstraints, cellConstraints)


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


class AtomicPositionConstraints():
	""" This will eventually be used to hold information on which atomic positions are fixed to initial value (including fixing things like bond lengths); currently just a stub

	"""
	def __init__(self):
		pass

	@property
	def constraintsPresent(self):
		return False

	@classmethod
	def initWithNoConstraints(cls):
		return cls()

