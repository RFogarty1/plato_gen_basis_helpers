
import copy

class CP2KBasisObjBase():
	"""Class containing all information CP2K needs for the basis set of one element

	"""

	@property
	def element(self):
		""" String, element symbol (e.g. "Mg")
		"""
		raise NotImplementedError("")

	@property
	def basis(self):
		""" String representing the name of the basis set in CP2K (e.g. DZV-GTH) """
		raise NotImplementedError("")

	@property
	def potential(self):
		""" String representing the name of the pseudopotential in CP2K (e.g. GTH-PBE-q2)"""
		raise NotImplementedError("")

	@property
	def basisFile(self):
		""" String representing the name of the basis file this basis set comes from (e.g. GTH_BASIS_SETS) """
		raise NotImplementedError("")

	@property
	def potFile(self):
		""" String representing the name of the potential file this psuedopotential comes from (e.g. GTH_POTENTIALS) """
		raise NotImplementedError("")

	@property
	def ghost(self):
		""" String representing whether to treat this as a ghost atom; with basis functions but no electrons or core """
		raise NotImplementedError("")

	@property
	def kind(self):
		""" String representing the label for this type of atom. Usually the same as element, but sometimes its neccesary to treat different atoms of the same element differently hence we need this"""
		raise NotImplementedError("")

class CP2KBasisObjStandard(CP2KBasisObjBase):

	def __init__(self, element=None, basis=None, potential=None, basisFile=None, potFile=None, ghost=False, kind=None):
		""" Initializer for class designed to hold all required information for the CP2K basis set/potential to use for one element
		
		Args (ALL are required, despite being keywords):
			element: (str) Element symbol (e.g. Mg)
			basis: (str) String representing the name of the basis set in CP2K (e.g. DZV-GTH)
			potential: (str) String representing the name of the pseudopotential in CP2K (e.g. GTH-PBE-q2)
			basisFile: (str) String representing the name of the basis file this basis set comes from (e.g. GTH_BASIS_SETS)
			potFile: (str) String representing the name of the potential file this psuedopotential comes from (e.g. GTH_POTENTIALS)
			ghost: (Bool, Optional, Default=False) Whether to treat this atom type ghost atoms with basis functions but no core/electrons
			kind: (str, Optional, Defaults to element if unset) Essentially the label for this atom type, often (not always) the same as element

		Raises:
			AttributeError: If ANY of the keyword arguments arent set (Except ghost)
		"""
		reqArgs = ["element", "basis", "potential", "basisFile", "potFile","ghost", "kind"]
		self.allArgs = list(reqArgs)

		self._element = element
		self._basis = basis
		self._potential = potential
		self._basisFile = basisFile
		self._potFile = potFile
		self._ghost = ghost
		self._kind = str(element) if kind is None else kind

		#Check that everything is set
		for arg in reqArgs:
			if getattr(self,arg) is None:
				raise AttributeError("{} is a required argument".format(arg))

	@property
	def element(self):
		return self._element

	@element.setter
	def element(self,val):
		self._element = val

	@property
	def basis(self):
		return self._basis

	@property
	def potential(self):
		return self._potential

	@potential.setter
	def potential(self,val):
		self._potential = val

	@property
	def basisFile(self):
		return self._basisFile

	@basisFile.setter
	def basisFile(self,val):
		self._basisFile = val

	@property
	def potFile(self):
		return self._potFile

	@potFile.setter
	def potFile(self,val):
		self._potFile = val

	@property
	def ghost(self):
		return self._ghost

	@ghost.setter
	def ghost(self,val):
		self._ghost = val

	@property
	def kind(self):
		return self._kind

	@kind.setter
	def kind(self, val):
		self._kind = val

	def toDict(self):
		return {k: getattr(self,k) for k in self.allArgs}

	@classmethod
	def fromDict(cls, inpDict):
		return cls(**inpDict)

	def __eq__(self,other):
		for arg in self.allArgs:
			if getattr(self,arg) != getattr(other,arg):
				return False
		return True


def getBasisObjsWithGhostVersionsIncluded(basisObjIter):
	""" Gets an iter of basis objects with ghost atom basis sets included for each atom type
	
	Args:
		basisObjIter: (iter of CP2KBasisObjStandard)
			 
	Returns
		outObjs: (iter of CP2KBasisObjStandard) Same as basisObjIter except for every non-ghost atom type in the input, there is a ghost atom version in this output list (i.e. this list includes basisObjIter objs)
 
	"""
	extraObjs = list()
	for obj in basisObjIter:
		if obj.ghost is False:
			newObj = getStandardGhostVersionOfBasisObj(obj)
			extraObjs.append( newObj )
	
	return [x for x in basisObjIter] + extraObjs

def getStandardGhostVersionOfBasisObj(basisObj):
	""" Gets a version of the input basis obj with ghost=True and ele key modded to have _ghost attached
	
	Args:
		basisObj: (CP2KBasisObjStandard)
			 
	Returns
		 outBasisObj: (CP2KBasisObjStandard), this is always a DIFFERENT object (even if the input is already a ghost type, in which case (basisObj is not outBasisObj) BUT (basisObj==outBasisObj)
 
	"""
	outBasisObj = copy.deepcopy(basisObj)
	if outBasisObj.ghost is True:
		return outBasisObj
	
	outBasisObj.element = outBasisObj.element + "_ghost"
	outBasisObj.kind = outBasisObj.kind + "_ghost"
	outBasisObj.ghost = True

	return outBasisObj
