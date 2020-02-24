


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



class CP2KBasisObjStandard(CP2KBasisObjBase):

	def __init__(self, element=None, basis=None, potential=None, basisFile=None, potFile=None):
		""" Initializer for class designed to hold all required information for the CP2K basis set/potential to use for one element
		
		Args (ALL are required, despite being keywords):
			element: (str) Element symbol (e.g. Mg)
			basis: (str) String representing the name of the basis set in CP2K (e.g. DZV-GTH)
			potential: (str) String representing the name of the pseudopotential in CP2K (e.g. GTH-PBE-q2)
			basisFile: (str) String representing the name of the basis file this basis set comes from (e.g. GTH_BASIS_SETS)
			potFile: (str) String representing the name of the potential file this psuedopotential comes from (e.g. GTH_POTENTIALS)
		
		Raises:
			AttributeError: If ANY of the keyword arguments arent set
		"""
		reqArgs = ["element", "basis", "potential", "basisFile", "potFile"]
		self.allArgs = list(reqArgs)

		self._element = element
		self._basis = basis
		self._potential = potential
		self._basisFile = basisFile
		self._potFile = potFile

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

	def __eq__(self,other):
		for arg in self.allArgs:
			if getattr(self,arg) != getattr(other,arg):
				return False
		return True
