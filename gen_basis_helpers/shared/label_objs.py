

class BaseLabel():

	@property
	def components(self):
		""" List of components contributing to identity of this label. """
		raise NotImplementedError("")

	@property
	def labelNames(self):
		""" List of the attrs used to define the label """
		raise NotImplementedError("")

	def __eq__(self, other):
		for attr in self.labelNames:
			try:
				if getattr(self,attr) != getattr(other,attr):
					return False
			except AttributeError: #If other is an inappropriate type (e.g an integer) this should trigger
				return False
		return True

	def __hash__(self):
		return hash( tuple([getattr(self,x) for x in self.labelNames]) )



class StandardLabel(BaseLabel):
	"""Label holding strings used to differentiate different objects

	Attributes:
		eleKey (str): key representing the element (or possible elements) involved.
		structKey (str): key representing the structure involved
		methodKey (str): key representing the calculation method used

	"""
	def __init__(self, eleKey=None, structKey=None, methodKey=None):

		setattr(self, "eleKey", eleKey)
		setattr(self, "structKey", structKey)
		setattr(self, "methodKey", methodKey)

		self.reqArgs = ["eleKey", "structKey", "methodKey"]
		for x in self.reqArgs:
			if getattr(self,x) is None:
				raise ValueError("{} is a required argument for StandardLabel".format(x))

	@property
	def labelNames(self):
		return self.reqArgs


	@property
	def components(self):
		""" List of components contributing to identity of this label. At time of writing this is eleKey, structKey and methodKey. Consistent ordering not garanteed.
		
		"""
		return [getattr(self,x) for x in self.labelNames]


	def __repr__(self):
		outDict = {k:getattr(self,k)  for k in self.labelNames}
		return str(outDict)



