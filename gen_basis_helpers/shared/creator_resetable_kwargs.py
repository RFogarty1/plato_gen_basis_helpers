
import contextlib

class CreatorWithResetableKwargsTemplate():
	"""Class represents a creator (Factory?) whereby its attributes are fed into the function create(). Any of these attributes can be overwritten temporarily by passing the relevant attribute into the optional arguments of create. When inheriting, createFromSelf() is what needs overwriting, rather than create

	"""

	registeredKwargs = set() #Neccesary to figure out valid Kwargs

	def __init__(self, **kwargs):
		""" General initializer (Template class docstring) for making a "creator" object. Kwargs can be found in self.registeredKwargs (can query the class itself). For each valid keyword setattr(self,keyword,value) will be used (i.e. these keywords set the attributes of the same name)
		
		Raises:
			KeyError: If an unregistered kwarg is passed
		"""
		#First initialise all arguments
		for key in self.registeredKwargs:
			setattr(self,key,None)

		#Now set all arguments where keyword is valid
		for key in kwargs:
			if key in self.registeredKwargs:
				setattr(self,key,kwargs[key])
			else:
				raise KeyError("{} is an invalid keyword.\n Available kwargs are {}".format(key , self.registeredKwargs))


	def create(self,**kwargs):
		""" Create desired object (this is the template docstring so cant be specific). Kwargs should all be present in instance.registeredKwargs; setting any one will cause the creation to occur as if that attribute was set on the object (but the state of the object will NOT change)
		
		Args:
			kwargs: see instance.registeredKwargs
				
		Returns
			The desired object
		
		Raises:
			KeyError: If an unregistered kwarg is passed
		"""
		with temporarilySetRegisteredAttrs(self,kwargs):
			self._checkArgsValidForConstruction()
			outObj = self._createFromSelf()
		return outObj


	def _createFromSelf(self):
		""" Hook. Create the object taking into account only the attributes currently set on this factory/creator instance. This method ALWAYS needs overwriting (Template class docstring here)
		
		Returns
			outObj: The object your trying to create with this creator class
		
			
		"""
		raise NotImplementedError("")


	def _checkArgsValidForConstruction(self):
		""" Hook- overriding is optional and probably rarely needed. This checks the arguments on self are valid and throws an error if not. This gets called during create when we've temporarily set the attributes tovalues determined by **kwargs
		
		
		Raises:
			ValueError: If any attribute is set to an invalid value
		"""
		pass







#This defined a context manger (with contextManger as x:....) that temporarily sets any registeredKwargs on a class  to the given values
#Use case is whenever you have a function which is determined by attributes on a class. Original example was creating a template class for creating a plot, but allowing some attributes (such as xLabel) to be overwritten at creation time.
@contextlib.contextmanager
def temporarilySetRegisteredAttrs(instance, modKwargDict):
	registeredKwargs = instance.registeredKwargs
	origKwargDict = {k:getattr(instance,k) for k in instance.registeredKwargs}
	for key in modKwargDict:
		if key in registeredKwargs:
			setattr(instance,key,modKwargDict[key])
		else:
			raise KeyError("{} is not a registeredKwargs. Registered Kwargs are {}".format(key,registeredKwargs))

	try:
		yield
	finally:
		for key in origKwargDict:
			setattr(instance, key, origKwargDict[key])




