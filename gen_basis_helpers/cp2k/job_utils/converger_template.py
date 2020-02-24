
import gen_basis_helpers.shared.calc_runners as calcRunners


class BaseConvergerStandardInputCreatorTemplate():
	""" Class should provide a template for creating a Standard input object for running calculations to converge a parameter in the CP2K input file

		Main function should be createStandardInput() [No arguments should be needed]

	"""

	def createStandardInput(self):
		""" Function creates a BaseStandardInputObj for running convergence calculations (and generating an output
		
		Returns
			standardInpt: Instance of a BaseStandardInputObj object. This has all the info required to run the convergence calculations and create an output object for analysis
		
		Raises:
			Errors
		"""
		raise NotImplementedError("")




class StandardConvergerStandardInputTemplate(BaseConvergerStandardInputCreatorTemplate):

	def __init__(self, convObjs, createBasicObjFunct, modCalcObjFunct, attrDict=None):
		""" Initializer for StandardConvergerStandardInputTemplate. Note that inheriting from this class and overwriting methods is another option.
		
		Args:
			convObjs: iter of objects describing the convergence parameters to test (e.g. grid cutoffs)
			createBasicObjFunct: f(self) function. See the class createBasicObject docstring for more info
			modCalcObjFunct: f(self,calcObj,convObj). See the modCalcObjWithConvergenceParams docstring for more info
			attrDict: (dict, optional). This provides a way of setting attributes on the converger, keys are attr names while values are what to set them to. For example you may want to set attrDict={"baseFolder":"some/path"}, in which case createBasicObjFunct and modCalcObjFunct will both be able to access the self.baseFolder attribute
	
		"""
		self._convObjs = convObjs
		self._createBasicObjFunct = createBasicObjFunct
		self._modCalcObjFunct = modCalcObjFunct
		if attrDict is not None:
			for key,val in attrDict.items():
				setattr(self, key, val)


	@property
	def convObjs(self):
		""" An iterable. The interface of the returned type is coupled to the "modCalcObjWithConvergenceParams" function
		"""
		return self._convObjs

	def createStandardInput(self):
		outCalcObjs = list()
		for x in self.convObjs:
			currCalcObj = self.createBasicObject()
			currCalcObj = self.modCalcObjWithConvergenceParams(currCalcObj,x)
			outCalcObjs.append(currCalcObj)
		return calcRunners.StandardInputObjComposite(outCalcObjs)

	def createBasicObject(self):
		""" Creates an object with all shared settings, such as geometry etc. This function is called and the output file modified by 
		
		Returns
			outObj: (CalcMethod Object) This needs to have a way of setting at least the input file paths on top of the standard interface. The implementation of this object is tightly coupled with the supplied modCalcObjWithConvergenceParams function, which needs to be able to make the neccesary modifications to the input method object.
		
		"""

		return self._createBasicObjFunct(self)


	def modCalcObjWithConvergenceParams(self, calcObj, convObj):
		""" Functions takes a basic calculation object and modifies it with the convergence properties defined in convObj. This should also be sorting out the input path to use (since that has to be different for ALL calculations)
		
		Args:
			calcObj: (CalcMethod Object) Contains the basic settings which are shared between all calculations (such as the geometry)
			convObj: Some kind of object that contains the information required to modify calcObj. E.g. it could contain a grid cutoff value (meaning it can be as simple as an integer really)
	
		Returns
			outCalcObj: The modified calcObj. Its fine to modify it in place and then return it, we use the return value to make it more general.
		
		"""
		return self._modCalcObjFunct(self,calcObj,convObj)

