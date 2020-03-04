

""" Purpose of these classes is to provide a standard way of running any calculation sets """


from . import misc_utils as misc

class BaseStandardInputObj():
	""" Base class for representing an "input object"; i.e. one used to manage running calculations. Goal is for this to make it easier to run/manage a load of separate calculations together, such as (for example) calculating a vacancy energy for a range of basis sets.

	"""

	@property
	def label(self):
		""" list of BaseLabel objects associated with input. """
		raise NotImplementedError("")

	@property
	def runComms(self):
		""" iter of run-commands which can be passed to subprocess.call """
		raise NotImplementedError("")

	def createOutputObj(self):
		""" Function to create an output object after runComms have been executed. The output object should implement the BaseStandardOutputObj interface """
		raise NotImplementedError("")


class BaseStandardOutputObj():
	""" Base class for representing a generic output object
	"""

	@property
	def data(self):
		""" An iterable where each entry represents the data for one leaf input object (generally will be one set of calculations of some kind). Format of each entries can vary with implementation (this is the base class doc-string)"""
		raise NotImplementedError("")

	@property
	def label(self):
		""" iter of BaseLabel objects """
		raise NotImplementedError("")


class StandardInputObj(BaseStandardInputObj):

	def __init__(self, workflow, label, mapFunction=None):
		""" Initializer to create a Standard Input Object
		
		Args:
			workflow: (plato_fit_integrals WorkFlowBase object) This workflow defines how to carry out the required calculations and get the neccesary info for the output objects
			label: (StandardLabel object) This contains strings for element, structure and method used, which should be sufficient to differentiate between input objs
			mapFunction: (f(inputObject), optional) Determines processing we carry out for the "data" we obtain in the output object. By default we just return workflow.output (i.e. the workflow itself determines the format). Note the workflow.run() needs to be done within this function. Effectively setting this variable is the same as overriding the createOutputObj function

		"""
		self.workflow = workflow
		self._label = label
		self._mapFunction = mapFunction 

	@property
	def runComms(self):
		return self.workflow.preRunShellComms

	@property
	def label(self):
		return [self._label]

	@property
	def mapFunction(self):
		""" (f(inputObject), optional) Determines processing we carry out for the "data" we obtain in the output object. By default we just return workflow.output (i.e. the workflow itself determines the format). Note the workflow.run() needs to be done within this function. Effectively setting this variable is the same as overriding the createOutputObj function
		"""
		return [self._mapFunction] #Needs to be a list to make it work for composite

	@mapFunction.setter
	def mapFunction(self,val):
		self._mapFunction = val

	def createOutputObj(self):
		#Hook to let user overwrite standard way and add extra processing or similar
		if self._mapFunction is not None:
			output = self._mapFunction(self)
		#Standard way
		else:
			self.workflow.run()
			output = self.workflow.output

		return StandardOutputObj(output, self._label) 


class StandardOutputObj(BaseStandardOutputObj):

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=False)
	def __init__(self, data, label):
		""" Initailizer to create a Standard Output object
		
		Args:
			data: Pretty much anything that represents the output data for one input object
			label: (StandardLabel object) This contains strings for element, structure and method used, which should be sufficient to differentiate between input objs

		"""
		self._data = data
		self._label = label

	@property
	def label(self):
		return [self._label]

	@property
	def data(self):
		return [self._data]


class StandardInputObjComposite(BaseStandardInputObj):

	runComms = misc.StandardComponentDescriptor("runComms")
	label = misc.StandardComponentDescriptor("label")

	def __init__(self, objs):
		self.objs = list(objs)

	def createOutputObj(self):
		outputObjs = list()
		for leaf in self.objs:
			outputObjs.append( leaf.createOutputObj() )
		return StandardOutputObjComposite(outputObjs)

	@property
	def mapFunction(self):
		""" (f(inputObject), optional) Determines processing we carry out for the "data" we obtain in the output object. By default we just return workflow.output (i.e. the workflow itself determines the format). Note the workflow.run() needs to be done within this function. Effectively setting this variable is the same as overriding the createOutputObj function

		Note: If Setting on a composite, you pass one function which is used for ALL leaf objects
		"""

		allMapFuncts = list()
		for branch in self.objs:
			allMapFuncts.extend( branch.mapFunction )
		return allMapFuncts

	@mapFunction.setter
	def mapFunction(self,value):
		for branch in self.objs:
			branch.mapFunction = value


class StandardOutputObjComposite(BaseStandardOutputObj):

	data = misc.StandardComponentDescriptor("data")
	label = misc.StandardComponentDescriptor("label")

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=True)
	def __init__(self, objs):
		self.objs = list(objs)


