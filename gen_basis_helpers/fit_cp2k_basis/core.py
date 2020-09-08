
import itertools as it

from ..shared import calc_runners as calcRunners
import plato_pylib.utils.job_running_functs as jobRunHelp

class CoeffsTransformer():
	""" Callable class implements __call__(self, coeffs), which returns the input coefficients in the same or different format. Reasons can be to restrict the values coeffs can take (e.g. by applying normalisation) or to get the full function that the coeffs are being optimised for

	"""
	def __call__(self, coeffs):
		raise NotImplementedError("")


class CoeffUpdaterStandard():
	"""Class Implements observer pattern to alter dependent classes whenever coefficients, which we optimise, change (i.e. every optimisation step). Can also apply a transformation to the coefficients it recieves. Coefficients are passed to this object using the __call__(self, coeffs) interface; they are then pushed to other objects via their updateCoeffs(coeffs) methods.

	"""
	def __init__(self, transformer=None, observers=None):
		""" Initialiser
		
		Args:
			transformer (Optional): (CoeffsTransformer) callable as new_coeffs=transformer(old_coeffs). This allows the coefficients input to the updater to be different to those transmitted; useful (for example) when trying to constrain the actual coefficients used (original pupose is to restrict to normalised version of the actual coefficients)
			observers (Optional): (iter of CoeffObserver) Each object must implement updateCoeffs(coeffs). Every object in this iter is told every time a new set of coefficients is passed to the CoeffUpdater
				 
		"""
		self.observers = list(observers) if observers is not None else list()
		self.transformer = transformer

	def _transformCoeffs(self,coeffs):
		if self.transformer is None:
			outVals = coeffs
		else:
			outVals = self.transformer(coeffs)
		return outVals

	def addObserver(self, observer):
		self.observers.append(observer)

	def __call__(self, coeffs):
		outCoeffs = self._transformCoeffs(coeffs)
		for x in self.observers:
			x.updateCoeffs(outCoeffs)



class CoeffObserver():
	"""Class can be an observer for CoeffUpdater; meaning it must implement an update function (coeff updater uses push-style updating

	"""

	def updateCoeffs(self, coeffs):
		""" Update coefficients on this object. This is called by the updater function/class whenever the relevant coefficients are updated
		
		Args:
			coeffs: (iter) Conains coefficients to update class with
				 
		"""
		raise NotImplementedError("")


class AdaptedStandardInput(calcRunners.StandardInputObj):
	""" Same as StandardInputObj but the createOutputObj() function should return an object(obj) which has property obj.data[0].objFunct

	"""
	pass

class ObjFunctCalculatorStandard():
	"""Objective function object; Callable as __call__(self, coeffs) and returns the objective function for the given set of coefficients

	"""
	def __init__(self, objs, coeffUpdater, nCores=1, weights=None, observers=None):
		""" Initializer
		
		Args:
			objs: (iter of AdaptedStandardInput objects) These contain run comms and each output an objective function
			coeffUpdater: (CoeffUpdaterStandard) Used to communicate the new set of coefficients
			nCores: (int) Number of cores to use for shell comms (generally meaning the running-jobs part)
			weights: (Optional, iter of floats) default is all ones. Multiple the objective values from objs by these weights
			observers: (Optional, list of ObjFunctObserver objects) List of observers. These are updated with the current objective function after every call to this objective function
 
		"""
		self.objs = list(objs)
		self.coeffUpdater = coeffUpdater
		self.nCores = nCores
		self.weights = list(weights) if weights is not None else [1 for x in self.objs]
		self.observers = list(observers) if observers is not None else list()

	def addObjValObserver(self, observer):
		self.observers.append(observer)

	def _updateObservers(self, objVal):
		for x in self.observers:
			x.updateObjVal(objVal)

	def _updateCoeffs(self, coeffs):
		self.coeffUpdater(coeffs)

	def _doPreRunShellComms(self):
		runComms = list()
		for x in self.objs:
			runComms.extend(x.runComms)
		jobRunHelp.executeRunCommsParralel(runComms, self.nCores, quiet=True, noCommsOk=True)

	def _calcTotalObjFunct(self):
		objFunctVals = self._getObjFunctVals()
		return self._combineObjFunctVals(objFunctVals)

	def _getObjFunctVals(self):
		outVals = list()
		for x in self.objs:
			currOutputObj = x.createOutputObj()
			outVals.extend( [x.objFunct for x in currOutputObj.data] )
		return outVals

	def _combineObjFunctVals(self, objFunctVals):
		weightedVals = [x*w for x,w in it.zip_longest(objFunctVals,self.weights)]
		return sum(weightedVals)

	def __call__(self, coeffs):
		self._updateCoeffs(coeffs)
		self._doPreRunShellComms()
		outVal = self._calcTotalObjFunct()
		self._updateObservers(outVal)
		return outVal

class ObjFunctObserver():
	"""Class can act as an observer for ObjFunctCalculatorStandard, since it implements an updateObjVal method

	"""
	def updateObjVal(self, objVal):
		raise NotImplementedError("")




