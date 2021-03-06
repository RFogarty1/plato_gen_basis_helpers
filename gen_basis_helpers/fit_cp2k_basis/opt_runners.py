
import types
from scipy.optimize import minimize

from . import core as coreHelp


def decoObjFunctCallToCatchCalledProcessError(objFunctObj, retVal):
	""" Decorates an instance of ObjFunctCalculatorStandard class such that when CalledProcessError is encountered the objective function returns retVal (instead of the program crashing)
	
	Args:
		objFunctObj: (ObjFunctCalculatorStandard instance)
		retVal: (float) Value to return when encountering CalledProcessError
			
	"""
	objFunctObj._doPreRunShellComms = _getWrapPreRunShellCommsToCatchCalledProcessError(objFunctObj._doPreRunShellComms)
	objFunctObj._calcTotalObjFunct = _getWrapCalcObjFunctValsToCatchCalledProcessError(objFunctObj._calcTotalObjFunct, retVal)

def _getWrapPreRunShellCommsToCatchCalledProcessError(inpFunct):
	def outFunct():
		try:
			inpFunct()
		except:
			pass
	return outFunct

def _getWrapCalcObjFunctValsToCatchCalledProcessError(inpFunct, retVal):
	def outFunct():
		try:
			outVal = inpFunct()
		except:
			outVal = retVal
		return outVal
	return outFunct

def carryOutOptimisationBasicOptions(objectiveFunct,startCoeffs, method=None, **kwargs):
	""" Driver to minimize an objective function using scipy.minimize
	
	Args:
		objectiveFunct: (f(coeffs)->val) Function we're trying to minimise
		startCoeffs: (iter) Starting values for the optimisation parameters
		method: (str) Passed to scipy.minimize as method=method
		kwargs: These are all passed to scipy.minimize
	
	Returns
		output: output.fitRes contains the output from the scipy.minimize function
	
	"""

	#Make sure we get the final coefficients in the input format (rather than what minimize sees/works with)
	transformedCoeffObserver = _TransformedCoeffsObserver()
	objectiveFunct.coeffUpdater.addObserver(transformedCoeffObserver)

	#Carry out the fit
	fitRes = minimize(objectiveFunct, startCoeffs,method=method, **kwargs)
	objectiveFunct(fitRes.x) #Run once more to get the optimised parameters. Should also writeTables as a side-effect
	output = types.SimpleNamespace(optRes=fitRes, transformedCoeffs=transformedCoeffObserver.coeffs )
	return output

class _TransformedCoeffsObserver(coreHelp.CoeffObserver):

	def __init__(self):
		pass

	def updateCoeffs(self, coeffs):
		self.coeffs = coeffs
