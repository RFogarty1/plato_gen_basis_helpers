
import types
from scipy.optimize import minimize



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
	fitRes = minimize(objectiveFunct, startCoeffs,method=method, **kwargs)
	objectiveFunct(fitRes.x) #Run once more to get the optimised parameters. Should also writeTables as a side-effect
	output = types.SimpleNamespace(optRes=fitRes)
	return output


