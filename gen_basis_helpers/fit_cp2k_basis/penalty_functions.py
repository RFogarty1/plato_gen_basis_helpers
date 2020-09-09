
import types
from ..workflows import base_flow as baseFlow
from ..shared import calc_runners as calcRunners
from . import core as coreHelp


class CoeffPenaltyFunctionBase():
	"""Class representation of a penalty function. Callable instance with signature __call__(self,coeffs)->penaltyVal where coeffs are the variables being fit 

	"""

	def __call__(self, coeffs):
		raise NotImplementedError("")

def getAdaptedStdInptFromPenaltyFunct(coeffPenaltyFunct):
	""" Get a AdaptedStandardInput object from a CoeffPenaltyFunctionBase object
	
	Args:
		coeffPenaltyFunct: (CoeffPenaltyFunctionBase object)
			 
	Returns
		 outAdaptedInpt: (AdaptedStandardInput object) This is the object used to calcualte partial objective functions. Note the .workflow implemented the CoeffObserver pattern; so thats what needs registering to the coeff updater
 
	"""
	stubLabel = None
	workflowRunFunct = CoeffPenaltyFunctionWorkflow(coeffPenaltyFunct)
	outObj = calcRunners.StandardInputObj(workflowRunFunct, stubLabel)
	return outObj

class CoeffPenaltyFunctionWorkflow(baseFlow.BaseWorkflow, coreHelp.CoeffObserver):

	def __init__(self, penaltyFunct):
		self.penaltyFunct = penaltyFunct
		self.coeffs = None

	def run(self):
		outVal = self.penaltyFunct(self.coeffs)
		self.output = types.SimpleNamespace(objFunct=outVal) #Should be list but seems to need to not be for this

	def updateCoeffs(self, coeffs):
		self.coeffs = list(coeffs)


class MaxAbsValCoeffPenaltyFunctionStandard(CoeffPenaltyFunctionBase):

	def __init__(self, maxValue, penaltyIfOver=None):
		""" Initializer
		
		Args:
			maxValue: (float) The maximum absolute value for the penalty function to return zero
			penaltyIfOver: (function, optional) f(maxVal,actVal) used to get the penalty function when coefficient goes out of bounds. Default is to use the absolute deviation [ (abs(maxVal)-abs(actVal)) ]
				 
		"""
		self.maxValue = maxValue
		defaultPenaltyFunct = lambda targVal,actVal: abs(targVal-actVal)
		self.penaltyIfOver = penaltyIfOver if penaltyIfOver is not None else defaultPenaltyFunct

	def _getValsThatAreOver(self, coeffs):
		outList = list()
		for x in coeffs:
			if abs(x) > abs(self.maxValue):
				outList.append(x)
		return outList

	def __call__(self, coeffs):
		valsOverMax = self._getValsThatAreOver(coeffs)
		outPenalties = [self.penaltyIfOver(self.maxValue,x) for x in valsOverMax]
		outVal = sum(outPenalties)
		return outVal
