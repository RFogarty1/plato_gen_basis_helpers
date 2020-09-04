
from . import base_flow as baseFlow
from ..gau_prod_theorem import get_ints_s_expansions as sIntHelp
from ..fit_cp2k_basis import core as coreFit
from ..fit_cp2k_basis import adapter_stdinp as stdInpAdapterHelp

import types


class BasisFunctSelfOverlapAtDistWorkflow(baseFlow.BaseLabelledWorkflow):
	""" Workflow for finding the overlap of a basis function with itself at a given distance along z

	"""

	def __init__(self, basisObj, angMom, dist):
		""" Initializer
		
		Args:
			basisObj: (plato_pylib GauPolyBasis object) This contains the Gaussian expansion for the basis function
			angMom: (int) The angular momentum for the basis function (s=0, p=1, d=2)
			dist: (float) The distance from origin of the second basis function (the first is at the origin)
				 
		"""
		self.basisObj = basisObj
		self.angMom = angMom
		self.dist = dist
		self._output = types.SimpleNamespace(**{k:None for k in self.namespaceAttrs})

	@property
	def namespaceAttrs(self):
		return ["overlap"]

	@property
	def output(self):
		return [self._output]

	def _checkParamsOkForRun(self):
		if self.angMom != 0:
			raise ValueError("angular momentum of {} currently not supported".format(self.angMom))
		assert self.basisObj.nPoly == 0, "Only r^0 terms allowed (nPoly=0) but found nPoly={}".format(self.basisObj.nPoly)
		

	def run(self):
		self._checkParamsOkForRun()
		overlapVal = sIntHelp.getSelfOverlapMcWedaWeightFromGauPolyBasis(self.basisObj, self.dist)
		self._output.overlap = overlapVal


class BasisFunctSelfOverlapCoeffUpdater(coreFit.CoeffObserver):
	""" Observer object; updates attached BasisFunctSelfOverlapAtDistWorkflow instance whenver it recieves an update on a set of coefficients being used 

	"""

	def __init__(self, workflow, coeffsToPolyBasMapper):
		""" Initializer
		
		Args:
			workflow: (BasisFunctSelfOverlapAtDistWorkflow workflow)
			coeffsToPolyBasMapper: (CoeffsTransformer object) Callable object/function which converts a set of input coefficients to the basis function required
				 
		"""
		self.workflow = workflow
		self.coeffsToPolyBasMapper = coeffsToPolyBasMapper

	def updateCoeffs(self, coeffs):
		newBasisObj = self.coeffsToPolyBasMapper(coeffs)
		self.workflow.basisObj = newBasisObj


class SelfOverlapWorkflowToObjFunctValStandard(stdInpAdapterHelp.WorkflowOutputToObjFunctValStandard):

	def __init__(self, targVals, targAndActValsToObjVal):
		""" Initializer
		
		Args:
			targVals: (iter of float) Target values for the overlap workflow; will usually be length 1
			targAndActValsToObjVal: f(targVals,actVals)->objFunctVal

		"""
		self.targVals = targVals
		self.targAndActValsToObjVal = targAndActValsToObjVal

	#Get the act values in same format as target
	def _getValuesFromOutput(self, output):
		outVals = [x.overlap for x in output]
		return outVals

	#convert a set of target/actual values to an objective function value
	def _getObjFunctValFromTargValsAndActVals(self, targVals, actVals):
		return self.targAndActValsToObjVal(targVals,actVals)


