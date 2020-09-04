
from . import base_flow as baseFlow
from ..gau_prod_theorem import get_ints_s_expansions as sIntHelp

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
