
from . import ref_conv_data_objs as refObjBase

""" NOTE: Just used the Mg values at the moment; and will stay that way until i run the convergence calculations for Zr """

def createConvObj():
	return ZrConvObj()

def createIntegGridConvObj():
	return IntegGridConvZr()

def createKPointGridConvObj():
	return ZrKPointConv()


class ZrConvObj(refObjBase.RefConvergenceDatabase):

	def __init__(self):
		self._integGridVals = ZrIntegGridConv()
		self._kptGridVals = ZrKPointConv()

	@property
	def integGridVals(self):
		return self._integGridVals

	@integGridVals.setter
	def integGridVals(self,value):
		self._integGridVals = value

	@property
	def kptGridVals(self):
		return self._kptGridVals


class ZrIntegGridConv(refObjBase.IntegGridConvergencePureElements):

	def __init__(self):
		pass

	def getPrimCellDft2AngularGrid(self,structKey):
		structToGrid = {"fcc": [50,50,50], "bcc": [60,50,50], "hcp": [50,50,50]}
#		structToGrid = {"fcc": [5,5,5], "bcc": [6,5,5], "hcp": [5,5,5]}
		return structToGrid[ structKey.lower() ]

	def getPrimCellDftFFTGrid(self,structKey):
		structToGrid = {key:0.2 for key in ["fcc","bcc","hcp"]} 
		return structToGrid[structKey]

	def getInterstitialDft2AngularGrid(self,structKey, dims):
		structToGrid = {"fcc": [50,40,40], "bcc": [50,40,40], "hcp": [50,40,40]}
		return structToGrid[structKey]

class ZrKPointConv(refObjBase.KPointConvergence):

	def __init__(self):
		pass

	def getKptsPrimCell(self,structKey):
		structToKpts = {"fcc": [20,20,20], "bcc": [20,20,20], "hcp": [20,20,12]}
#		structToKpts = {"fcc": [1,1,1], "bcc": [1,1,1], "hcp": [1,1,1]}
		return structToKpts[ structKey.lower() ]
		
	def getKPtsSuperCell(self,structKey,dims):
		if tuple(dims) == (1,1,1):
			return self.getKptsPrimCell(structKey)
		structToFunct = {"hcp": _getHcpKPtsForSuperCells}
		return structToFunct[structKey](dims)




def _getHcpKPtsForSuperCells(dims):
	dimKey = tuple(dims)
	dimToKpts = {(4,4,3):[4,4,3],
	             (3,3,2):[5,5,5]}
	return dimToKpts[dimKey]

