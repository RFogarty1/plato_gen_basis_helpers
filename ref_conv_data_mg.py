
import ref_conv_data_objs as refObjBase


def createConvObj():
	return MgConvObj()

def createIntegGridConvObj():
	return IntegGridConvMg()

def createKPointGridConvObj():
	return MgKPointConv()


class MgConvObj(refObjBase.RefConvergenceDatabase):

	def __init__(self):
		self._integGridVals = MgIntegGridConv()
		self._kptGridVals = MgKPointConv()

	@property
	def integGridVals(self):
		return self._integGridVals

	@property
	def kptGridVals(self):
		return self._kptGridVals


class MgIntegGridConv(refObjBase.IntegGridConvergencePureElements):

	def __init__(self):
		pass


	def getPrimCellDft2AngularGrid(self,structKey):
		structToGrid = {"fcc": [50,50,50], "bcc": [60,50,50], "hcp": [50,50,50]}
		return structToGrid[ structKey.lower() ]

	def getPrimCellDftFFTGrid(self,structKey):
		structToGrid = {key:0.15 for key in ["fcc","bcc","hcp"]} 
		return structToGrid[structKey]

	def getInterstitialDft2AngularGrid(self,structKey, dims):
		structToGrid = {"fcc": [50,40,40], "bcc": [50,40,40], "hcp": [50,40,40]}
		return structToGrid[structKey]

class MgKPointConv(refObjBase.KPointConvergence):

	def __init__(self):
		pass

	def getKptsPrimCell(self,structKey):
		structToKpts = {"fcc": [20,20,20], "bcc": [20,20,20], "hcp": [20,20,12]}
		return structToKpts[ structKey.lower() ]
		
	def getKPtsSuperCell(self,structKey,dims):
		if tuple(dims) == (1,1,1):
			return self.getKptsPrimCell(structKey)
		structToFunct = {"hcp": _getHcpKPtsForSuperCells}
		return structToFunct[structKey](dims)




def _getHcpKPtsForSuperCells(dims):
	dimKey = tuple(dims)
	dimToKpts = {(4,4,3):[4,4,3],
	             (3,3,2):[8,8,6]}
	return dimToKpts[dimKey]



