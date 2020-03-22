


from . import ref_conv_data_objs as refObjBase


def createConvObj():
	return MgOConvObj()

def createKPointGridConvObj():
	return MgOKPointConv()


class MgOConvObj(refObjBase.RefConvergenceDatabase):

	def __init__(self):
		self._kPtGridVals = MgOKPointConv()

	@property
	def kptGridVals(self):
		return self._kPtGridVals

	@property
	def integGridVals(self):
		raise NotImplementedError("")


class MgOKPointConv(refObjBase.KPointConvergence):

	def __init__(self):
		pass

	def getKptsPrimCell(self,structKey):
		structToKpts = {"rocksalt": [10,10,10]}
		return structToKpts[ structKey.lower() ]
		
	def getKPtsSuperCell(self,structKey,dims):
		raise NotImplementedError("")


