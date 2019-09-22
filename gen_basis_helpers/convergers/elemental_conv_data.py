
''' Centralised place for parameters to use when carrying out convergence calcs '''

from . import total_energy_conv as eConvergers

def createElementConfigStructs(elementList):
	lowerCaseList = [x.lower() for x in elementList]
	outStructs = dict()
	if "mg" in lowerCaseList:
		outStructs["mg"] = MgElementConfig()
	if "zr" in lowerCaseList:
		outStructs["zr"] = ZrElementConfig()
	return outStructs


# Options for Mg
class MgElementConfig(eConvergers.ElementConfigEqmStructs):
	
	def __init__(self):
		pass
	
	def getConvVariableValues(self,structType, varyType):
		if varyType == "angularGridConv":
			return self._getVaryAngularPoints(structType)
		elif varyType == "fftGridConv":
			return self._getFFTVaryPoints(structType)
		elif varyType == "kPtConv":
			return self._getVaryKPoints(structType)
		elif varyType == "radialGridConv":
			return self._getVaryRadialPoints(structType)
		else:
			self._reportInvalidVaryType(varyType)
	
	def getDefaultMeshSpacing(self,struct,varyType):
		if (varyType=="angularGridConv") or (varyType=="radialGridConv") or (varyType=="kPtConv"):
#			return [5,5,5]
			return [70,40,40]
		elif varyType == "fftGridConv":
			return 0.8
		else:
			self._reportInvalidVaryType(varyType)
	
	def getDefKPoints(self, structStr):
		structToKPoints = {"hcp":[10,10,6],
						   "bcc":[10,10,10],
						   "fcc":[10,10,10]}
		
		return structToKPoints[structStr]
	
	def _getVaryAngularPoints(self,structType):
#		return [5,10,15]
		return [x for x in range(30, 70, 5)]

	def _getVaryRadialPoints(self,structType):
		return [x for x in range(40,90,10)]

	def _getFFTVaryPoints(self,structType):
		return [0.7,0.6,0.5] #TODO: Set proper values
		
	def _getVaryKPoints(self,structType):
		singleAtomVals = [ [x,x,x] for x in [2,6,8,10,14,20,24,32] ]
		hcpVals = [ [2,2,1],   [6,6,4],	[8,8,5]   ,   [10,10,6],
					[14,14,9], [20,20,12], [24,24,15],   [32,32,20]  ]
		structToVaryPoints = { "hcp":hcpVals,
							   "bcc":singleAtomVals,
							   "fcc":singleAtomVals }
	
		return structToVaryPoints[structType]
	
	def _reportInvalidVaryType(self,varyType):
		raise ValueError("{} is an invalid option for varyType".format(varyType))



class ZrElementConfig(eConvergers.ElementConfigEqmStructs):
	
	def __init__(self):
		pass
	
	def getConvVariableValues(self,structType, varyType):
		if varyType == "angularGridConv":
			return self._getVaryAngularPoints(structType)
		elif varyType == "fftGridConv":
			return self._getFFTVaryPoints(structType)
		elif varyType == "kPtConv":
			return self._getVaryKPoints(structType)
		elif varyType == "radialGridConv":
			return self._getVaryRadialPoints(structType)
		else:
			self._reportInvalidVaryType(varyType)
	
	def getDefaultMeshSpacing(self,struct,varyType):
		if (varyType=="angularGridConv") or (varyType=="radialGridConv") or (varyType=="kPtConv"):
#			return [5,5,5]
			return [70,40,40]
		elif varyType == "fftGridConv":
			return 0.8
		else:
			self._reportInvalidVaryType(varyType)
	
	def getDefKPoints(self, structStr):
		structToKPoints = {"hcp":[10,10,6],
						   "bcc":[10,10,10],
						   "fcc":[10,10,10]}
		
		return structToKPoints[structStr]
	
	def _getVaryAngularPoints(self,structType):
#		return [5,8,13]
		return [x for x in range(30, 70, 5)]

	def _getVaryRadialPoints(self,structType):
#		return [10,15,20]
		return [x for x in range(40,90,10)]

	def _getFFTVaryPoints(self,structType):
		return [0.7,0.6,0.5] #TODO: Set proper values
		
	def _getVaryKPoints(self,structType):
#		return [ [1,1,1], [2,2,2] ]
		singleAtomVals = [ [x,x,x] for x in [2,6,8,10,14,20,24,32] ]
		hcpVals = [ [2,2,1],   [6,6,4],	[8,8,5]   ,   [10,10,6],
					[14,14,9], [20,20,12], [24,24,15],   [32,32,20]  ]
		structToVaryPoints = { "hcp":hcpVals,
							   "bcc":singleAtomVals,
							   "fcc":singleAtomVals }
	
		return structToVaryPoints[structType]
	
	def _reportInvalidVaryType(self,varyType):
		raise ValueError("{} is an invalid option for varyType".format(varyType))



