
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
			if struct in ["bcc","fcc","hcp"]:
				return [70,40,40]
			elif struct in ["hcp_inter_octa_relaxed_constant_pressure_3_3_2"]:
				return [60,35,35]
			else:
				self._reportInvalidStructKey(struct)

		elif varyType == "fftGridConv":
			return 0.8
		else:
			self._reportInvalidVaryType(varyType)
	
	def getDefKPoints(self, structStr):
		structToKPoints = {"hcp":[10,10,6],
						   "bcc":[10,10,10],
						   "fcc":[10,10,10],
		                   "hcp_inter_octa_relaxed_constant_pressure_3_3_2": [8,8,6]}
		
		return structToKPoints[structStr]
	
	def _getVaryAngularPoints(self,structType):
		if structType in ["bcc","fcc","hcp"]:
			return [x for x in range(30, 65, 5)] #60 angular is the limit for rc=8.0 basis
		elif structType in ["hcp_inter_octa_relaxed_constant_pressure_3_3_2"]:
			return [x for x in range(25,45,5)]
		else:
			raise KeyError("structType = {} not recognised".format(structType))

	def _getVaryRadialPoints(self,structType):
		if structType in ["bcc","fcc","hcp"]:
			return [x for x in range(40,90,10)]
		elif structType in ["hcp_inter_octa_relaxed_constant_pressure_3_3_2"]:
			return [x for x in range(40,80,10)]
		else:
			raise KeyError("structType = {} not recognised".format(structType))

	def _getFFTVaryPoints(self,structType):
		return [0.7,0.6,0.5] #TODO: Set proper values
		
	def _getVaryKPoints(self,structType):
		singleAtomVals = [ [x,x,x] for x in [2,6,8,10,14,20,24,32] ]
		hcpVals = [ [2,2,1],   [6,6,4],	[8,8,5]   ,   [10,10,6],
					[14,14,9], [20,20,12], [24,24,15],   [32,32,20]  ]

		hcp332_vals = [ [2,2,1], [4,4,3], [6,6,4], [8,8,6], [10,10,6] ]
		structToVaryPoints = { "hcp":hcpVals,
							   "bcc":singleAtomVals,
							   "fcc":singleAtomVals,
		                       "hcp_inter_octa_relaxed_constant_pressure_3_3_2": hcp332_vals}
	
		return structToVaryPoints[structType]
	
	def _reportInvalidVaryType(self,varyType):
		raise ValueError("{} is an invalid option for varyType".format(varyType))

	def _reportInvalidStructKey(self,structKey):
		raise ValueError("{} is an invalid option for structKey".format(structKey))

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
			if struct in ["bcc","fcc","hcp"]:
				return [70,40,40]
			elif struct in ["hcp_inter_octa_relaxed_constant_pressure_3_3_2"]:
				return [60,35,35]
			else:
				self._reportInvalidStructKey(struct)
		elif varyType == "fftGridConv":
			return 0.8
		else:
			self._reportInvalidVaryType(varyType)
	
	def getDefKPoints(self, structStr):
		structToKPoints = {"hcp":[10,10,6],
						   "bcc":[10,10,10],
						   "fcc":[10,10,10],
		                   "hcp_inter_octa_relaxed_constant_pressure_3_3_2": [8,8,6]}

		return structToKPoints[structStr]
	
	def _getVaryAngularPoints(self,structType):
		if structType in ["bcc","fcc","hcp"]:
			return [x for x in range(30, 65, 5)] #60 angular is the limit for rc=8.0 basis
		elif structType in ["hcp_inter_octa_relaxed_constant_pressure_3_3_2"]:
			return [x for x in range(25,45,5)]
		else:
			raise KeyError("structType = {} not recognised".format(structType))


	def _getVaryRadialPoints(self,structType):
		if structType in ["bcc","fcc","hcp"]:
			return [x for x in range(40,90,10)]
		elif structType in ["hcp_inter_octa_relaxed_constant_pressure_3_3_2"]:
			return [x for x in range(40,80,10)]
		else:
			raise KeyError("structType = {} not recognised".format(structType))


	def _getFFTVaryPoints(self,structType):
		return [0.7,0.6,0.5] #TODO: Set proper values
		
	def _getVaryKPoints(self,structType):
		singleAtomVals = [ [x,x,x] for x in [2,6,8,10,14,20,24,32] ]
		hcpVals = [ [2,2,1],   [6,6,4],	[8,8,5]   ,   [10,10,6],
					[14,14,9], [20,20,12], [24,24,15],   [32,32,20]  ]

		hcp332_vals = [ [2,2,1], [4,4,3], [6,6,4], [8,8,6], [10,10,6] ]
		structToVaryPoints = { "hcp":hcpVals,
							   "bcc":singleAtomVals,
							   "fcc":singleAtomVals,
		                       "hcp_inter_octa_relaxed_constant_pressure_3_3_2": hcp332_vals}
		return structToVaryPoints[structType]

	def _reportInvalidVaryType(self,varyType):
		raise ValueError("{} is an invalid option for varyType".format(varyType))

	def _reportInvalidStructKey(self,structKey):
		raise ValueError("{} is an invalid option for structKey".format(structKey))

