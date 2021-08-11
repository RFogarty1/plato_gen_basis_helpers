
import plato_pylib.shared.unit_convs as uConvHelp


class NudgedBandReplicas():

	#lenConvFactor=1
	def __init__(self, startCoords, endCoords, interCoords=None):
		""" Initializer
		
		Args:
			startCoords: (iter of len-3 iters) The co-ordinates of the initial structure
			endCoords: (iter of len-3 iters) The co-ordinates of the final structure
			interCoords: (iter of iter of len-3 iters) Each element is a set of co-ordinates for an intermediate
 
		NOTE:
			Its fine to add element symbols to co-ords (at the end) but it wont change the output

		"""
		self.startCoords = startCoords
		self.endCoords = endCoords
		self.interCoords = interCoords

	#REALLY this is what needs to work properly; can easily make an alternative class with diff options
	#(e.g. to include ability to specify velocities as well)
	@property
	def optDict(self):
		outCoords = list()
		outCoords.append([x[:3] for x in self.startCoords])

		if self.interCoords is not None:
			for currCoords in self.interCoords:
				outCoords.append( [x[:3] for x in currCoords] )

		outCoords.append([x[:3] for x in self.endCoords])
		return {"nudgedband_replica_coords":outCoords}


class NudgedBandOptsStd():

	def __init__(self, numbReplicas=None, procsPerReplicaEnv=None, springConstant=None, nebType=None, alignFrames=True, rotateFrames=True, printInitConfigInfo=True, optEndPoints=False):
		""" Initializer
		
		Args (Leaving any as None should lead to CP2K defaults being used):
			numbReplicas: (int) Number of replicas (images) to use in the NEB calc
			procsPerReplicaEnv: (int) Number of processors per replica environment. I THINK you need to set this such that floor(nProcs/procesPerReplicaEnv) = numbReplicas or CP2K will calculate extra geometries + be generally super ineficient
			springConstant: (float) Value of the spring constant linking replicas. Units are whatever cp2k uses
			nebType: (str) Type of nudged elastic band calculation to use. IT-NEB (Improved tangent-NEB) or CI-NEB (Climbing image NEB) are the two most likely options
			alignFrames: (Bool) Corresponds to the align_frames keyword. Probably usually safeist to set to False (it can mess with the final geometry)
			rotateFrames: (Bool) Corresponds to the rotate_frames keyword. Probably usually safeist to set to False (it can mess with the final geometry)
			printInitConfigInfo: (Bool) Whether to print a file with the initial configuration information (Program_run_info/INITIAL_CONFIGURATION_INFO keyword in CP2K)
			optEndPoints: (Bool) Whether to optimise the end points as part of the NEB calculation
 
		"""
		self.numbReplicas = numbReplicas
		self.procsPerReplicaEnv = procsPerReplicaEnv
		self.springConstant = springConstant
		self.nebType = nebType
		self.alignFrames = alignFrames
		self.rotateFrames = rotateFrames
		self.printInitConfigInfo = printInitConfigInfo
		self.optEndPoints = optEndPoints

	@property
	def optDict(self):
		attrToKeyMap = {"numbReplicas":"nudgedband_numbReplicas",
		                "procsPerReplicaEnv": "nudgedband_procsPerReplicaEnv",
		                "springConstant": "nudgedBand_springConstant",
		                "nebType": "nudgedband_type",
		                "alignFrames": "nudgedband_alignFrames",
		                "rotateFrames":"nudgedband_rotateFrames",
		                "printInitConfigInfo":"nudgedband_printInitConfigInfo",
		                "optEndPoints":"nudgedband_optEndPoints"}
		outDict = dict()

		for key in attrToKeyMap.keys():
			if getattr(self,key) is not None:
				outDict[ attrToKeyMap[key] ] = getattr(self,key)

		return outDict




