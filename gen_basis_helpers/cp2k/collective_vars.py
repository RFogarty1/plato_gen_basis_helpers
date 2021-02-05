

class MetaVarStandard():
	""" Class representing a metavariable """

	def __init__(self, index=1, scale=1):
		""" Initializer
		
		Args:
			index: (int) The index for the metavariable
			scale: (float) See CP2K for the scale keyword; it seems somehow related to the width of hills but i havnt really looked at the maths

		"""
		self.index = index
		self.scale = scale

	def addMetaVarToPyCp2kObj(self, pyCp2kObj):
		""" Adds the collective variable to the metadynamics section of the pycp2k obj
		
		Args:
			pyCp2kObj: Backend object used to generate input files
				 
		Returns
			Nothing but adds the metavariable to the pycp2k obj
	 
		Raises:
			Errors
		"""
		metaDynSection = pyCp2kObj.CP2K_INPUT.MOTION.FREE_ENERGY.METADYN
		metaDynSection.METAVAR_add()
		mVar = metaDynSection.METAVAR_list[-1]
		mVar.Colvar = self.index
		mVar.Scale = self.scale

class CollectiveVarStandard():
	""" Class representing a collective variable for CP2K. Main function is self.addColVarToSubsys """

	def addColVarToSubsys(self, pyCp2kObj):
		""" Adds the collective variable to the pycp2k obj
		
		Args:
			pyCp2kObj: Backend object used to generate input files. MUST have force_eval_list sorted for now 
				 
		Returns
			Nothing but adds the collective variable to the pycp2k obj
	 
		"""
		raise NotImplementedError("")



class DistancePointPlaneColVar(CollectiveVarStandard):

	def __init__(self, atomPlaneIndices, atomPointIndex):
		""" Initializer
		
		NOTE:
			Indices start at zero for this object. Mapping to cp2k index system (which starts at one) is handled internally

		Args:
			atomPlaneIndices: (len-3 iter). Indices for atoms which define the plane
			atomPointIndex: (int) Index for atom defining the point

		"""
		self.atomPlaneIndices = atomPlaneIndices
		self.atomPointIndex = atomPointIndex

	def addColVarToSubsys(self, pyCp2kObj):
		pyCp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS.COLVAR_add()
		colVar = pyCp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS.COLVAR_list[-1].DISTANCE_POINT_PLANE
		colVar.Atom_point = self.atomPointIndex + 1
		colVar.Atoms_plane = [x+1 for x in self.atomPlaneIndices]

