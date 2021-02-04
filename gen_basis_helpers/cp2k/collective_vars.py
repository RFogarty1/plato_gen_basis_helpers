

class CollectiveVarStandard():
	""" Class representing a collective variable for CP2K. Main function is self.addColVarToSubsys """

	def addColVarToSubsys(self, pyCp2kObj):
		""" Adds the collective variable to the pycp2k obj
		
		Args:
			pyCp2kObj: Backend object used to generate input files. MUST have force_eval_list sorted for now 
				 
		Returns
			Nothing but adds the collective variable to the pycp2k obj
	 
		Raises:
			Errors
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

