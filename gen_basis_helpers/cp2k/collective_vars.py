
import itertools as it
from ..shared import simple_vector_maths as vectHelp

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


class DistancePointPlaneColVar_att2(CollectiveVarStandard):

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

		#Add the points for the plane
		for idx in self.atomPlaneIndices:
			colVar.POINT_add()
			colVar.POINT_list[-1].Atoms = [idx+1]

		#Add the points for the single atom point
		colVar.POINT_add()
		colVar.POINT_list[-1].Atoms = [self.atomPointIndex + 1]

		#Tell CP2K which points are the plane and which are the single atom point
		colVar.Atom_point = 4
		colVar.Atoms_plane = [1,2,3]

class DistancePointPlaneColVar_fixedPointsForPlane(CollectiveVarStandard):

	def __init__(self, atomPointIndex, planeXyzVals):
		""" Initializer

		NOTE:
			Indices start at zero for this object. Mapping to cp2k index system (which starts at one) is handled internally
	
		Args:
			atomPointIndex: (int) Index for atom defining the point
			planeXyzVals: (iter of len-3 floats) Each is [x,y,z] and together they should define the plane

		"""
		self.atomPointIndex = atomPointIndex
		self.planeXyzVals = planeXyzVals

	def addColVarToSubsys(self, pyCp2kObj):
		pyCp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS.COLVAR_add()
		colVar = pyCp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS.COLVAR_list[-1].DISTANCE_POINT_PLANE

		#Sort out points defining the plane
		for idx in range(len(self.planeXyzVals)):
			colVar.POINT_add()
			colVar.POINT_list[-1].Type = "FIX_POINT"
			colVar.POINT_list[-1].Xyz = self.planeXyzVals[idx]

		#Add the points for the single atom point
		colVar.POINT_add()
		colVar.POINT_list[-1].Atoms = [self.atomPointIndex + 1]

		#Tell CP2K which points are the plane and which are the single atom point
		colVar.Atom_point = 4
		colVar.Atoms_plane = [1,2,3]

def getXyzValsFromSurfacePlaneEquationAndInpCellStandard(planeEqn, inpCell, unitLength=True):
	""" Takes plane Equation (defining a slice through a surface, orthog to a/b) and input cell and returns 3 points in the plane
	
	Args:
		planeEqn: (ThreeDimPlaneEquation object) Normal vector should be orthogonal to a/b cell vectors
		inpCell: (plato_pylib UnitCell object) Used to get a/b vectors; which should span the plane
		unitLength: (Bool) If True then unit translations are used to map generate two fo the three points used

	Returns
		outXyz: (iter of len-3 float iters)
 
	"""

	#Figure out the two translation vectors (which should lie in the plane)
	vectA = inpCell.getLattVects()[0] 
	vectB = inpCell.getLattVects()[1]
	if unitLength:
		vectA = vectHelp.getUnitVectorFromInpVector(vectA)
		vectB = vectHelp.getUnitVectorFromInpVector(vectB)

	pointA = planeEqn.getPointClosestToOrigin()
	pointB = [ a+v1 for a,v1 in it.zip_longest(pointA,vectA) ]
	pointC = [ a+v2 for a,v2 in it.zip_longest(pointA,vectB) ]

	return [pointA, pointB, pointC]

