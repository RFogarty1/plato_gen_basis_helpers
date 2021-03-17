
import copy
import itertools as it

from ..analyse_md import get_indices_from_geom_impl as getIndicesImplHelp
from ..shared import plane_equations as planeEqnHelp
from ..shared import cart_coord_utils as cartHelp
from ..shared import simple_vector_maths as vectHelp


def getExtraSurfacePlaneOfAtomsByCopyAndTransAdjacentPlane(inpCell, surfEles=None, top=True, distTol=1e-1):
	""" Get co-ordinates required for an extra surface plane. This is accomplished by using the current outer plane as a mirror plane for reflecting the plane adjacent to the current outer plane
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		surfEles: (iter of str) Contains elements that can be found on the surface. Default is to treat all elements present as surface elements
		top: (Bool) If True we look for the top surface (atoms with highest c values), if False we get bottom surface indices (most -ve c values)
		distTol: (float) Distance tolerance for considering an atom to be in a plane. Probably want to increase this A LOT of the time
			 
	Returns
		outCoords: (iter of len-4 iters) [[x1,y1,z2,"chemSymbol]...[xn,yn,zn,"chemSymbol"]] Cartesian co-ordinates for a new layer
 
	"""
	#Step 0: deal with default args
	if surfEles is None:
		surfEles = list(set( [x[-1] for x in inpCell.fractCoords] ) )

	#Step 1 = get indices of the outer plane and its adjacent plane
	outerPlaneIndices, adjPlaneIndices = getIndicesImplHelp.getIndicesGroupedBySurfLayerForFirstNLayers(inpCell, surfEles, top=top, distTol=distTol, nLayers=2)

	#Step 2 = get central plane equation for each set
	surfPlaneEqn = cartHelp.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(inpCell)

	outerPlaneEqn = cartHelp.getAveragePlaneEqnForAtomIndices(inpCell, outerPlaneIndices, surfPlaneEqn)
	adjPlaneEqn = cartHelp.getAveragePlaneEqnForAtomIndices(inpCell, adjPlaneIndices, surfPlaneEqn)

	#Step 3 = get a translation vector to map adjacent plane to new outer plane. This is -2* the vecto mapping the outer and adjacent planes
	outerToAdjVector = planeEqnHelp.getVectorToMoveFromParallelPlanesAToB(outerPlaneEqn, adjPlaneEqn)
	outTransVector = [-2*x for x in outerToAdjVector]

	#Step 4 = get translated adjCoords; this gets us our next layer. The copy.deepcopy should be unnecesary but just being extra defensive
	adjCoords =   [x for idx,x in enumerate(inpCell.cartCoords) if idx in adjPlaneIndices]
	outCoords = list()
	for coord in copy.deepcopy(adjCoords):
		newCoord = [x+t for x,t in it.zip_longest(coord[:3],outTransVector)] + [coord[-1]]
		outCoords.append(newCoord)

	return outCoords




