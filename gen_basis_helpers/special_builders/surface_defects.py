
import copy

from ..analyse_md import get_indices_from_geom_impl as getIndicesImplHelp
from ..shared import cart_coord_utils as cartHelp

from . import add_surface_layers as addSurfLayerHelp

def getAdatomGeometryFromInpCleanSurfaceGeomSimpleImpl(inpGeom, adAtomIdx=None, distTol=1, top=True, surfEles=None):
	""" Takes a clean surface geometry and adds an adatom by removing the top surface (except 1 atom) and adding a layer of surface on the bottom of the surface. This has the NET effect of adding an adatom. NOTE: NOT UNIT TESTED AT TIME OF WRITING
	
	Args:
		inpGeom: (plato_pylib UnitCell object)
		distTol: (float): Tolerance for considering whether an atom is in the same plane as another
		top: (Bool) If True we create the adatom on the top of the surface, if False we add it to the bottom of the surface
		surfEles: (iter of str) Can optionally pass to restrict search for surface to these elements (useful if trying to make an adatom for a surface/water system for example)
		adAtomIdx: (int) Optionally can pass the index of the adatom in the input cell
			 
	Returns
		outGeom: (plato_pylib UnitCell object) Geometry with the adatom added
 
	"""
	outGeom = copy.deepcopy(inpGeom)
	topLayerIndices = getIndicesImplHelp.getIndicesGroupedBySurfLayerForFirstNLayers(inpGeom,["Mg"], top=top, distTol=distTol, nLayers=1)[0]
	if adAtomIdx is None:
		adAtomIdx = cartHelp.getMostCentralIdxFromList(inpGeom, topLayerIndices)
	else:
		assert adAtomIdx in topLayerIndices

	outGeom.cartCoords = [x for idx,x in enumerate(outGeom.cartCoords) if (idx not in topLayerIndices) or (idx==adAtomIdx)]
	newSurfLayer = addSurfLayerHelp.getExtraSurfacePlaneOfAtomsByCopyAndTransAdjacentPlane(outGeom, top=not(top), distTol=distTol)
	outGeom.cartCoords = outGeom.cartCoords + newSurfLayer
	return outGeom

