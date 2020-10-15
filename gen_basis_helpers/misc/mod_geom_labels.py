import copy
from ..shared import cart_coord_utils as cartHelp

""" Used to modify geoms with original goal of excluding certain elements for the Grimme-D3 correction """

def modSurfaceAtomLabels(inpGeom, excludeElements=None, modFunct=None):
	""" Modifies (in place) labels for any atoms in inpGeom that are considered to be on a surface
	
	Args:
		inpGeom: (plato_pylib UnitCell object) Contains the geometry
		excludeElements: (Optional, str iter)  Will limit to only look at elements NOT given in this list. Default is to include ALL elements
		modFunct: (function) f(currLabel)->newLabel. Default is to just add "_surface" to whatever the current label is
			 
	"""
	#-1) Deal with input args
	defFunct = lambda x: x+"_surface"
	modFunct = defFunct if modFunct is None else modFunct

	excludeElements = list() if excludeElements is None else excludeElements

	cartCoords = inpGeom.cartCoords
	rawCoords = [x[:3] for x in cartCoords]


	#0) Get the surface planes (ignoring atoms in excludeElements)
	filteredCoords = [x for x in cartCoords if x[-1] not in excludeElements]
	tempCell = copy.deepcopy(inpGeom)
	tempCell.cartCoords = filteredCoords

	#1) Get the surface planes
	topSurfacePlane = cartHelp.getPlaneEqnForOuterSurfaceAtoms(tempCell, top=True)
	bottomSurfacePlane = cartHelp.getPlaneEqnForOuterSurfaceAtoms(tempCell, top=False)

	#2) Get indices within these planes
	indicesTop = cartHelp.getFilteredIndicesForCoordsInInputPlane(rawCoords, topSurfacePlane)
	indicesBottom = cartHelp.getFilteredIndicesForCoordsInInputPlane(rawCoords, bottomSurfacePlane)
	indicesToMod = indicesTop + indicesBottom #duplicates should be harmless

	#3) Modify based on this
	for idx,coords in enumerate(cartCoords):
		if (idx in indicesToMod) and (cartCoords[idx][-1] not in excludeElements):
			cartCoords[idx][-1] = modFunct(cartCoords[idx][-1])

	#4) modify inpGeom
	inpGeom.cartCoords = cartCoords


def getKindIdxForLabelFromCoords(inpCoords, kindLabel):
	""" Get idx of kindLabel in inpCoords (explained in args). Original purpose is for figuring out indices to exclude in D3 correction in CP2K
	
	Args:
		inpCoords: (len-4 iter) Each element is [x,y,z,label]. Generally will come from UnitCell.cartCoords or UnitCell.fractCoords
		kindLabel: (str) Label we're looking for in inpCoords
			 
	Returns
		 outIdx: (int) kindLabel is the outIdx-th label found. E.g. H in ["O","Mg","O", "O","H","Ne","K","H"] would be outIdx=3 since its the 3rd label to appear at all. NOTE: WE'RE USING INDEXING STARTING AT 1 HERE (NOT ZERO)
 
	Raises:
		 KeyError: if kindLabel not found in inpCoords
	"""
	foundLabels = list()
	currIdx = 1
	outIdx = None
	for x in inpCoords:
		currLabel = x[-1]
		if currLabel==kindLabel:
			outIdx = currIdx
			break
		elif x[-1] not in foundLabels:
			foundLabels.append(x[-1])
			currIdx += 1

	if outIdx is None:
		raise KeyError("{} not found in inpCoords".format(kindLabel))

	return outIdx
