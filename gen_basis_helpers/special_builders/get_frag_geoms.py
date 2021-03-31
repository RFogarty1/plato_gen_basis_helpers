

import copy


def getIsolatedFragGeomsFromFragmentIndices(inpGeom, fragIndices):
	""" Returns iter of geometries created from fragments of inpGeom
	
	Args:
		inpGeom: (plato_pylib UnitCell instance)
		fragIndices: (iter of iter of ints) List of indices to use in each output geometry
			 
	Returns
		outGeoms: (iter of UnitCell instance) All have the same cell, but different contents, as inpGeom. The atoms placed in each are determined by fragIndices
 
	"""

	cartCoords = inpGeom.cartCoords
	
	outGeoms = list()
	
	#get isolated fragments
	for frag in fragIndices:
		currCell = copy.deepcopy(inpGeom)
		currCell.cartCoords = [x for idx,x in enumerate(currCell.cartCoords) if idx in frag]
		outGeoms.append(currCell)
	
	
	return outGeoms
	

def getFragPairGeomsFromIndices(inpGeom, fragIndices, maxSep=None):
	""" Description of function
	
	Args:
		inpGeom: (plato_pylib UnitCell object)
		fragIndices: (iter of iter of ints) List of indices to use in each output geometry
		maxSep: (int) Maximum separation for creating pair-geom. For example geoms [A,B,C] with maxSep=1 means AB and BC will be returned, but not AC. The default (None) means return ALL pairs

	Returns
		outGeoms: (iter of UnitCell instance) All have the same cell, but different contents, as inpGeom. Order for simple example [A,B,C] is [AB,AC,BC]
 
	"""
	outGeoms = list()
	for idxA,fragA in enumerate(fragIndices):
		for idxB,fragB in enumerate(fragIndices[idxA+1:], start=idxA+1):
			inclThisFrag = True

			if (maxSep is not None):
				if (idxB-idxA) > maxSep:
					inclThisFrag = False
			if inclThisFrag:
				combinedIndices = fragA + fragB
				currCell = copy.deepcopy(inpGeom)
				currCell.cartCoords = [x for idx,x in enumerate(currCell.cartCoords) if idx in combinedIndices]
				outGeoms.append(currCell)

	return outGeoms


