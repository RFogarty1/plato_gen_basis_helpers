

#import plato_pylib.shared.ucell_class as uCellHelp
import itertools as it
import numpy as np
import MDAnalysis.lib.pkdtree as pkdTreeHelp

from . import mdanalysis_interface as mdAnalInter

#TODO: Its possible to use search() instead of search_pairs() to get a more limited set of neighbour lists
#Will probably need to make that optional@ some point
def getNeighbourListsForInpCell_imagesMappedToCentral(inpCell, cutoff):
	""" Gets neighbour lists for all atoms in inpCell. PBCs are taken into account but only the central atom versions are returned (e.g. if an O image is a neighbour, we map that index back to the relevant central cell)
	
	Args:
		inpCell: (UnitCell object)
		cutoff: (float)
 
	Returns
		nebLists: [iter of iters] Each element corresponds to the list of neighbours for the relevant index in inpCell.cartCoords
 
	"""
	#Do data conversions
	boxDims = mdAnalInter.getMDAnalysisDimsFromUCellObj(inpCell)
	coords = np.array( [x[:3] for x in inpCell.cartCoords] )

	#Build the neighbour lists
	treeObj = pkdTreeHelp.PeriodicKDTree(box=boxDims)
	treeObj.set_coords(coords, cutoff=cutoff)
	allPairs = treeObj.search_pairs(cutoff)

	#Get neighbour lists from pairs
	outList = [list() for x in range(len(coords))]
	for idx,pair in enumerate(allPairs):
		iA,iB = pair
		outList[iA].append(iB)
		outList[iB].append(iA)

	#Sort lists and remove duplicates (likely actually unneccesary; no test fails if i dont)
	for idx,unused in enumerate(outList):
		outList[idx] = sorted(list(set(outList[idx])))

	return outList 

