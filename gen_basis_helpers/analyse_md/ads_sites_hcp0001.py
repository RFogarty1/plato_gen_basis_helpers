
import copy

import numpy as np

from . import ads_sites_impl as adsSitesHelp
from . import calc_dists as calcDistHelp
from . import get_neb_lists as nebListHelp


def getHcp0001TopAdsSites(inpGeom, surfIndices, siteNames="top"):
	""" Simple function that gets FixedIndicesAdsSiteBase from a set of surface atom indices. This is sorta trivial, but the function exists mainly so that we have similar interfaces to get all sites
	
	Args:
		inpGeom: (plato_pylib UnitCell object)
		surfInfices: (iter of ints) Indices of atoms you want to consider surface atoms
			 
	Returns
		topSites: (iter of FixedIndicesAdsSiteBase objects) 
 
	"""
	outObjs = [adsSitesHelp.TopStandard(idx,siteName=siteNames) for idx in surfIndices]
	return outObjs

def getHcp0001BridgeAdsSites(inpGeom, surfIndices, distTol=1, maxDist=None, siteNames="bridge", otherSurfIndices=None):
	""" Figures out pairs of atoms that form bridge sites for a hcp0001 surface. Generally each surface atom will have 3 neighbours which it forms a bridge site with
	
	Args:
		inpGeom: (plato_pylib UnitCell object)
		surfIndices: (iter of ints) indices of atoms you want to consider to be surface atoms (passing explicitly can make things easier in cases such as a fully hydroxylated surface)
		siteNames: (str) What to name the name sites
		distTol: (float) If the VERY nearest neighbour distance is x, then neighbours within x+distTol are considered as "nearest neighbours" which can form a bridge site
		maxDist: (float) The maximum distance to consider two atoms to be neighbours. STRONGLY SUGGEST THIS IS SET; at least for efficiency reasons
		otherSurfIndices: (Optional, iter of ints) Default is the same as for surfIndices. Each bridge site is defined by exactly ONE index from surfIndices and ONE index from otherSurfIndices
 
	Returns
		bridgeSites: (iter of FixedIndicesAdsSiteBase objects) 
 
	"""
	#Sort out args
	maxDist = np.inf if maxDist is None else maxDist

	#Sorting should just make everything less messy
	useIndices = sorted(surfIndices)
	otherUseIndices = sorted(otherSurfIndices) if otherSurfIndices is not None else useIndices

	comboUseIndices = sorted( list(set(useIndices + otherUseIndices)) )

	#Sort out the cell we use
	useCell = copy.deepcopy(inpGeom)
	_tempCartCoords = inpGeom.cartCoords
	cartCoords = [_tempCartCoords[idx] for idx in comboUseIndices]
	useCell.cartCoords = cartCoords

	#Get neighbour lists for all atoms + a distance matrix
	nebLists = nebListHelp.getNeighbourListsForInpCell_imagesMappedToCentral(useCell, maxDist)
	distMatrix = calcDistHelp.calcDistanceMatrixForCell_minImageConv(useCell)

	#Remove self from neighbour lists (only needed for cutoff=np.inf i think?)
	for idx, vals in enumerate(nebLists):
		try:
			nebLists[idx].pop( nebLists[idx].index(idx) )
		except ValueError:	
			pass

	numbIndices = len(comboUseIndices)
	unMappedIdxPairs = list()

	for currIdx in range(numbIndices):
		#find any neighbours within r(nearestNeb) + distTol
		currNebList = nebLists[currIdx]

		if len(currNebList) > 0:

			minDist = min( [distMatrix[currIdx][nebIdx] for nebIdx in currNebList] )
			currOutNebs = list()
			for nebIdx in currNebList:
				if (distMatrix[currIdx][nebIdx] <= minDist+distTol):
					currOutNebs.append(nebIdx)

			#Get any pairs where both indices >= currIdx
			currOutPairs = [ [currIdx,nebIdx] for nebIdx in currOutNebs if nebIdx>currIdx]
			unMappedIdxPairs.extend(currOutPairs)


	#Map the indices back to the original cell indices
	mappedIdxPairs = list()
	for idxPair in unMappedIdxPairs:
		mappedIdxPairs.append( [ comboUseIndices[idxPair[0]],comboUseIndices[idxPair[1]] ] )

	#Filter out any indices where two are in ONLY the same list (e.g. both in surfIndices but neither in otherSurfIndices)
	mappedAndFilteredIdxPairs = list()
	for idxPair in mappedIdxPairs:
		idxABools = [idxPair[0] in useIndices, idxPair[0] in otherUseIndices]
		idxBBools = [idxPair[1] in useIndices, idxPair[1] in otherUseIndices]
		comboBools = [ boolA or boolB for boolA,boolB in zip(idxABools, idxBBools)]
		if all(comboBools):
			mappedAndFilteredIdxPairs.append(idxPair)


	#Create the output object from the index pairs
	outObjs = list()
	for idxPair in mappedAndFilteredIdxPairs:
		currObj = adsSitesHelp.BridgeStandard(idxPair, siteName=siteNames)
		outObjs.append(currObj)

	return outObjs


