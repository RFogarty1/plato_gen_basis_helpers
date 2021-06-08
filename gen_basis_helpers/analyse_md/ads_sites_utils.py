
import copy
import itertools as it

import numpy as np

from . import calc_dists as calcDistHelp
from ..shared import cart_coord_utils as cartHelp
from ..shared import plane_equations as planeEqnHelp

class AddAdsSitesToGeomsStandard():
	""" Class for adding adsorbate sites to geometries """

	def __init__(self, adsSitesObjs):
		""" Initializer
		
		Args:
			adsSitesObjs: (iter of FixedIndicesAdsSiteBase objects). Need positionFromGeom function and .siteName attribute
 
		"""
		self.adsSitesObjs = adsSitesObjs

	def addToUnitCell(self, inpCell):
		cartCoords = inpCell.cartCoords
		newCoords = list()
		for siteObj in self.adsSitesObjs:
			currCoords = siteObj.positionFromGeom(inpCell, inpCartCoords=cartCoords)
			currCoords += [siteObj.siteName]
			newCoords.append( currCoords )

		inpCell.cartCoords = cartCoords + newCoords


def getAssignedAdsIndiceForTrajectory(inpTraj, adsSiteObjs, inpIndices, maxHozDist=2, maxTotDist=None):
	""" Runs "assignAdsIndicesToIndividualAdsorptionSites" for each step in inpTraj and returns the output
	
	Args:
		inpTraj: (TrajectoryInMemory) Contains all info for a simulation
		other args: See assignAdsIndicesToIndividualAdsorptionSites for description

	Returns
		assignedAdsIndicesForTraj: (One output per step) Each element is an iter of iters; For each step each element corresponds to one index in adsSiteObjs. The value of each of these elements is a list of indices showing the adsorbates assigned to that site

 
	"""
	sharedArgs = [adsSiteObjs, inpIndices]
	sharedKwargs = {"maxHozDist":maxHozDist, "maxTotDist":maxTotDist}
	outObjs = list()
	for step in inpTraj:
		currCell = step.unitCell
		currData = assignAdsIndicesToIndividualAdsorptionSites(currCell, *sharedArgs, **sharedKwargs)
		outObjs.append(currData)
	return outObjs


def getPlotDataFromAssignedIndices_adsSiteOccupiedOverTime(inpTimes, assignedAdsIndices):
	""" Function to get data to plot showing when adsorbate sites are occupied (xVals are site indices, y-values are times)
	
	Args:
		inpTimes: (iter of floats) The time for each snapshot. Generally get from [x.time for x in traj]
		assignedAdsIndices: (iter of iter of iters) The output of "getAssignedAdsIndiceForTrajectory". Each element is for one step of the trajectory, each element of that contains a list of adsorbate atoms assigned to the adosrbate site of that index. e.g. assignedAdsIndices[stepIdx][adsSiteIdx] contains the atomic indices assigned to adsSiteIdx and time corresponding to stepIdx
			 
	Returns
		plotData: (iter of iter of len-2 iters) Each element contains data to plot for one of the adsorbate sites. In each case this means nx2 data; where 1st column is just the index of the adsorption site (i.e. the same for each data point) while the 2nd column contains either the time of this step (if occupied) or np.nan (if not occupied)
 
	Raises:
		AssertionError: If len(inpTimes)!=len(assignedAdsIndices)

	"""

	assert len(inpTimes)==len(assignedAdsIndices), "len(inpTimes)={}, len(assignedAdsIndices)={}".format( len(inpTimes),len(assignedAdsIndices) )

	outData = [list() for x in range(len(assignedAdsIndices[0]))]
	for time, assignedData in it.zip_longest(inpTimes, assignedAdsIndices):
		for adsSiteIdx,unused in enumerate(assignedData):
			if len(assignedData[adsSiteIdx]) == 0:
				outData[adsSiteIdx].append( [adsSiteIdx, np.nan] )
			elif len(assignedData[adsSiteIdx]) > 0:
				outData[adsSiteIdx].append( [adsSiteIdx,time] )
			else:
				raise ValueError("Shouldnt ever trigger this")

	return outData


def getPlotDataFromAssignedIndices_adsorbatesAdsorbedOverTime(inpTimes, assignedAdsIndices):
	""" Function to get data to plot showing when adsorbates are adsorbed at ANY site (x-values are adsorbate indices, y-values are times)
	
	Args:
		inpTimes: (iter of floats) The time for each snapshot. Generally get from [x.time for x in traj]
		assignedAdsIndices: (iter of iter of iters) The output of "getAssignedAdsIndiceForTrajectory". Each element is for one step of the trajectory, each element of that contains a list of adsorbate atoms assigned to the adosrbate site of that index. e.g. assignedAdsIndices[stepIdx][adsSiteIdx] contains the atomic indices assigned to adsSiteIdx and time corresponding to stepIdx
			 
	Returns
		plotData: (iter of iter of len-2 iters) Each element contains data to plot for one of the adsorbate sites. In each case this means nx2 data; where 1st column is just the index of the adsorbate (i.e. the same for each data point) while the 2nd column contains either the time of this step (if adsorbed) or np.nan (if not adsorbed)

	"""
	mappedData = _getAssignedAdsIndicesInAdsorbingAtomCentricForm(assignedAdsIndices)
	outData = getPlotDataFromAssignedIndices_adsSiteOccupiedOverTime(inpTimes, mappedData)
	return outData

def _getAssignedAdsIndicesInAdsorbingAtomCentricForm(assignedAdsIndices):

	#First figure out ALL the inpIndices that ever adsorb to something
	allIndices = set()
	for assIndices in assignedAdsIndices:
		currAdsorbateIndices = set(it.chain(*assIndices))
		allIndices.update(currAdsorbateIndices)

	#Now map the original indices to the output indices
	outputIndices = [idx for idx in range(len(allIndices))]
	mapInpToOutputIndices = { key:idx for idx,key in enumerate( sorted(list(allIndices)) ) }

	#Now create the output data
	outData = list()
	for stepIdx, assIndices in enumerate(assignedAdsIndices):
		currData = [list() for x in range(len(allIndices))]
		for adsSiteIdx, adsIndices in enumerate( assIndices ):
			mappedAdsIndices = [mapInpToOutputIndices[x] for x in adsIndices] #len-0 or len-1 generally
			for mappedIdx in mappedAdsIndices:
				currData[mappedIdx].append(adsSiteIdx)
		outData.append(currData)
	
	return outData


#Mainly test this by calling from "assignAdsIndicesToAdsorptionSites"
def assignAdsIndicesToIndividualAdsorptionSites(inpCell, adsSiteObjs, inpIndices, maxHozDist=2, maxTotDist=None):
	""" For each adsSite, assign a list of inpIndices which can be assigned to that site (length=0 or 1 generally expected)
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		adsSiteObjs: (iter of FixedIndicesAdsSiteBase objects) Each of these represents an adsorption site
		inpIndices: (iter of ints) Indices of adsorbate atoms in inpCell (we want to assign THESE to adsorption sites)
		maxHozDist: (float) Maximum horizontal (in-plane) distance from adsorption site to adsorabte atom
		maxTotDist: (float, Optional) Maximum total distance from adsorption site to adsorbate atom. Setting to None means maximum distance not taken into account
			 
	Returns
		 siteVsAds: (iter of iters) Each element corresponds to one index in adsSiteObjs. The value of each element is a list of indices showing the adsorbates assigned to that site
 
	"""

	#Get (PBC-aware) total distance and horizontal distance matrices between adsSiteObjs and inpIndices
	inpCellCartCoords = inpCell.cartCoords #want to save accesing this multiple times
	useCell = copy.deepcopy(inpCell)
	nCoords = len(useCell.fractCoords)
	adsSiteCoords = [ site.positionFromGeom(useCell, inpCartCoords=inpCellCartCoords) + ["ads_site"] for site in adsSiteObjs ]
	useCell.cartCoords = useCell.cartCoords + adsSiteCoords
	firstAdsSiteIdx = len(inpCell.cartCoords)
	adsSiteIndices = [x for x in range(firstAdsSiteIdx,firstAdsSiteIdx+len(adsSiteCoords))]

	distMatrix = calcDistHelp.calcDistanceMatrixForCell_minImageConv(useCell, indicesA=inpIndices, indicesB=adsSiteIndices)
	hozDistMatrix = calcDistHelp.calcHozDistMatrixForCell_minImageConv(useCell, indicesA=inpIndices, indicesB=adsSiteIndices)	


	#
	cartCoords = useCell.cartCoords
	surfPlaneEqn = cartHelp.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(useCell)
	outLists = [list() for x in adsSiteObjs]


	for mappedInpIdx,inpIdx in enumerate(inpIndices):
		currPos = cartCoords[inpIdx][:3]
		minHozDist, currSiteIdx = maxHozDist, None
		for adsIdx,adsPos in enumerate(adsSiteCoords):
			currHozDist = hozDistMatrix[mappedInpIdx][adsIdx]
			if currHozDist < minHozDist:
				#Check maximum distance criterion is satisfied
				if maxTotDist is None:
					currSiteIdx = adsIdx
					minHozDist = currHozDist
				else:
					totDist = distMatrix[mappedInpIdx][adsIdx]
					if totDist <= maxTotDist:
						currSiteIdx = adsIdx
						minHozDist = currHozDist

		if currSiteIdx is not None:
			outLists[currSiteIdx].append(inpIdx)

	return outLists


def assignAdsIndicesToAdsorptionSites(inpCell, adsSiteObjs, inpIndices, maxHozDist=2, maxTotDist=None):
	""" For each unique sitePosition label (defined in adsObjs) this function returns a list of atom indices which can be associated with that site
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		adsSiteObjs: (iter of FixedIndicesAdsSiteBase objects) Each of these represents an adsorption site
		inpIndices: (iter of ints) Indices of adsorbate atoms in inpCell (we want to assign THESE to adsorption sites)
		maxHozDist: (float) Maximum horizontal (in-plane) distance from adsorption site to adsorabte atom
		maxTotDist: (float, Optional) Maximum total distance from adsorption site to adsorbate atom. Setting to None means maximum distance not taken into account

	Returns
		outDict: (dict) Keys are sitePositions in "adsObjs" and "None" [for cases where an adsorbate atom isnt close enough to ANY of the sites]. Values are lists of indices associated with the relevant type of site

	NOTES:
		Each index is assigned to at most one site. This is determined by the site with lowest hozDist (rather than totDist)
 
	"""
	#Initialize
	outDict = { key:list() for key in set([x.siteName for x in adsSiteObjs]) } 
	assert "None" not in outDict.keys()
	outDict["None"] = list()

	#Assign vals into indices
	adsorbatesForIdx = assignAdsIndicesToIndividualAdsorptionSites(inpCell, adsSiteObjs, inpIndices, maxHozDist=maxHozDist, maxTotDist=maxTotDist)
	for idx, vals in enumerate(adsorbatesForIdx):
		currSite = adsSiteObjs[idx].siteName
		for val in vals:
			outDict[currSite].append(val)

	#Figure out the ones with None
	assignedIndices = [ x for x in it.chain( *[outDict[key] for key in outDict.keys()] ) ]
	outDict["None"] = [ x for x in inpIndices if x not in assignedIndices ]

	return outDict



