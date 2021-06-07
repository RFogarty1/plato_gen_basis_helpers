
import copy

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
#	adsSitePositions = [x.positionFromGeom(useCell, inpCartCoords=cartCoords) for x in adsSiteObjs]
	for mappedInpIdx,inpIdx in enumerate(inpIndices):
		currPos = cartCoords[inpIdx][:3]
		minHozDist, currSiteName = maxHozDist, "None"
		for adsIdx,adsPos in enumerate(adsSiteCoords):
			currHozDist = hozDistMatrix[mappedInpIdx][adsIdx]
			if currHozDist < minHozDist:
				#Check maximum distance criterion is satisfied
				if maxTotDist is None:
					currSiteName = adsSiteObjs[adsIdx].siteName
					minHozDist = currHozDist
				else:
					totDist = distMatrix[mappedInpIdx][adsIdx]
					if totDist <= maxTotDist:
						currSiteName = adsSiteObjs[adsIdx].siteName
						minHozDist = currHozDist

		outDict[currSiteName].append(inpIdx)

	return outDict


