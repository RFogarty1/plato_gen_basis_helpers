
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

	#
	cartCoords = inpCell.cartCoords
	surfPlaneEqn = cartHelp.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(inpCell)
	adsSitePositions = [x.positionFromGeom(inpCell, inpCartCoords=cartCoords) for x in adsSiteObjs]
	for idx in inpIndices:
		currPos = cartCoords[idx][:3]
		minHozDist, currSiteName = maxHozDist, "None"
		for adsIdx,adsPos in enumerate(adsSitePositions):
			currDist = planeEqnHelp.getOutOfPlaneDistTwoPoints(currPos, adsPos, surfPlaneEqn)
			if currDist<minHozDist:
				#Check the maximum distance criterion is satisfied
				if maxTotDist is None:
					currSiteName = adsSiteObjs[adsIdx].siteName
					minHozDist = currDist
				else:
					totDist = calcDistHelp.calcSingleDistBetweenCoords_minImageConv(inpCell, adsPos, currPos)
					if totDist <= maxTotDist:
						currSiteName = adsSiteObjs[adsIdx].siteName
						minHozDist = currDist

		outDict[currSiteName].append(idx)

	return outDict


