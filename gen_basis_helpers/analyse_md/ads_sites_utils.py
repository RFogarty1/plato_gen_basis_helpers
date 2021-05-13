



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

