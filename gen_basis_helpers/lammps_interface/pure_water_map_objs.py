

class GetPureWaterBondsFromInpGeomCentralCellOnly():

	def __init__(self, adsDetector, bondTypeIdx=1):
		""" Initializer 
		
		Args:
			adsDetector: (DetectH2OAdsorbatesFromInpGeomStandard object). This has a .getCentralCentralAdsorbateIndices function
				 
		"""
		self.adsDetector = adsDetector
		self.bondTypeIdx = bondTypeIdx

	def getBondsFromInpGeom(self, inpGeom):
		cartCoords = inpGeom.cartCoords
		waterIndices = _getWaterIndicesFromInpGeomAndAdsDetector(inpGeom, self.adsDetector)
		outIndices = list()
		currIdx = 1
		for waterMolIndices in waterIndices:
			currIndices = _getWaterBondsFromCartCoordsAndWaterMoleculeIndices(cartCoords, waterMolIndices)
			for indices in currIndices:
				outIndices.append( [currIdx, self.bondTypeIdx] + indices )  
				currIdx += 1

		return outIndices

	def __call__(self, inpGeom):
		return self.getBondsFromInpGeom(inpGeom)


class GetPureWaterAnglesFromInpGeomCentralCellOnly():

	def __init__(self, adsDetector, angleTypeIdx=1):
		""" Initializer 
		
		Args:
			adsDetector: (DetectH2OAdsorbatesFromInpGeomStandard object). This has a .getCentralCentralAdsorbateIndices function
				 
		"""
		self.adsDetector = adsDetector
		self.angleTypeIdx = angleTypeIdx

	def getAnglesFromInpGeom(self, inpGeom):
		cartCoords = inpGeom.cartCoords
		waterIndices = _getWaterIndicesFromInpGeomAndAdsDetector(inpGeom, self.adsDetector)
		outIndices = list()
		currIdx = 1
		for waterMolIndices in waterIndices:
			currIndices = _getWaterAngleIndicesFromCartCoordsAndWaterMoleculeIndices(cartCoords, waterMolIndices)
			outIndices.append( [currIdx, self.angleTypeIdx] + currIndices )
			currIdx+=1
		return outIndices

	def __call__(self, inpGeom):
		return self.getAnglesFromInpGeom(inpGeom)


class GetPureWaterMoleculeIndicesFromInpGeomCentralCellOnly():

	def __init__(self, adsDetector):
		""" Initializer 
		
		Args:
			adsDetector: (DetectH2OAdsorbatesFromInpGeomStandard object). This has a .getCentralCentralAdsorbateIndices function
				 
		"""
		self.adsDetector = adsDetector

	def getMoleculeIndicesFromGeom(self, inpGeom):
		cartCoords = inpGeom.cartCoords
		waterIndices = _getWaterIndicesFromInpGeomAndAdsDetector(inpGeom, self.adsDetector)

		outIndices = [-1 for x in cartCoords]

		for molIdx,waterAdsIndices in enumerate(waterIndices,start=1):
			for atomIdx in waterAdsIndices:
				outIndices[atomIdx] = molIdx

		return outIndices
	
	def __call__(self, inpGeom):
		return self.getMoleculeIndicesFromGeom(inpGeom)

def _getWaterIndicesFromInpGeomAndAdsDetector(inpGeom, adsDetector):
	cartCoords = inpGeom.cartCoords
	waterIndices = adsDetector.getCentralCentralAdsorbateIndices(inpGeom)
	#Check we dont get more water molecules than oxy atoms
	numbWater = len(waterIndices)
	numbOxygen = len( [x for x in cartCoords if x[-1].upper()=="O"] )
	assert numbWater==numbOxygen, "Number water should equal {} (number of oxygen) but {} found".format(numbOxygen, numbWater)
	return waterIndices

def _getWaterBondsFromCartCoordsAndWaterMoleculeIndices(cartCoords, waterIndices):
	oIndices = [x for x in waterIndices if cartCoords[x][-1].upper()=="O"]
	hIndices = [x for x in waterIndices if cartCoords[x][-1].upper()=="H"]

	assert len(oIndices)==1
	outList = [ [oIndices[0]+1,x+1] for x in hIndices]

	return outList

def _getWaterAngleIndicesFromCartCoordsAndWaterMoleculeIndices(cartCoords, waterIndices):
	oIndices = [x for x in waterIndices if cartCoords[x][-1].upper()=="O"]
	hIndices = [x for x in waterIndices if cartCoords[x][-1].upper()=="H"]
	assert len(oIndices)==1
	assert len(hIndices)==2
	outVals = [hIndices[0],oIndices[0],hIndices[1]]
	return [x+1 for x in outVals]


