
import copy
from . import lammps_geom as lammpsGeomHelp

class PureWaterCellToLammpsGeomMapperCentralOnly():
	""" Callable class (see getLammpsGeom for interface) to get a LAMMPS geometry object from a UnitCell object containing a pure-water system. NOTE: Won't work if any water molecules cross periodic boundaries

	"""
	def __init__(self, waterDetector, eleToTypeIdx, eleToCharge=None, eleToMass=None, convBohrToAng=True, modTiltFactors=None, atomStyle="full"):
		""" Initializer
		
		Args:
			waterDetector: (DetectH2OAdsorbatesFromInpGeomStandard object). This has a .getCentralCentralAdsorbateIndices function. NOTE: This needs setting up such that it works for geometries in ANGSTROMS
			eleToTypeIdx: (Dict) Keys are element strs while values are the type indices for this element
			eleToCharge: (Dict) Keys are element strs while values are the charge to associate with it
			eleToMass: (Dict) Keys are element strs while values are the mass to associate with it. Default value almost always will be fine
			convBohrToAng: (Bool) If True then the output cell will be converted to angstrom (from bohr)
			modTiltFactors: (list) Factors applied to tilt factors for lammps. This is mainly for hexgonal cells, where tilt factors include 0.5, which can become slightly greater than 0.5 (and throw an error in lammps) due to floats imprecision (or similar). Thus recommended to set to [0.999,0.999,0.999] in this situation. Effectively defaults to [1,1,1]. 
			atomStyle:(str) Only "full" currently supported, so just dont change from default
 
		"""
		self.waterDetector = waterDetector
		self.eleToTypeIdx = eleToTypeIdx
		self.eleToCharge = eleToCharge
		self.eleToMass = eleToMass
		self.convBohrToAng = convBohrToAng
		self.modTiltFactors = modTiltFactors
		self.atomStyle = atomStyle

	def getLammpsGeom(self, inpGeom):
		""" Get a LammpsGeom object from a plato_pylib UnitCell
		
		Args:
			inpCell: (plato_pylib UnitCell object) Must only have water present, and water bonds should not cross PBCs (MIGHT be ok if some atoms are outside the cell boundaries, definitely will not work )
				 
		Returns:
			outObj: (LammpsGeom) Lammps geom object for input geometry
	 
		"""
		if self.atomStyle=="full":
			getDataDictFunct = lammpsGeomHelp.GetDataDictFromLammpsGeomAtomStyleFull(modTiltFactors=self.modTiltFactors)
		else:
			raise ValueError("{} is an invalid atomStyle".format(self.atomStyle))

		#Create the map functions
		geomToBondInfo  = GetPureWaterBondsFromInpGeomCentralCellOnly(self.waterDetector)
		geomToAngleInfo = GetPureWaterAnglesFromInpGeomCentralCellOnly(self.waterDetector)
		geomToMoleculeIDs = GetPureWaterMoleculeIndicesFromInpGeomCentralCellOnly(self.waterDetector)

		#Modify the units on inpCell if requested
		outGeom = copy.deepcopy(inpGeom)
		if self.convBohrToAng:
			outGeom.convBohrToAng()

		#Create the output object
		outKwargs = {"eleToTypeIdx":self.eleToTypeIdx, "eleToCharge":self.eleToCharge, "eleToMass":self.eleToMass,
		             "geomToBondInfo":geomToBondInfo, "geomToAngleInfo":geomToAngleInfo, "geomToMoleculeIDs":geomToMoleculeIDs,
		             "getDataDictFunct":getDataDictFunct}
		outObj = lammpsGeomHelp.LammpsGeom(outGeom, **outKwargs)
		return outObj

	def __call__(self, inpGeom):
		return self.getLammpsGeom(inpGeom)


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

	def __eq__(self, other):
		if self.bondTypeIdx != other.bondTypeIdx:
			return False
		if self.adsDetector != other.adsDetector:
			return False
		return True

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

	def __eq__(self, other):
		if self.angleTypeIdx != other.angleTypeIdx:
			return False
		if self.adsDetector != other.adsDetector:
			return False
		return True

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

	def __eq__(self, other):
		if self.adsDetector != other.adsDetector:
			return False
		return True
	
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


