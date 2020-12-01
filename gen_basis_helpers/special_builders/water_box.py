
import copy
import itertools as it
import types

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.shared.unit_convs as uConvHelp

from ..adsorption import water_adsorbate as waterAdsHelp
from ..shared import cart_coord_utils as cartHelp

from . import fill_box as boxFillHelp


class GetBulkWaterBoxesFromSurfWithAdsorbatesBothSidesStandard():
	""" Callable class for getting boxes of bulk water relevant to MD for a given surface+adsorbates(both sides) geometry
	"""
	def __init__(self, getEmptyBox, fillEmptyBox):
		""" Initializer
		
		Args:
			getEmptyBox: (GetBulkWaterEmptyBoxFromSurfWithAdsorbatesBothSidesStandard object)
			fillEmptyBox: (GetWaterBoxForMDFromEmptyBoxStandard object)
 
		"""
		self.getEmptyBox = getEmptyBox
		self.fillEmptyBox = fillEmptyBox


	def getBulkWaterBoxes(self, inpCell, nWater, targDensity=None, extraHeight=None, waterAdsGap=None):
		""" Gets boxes relevant to simulating the bulk water region for a surface with water adsorbed on both sides
		
		Args:
			inpCell: (plato_pylib UnitCell object) Cell containing surface geometry with adsorbates on both sides
			nWater: (int) Total number of water molecules required for the simulation

		Kwargs:
			These are all attributes of either self.getEmptyBox or self.fillEmptyBox. Thus, will default to those values if not set. See the relevant docstrings on self.getEmptyBox and self.fillEmptyBox for details
 
		Returns
			outCells: (SimpleNamespace) Contains multiple cells explained below
	 
		outCells:
			cellForMD: (plato_pylib UnitCell) This is the cell that should be used to get starting co-ordinates for the bulk water region
			emptyCell: (plato_pylib UnitCell) This is the cell without any water present AND without a water-Ads gap taken into account. Therefore it may be sligtly larger than cellForMD
			waterAdsGap: (float) This is the gap enforced between the bulk water region and the adsorbate region. Its enforced by making cellForMD smaller than needed based on the density target

		"""
		emptyBox, nWaterBulk = self.getEmptyBox(inpCell, nWater, targDensity=targDensity, extraHeight=extraHeight)
		cellForMD = self.fillEmptyBox(emptyBox, nWaterBulk, waterAdsGap=waterAdsGap)
		outAdsGap = self.fillEmptyBox.waterAdsGap if waterAdsGap is None else waterAdsGap
		outDict = {"emptyCell":emptyBox, "cellForMD":cellForMD, "waterAdsGap":outAdsGap}
		return types.SimpleNamespace(**outDict)

	def __call__(self, inpCell, nWater, targDensity=None, extraHeight=None, waterAdsGap=None):
		return self.getBulkWaterBoxes(inpCell, nWater, targDensity=targDensity, extraHeight=extraHeight, waterAdsGap=waterAdsGap)


class GetWaterBoxForMDFromEmptyBoxStandard():
	""" Callable class for getting a unitCell filled with water from an empty cell and number of water molecules, see getMDCell for interface

	"""

	def __init__(self, primCell, waterAdsObjs, waterAdsGap=0):
		""" Initializer
		
		Args:
			waterAdsGap: (float) Gap between an adsorbate layer and bulk water. This gets subtraced from BOTH ends of the cell (so set to 0.5*excludeDist where excludeDist is the total amount you want the cell shortened by)
			primCell: (plato_pylib UnitCell object) This represents the smallest unit cell of the bulk system with each co-ordinate representing a lattice site (symbol is irrelevant). 
			waterAdsObjs: (list of Adsorbate objects) The length of this needs to match the length of singleHexCell

		"""
		self.primCell = primCell
		self.waterAdsObjs = waterAdsObjs
		self.waterAdsGap = waterAdsGap


	def getMDCell(self, emptyCell, numbWater, waterAdsGap=None):
		cellFiller = boxFillHelp.CellFillerStandard(self.primCell)
		outCell = self._getEmptyCellToUseFromInpCell(emptyCell, waterAdsGap=waterAdsGap)
		lattSites = cellFiller(numbWater, outCell)
		outCoords = self._getCartCoordsForAdsorbatesAtLattSites(lattSites)
		outCell.cartCoords = outCoords
		return outCell

	def _getCartCoordsForAdsorbatesAtLattSites(self, lattSiteCoords):
		adsObjs = it.cycle(self.waterAdsObjs)
		outCoords = list()
		for currSite,adsObj in zip(lattSiteCoords,adsObjs):
			for coord in adsObj.geom:
				currXYZ = [x+t for x,t in zip(coord,currSite[:3])]
				currXYZ.append(coord[-1])
				outCoords.append(currXYZ)
		return outCoords

	def _getEmptyCellToUseFromInpCell(self, inpCell, waterAdsGap=None):
		waterAdsGap = self.waterAdsGap if waterAdsGap is None else waterAdsGap
		outLattParams = [x for x in inpCell.getLattParamsList()]
		outLattParams[-1] -= self.waterAdsGap*2
		outLattAngles = inpCell.getLattAnglesList()
		outCell = uCellHelp.UnitCell(lattParams=outLattParams,lattAngles=outLattAngles)
		return outCell

	def __call__(self, emptyCell, numbWater, waterAdsGap=None):
		return self.getMDCell(emptyCell, numbWater, waterAdsGap=waterAdsGap)



class GetBulkWaterEmptyBoxFromSurfWithAdsorbatesBothSidesStandard():
	""" Callable class (see self.getBoxAndNumberOfWaterToAdd for interface) for getting an empty box for a surface with water on both sides. The box height is set such that adding (N_tot-N_ads) water molecules to it will produce the correct density between the surface layers

		WARNING: This wont work if your adsorbates cross periodic boundaries in the surface plane.

	"""

	def __init__(self,adsDetectorTop, adsDetectorBot, targDensity, extraHeight=0, massDict=None):
		""" Initializer
		
		Args:
			adsDetectorTop: (AdsorbatesFromInpGeom object) Used to get adsorbate objects above the surface from an input UnitCell object
			adsDetectorBot: (AdsorbatesFromInpGeom object) Used to get adsorbate objects below the surface from an input UnitCell object
			targDensity: (float) Target density for the combined bulk and adsorbed water region. This region is h_adsA, h_adsB + h_bulk, where h_adsA/B are interplane spacings between the top/bottom of each adsorbate layer.
			extraHeight: (float, Optional) Value added to the water-region height. Generally passed at runtime(rather than initiation) and used to generate a series of cells with slightly different density values (the stress-tensor for runs on these may then be calculated to find which best represents the NVT ensemble at a given pressure)
			massDict: (dict, Optional) The dictionary used to translate element symbols to masses. Default will likely always be fine though
 
		"""
		self.adsDetectorTop = adsDetectorTop
		self.adsDetectorBot = adsDetectorBot
		self.targDensity = targDensity
		self.extraHeight = extraHeight
		self.massDict = massDict

	#I MAY add another function that returns both the empty box AND the number of water molecules that need adding to it
	def getBox(self, inpCell, nWaterTotal, targDensity=None, extraHeight=None):
		""" Get a UnitCell object for the bulk water region of an MD simulation
		
		Args:
			inpCell: (UnitCell object) Geometry of the surface with adsorbates.
			nWaterTotal: (int) The total number of water molecules you want in your surface+water simulation. This must be at least as many water molecules as found in inpCell
			extraHeight: (float, Optional) Value added to the water-region height. Default is that passed by initialiation of this instance
			targDensity: (float,Optional): Target density in g/units**3 where units are the length units used in inpCell. Default value is determined on initliations; also see the __init__ docstring for a more precise meaning of density
 
		Returns
			emptyBox: (UnitCell object) The cell used to contain the bulk water
	 
		"""
		return self.getBoxAndNumberOfWaterToAdd(inpCell,nWaterTotal, targDensity=targDensity, extraHeight=extraHeight)[0]

	def getBoxAndNumberOfWaterToAdd(self, inpCell, nWaterTotal, targDensity=None, extraHeight=None):
		""" See getBox for details. Only difference is that this function also returns the number of water molecules to add to the box as the second returned parameter
		"""
		targDensity = self.targDensity if targDensity is None else targDensity
		extraHeight = self.extraHeight if extraHeight is None else extraHeight

		#Get all adsorbates above/below
		aboveAdsorbates = self.adsDetectorTop(inpCell)
		belowAdsorbates = self.adsDetectorBot(inpCell)
		heightTopAds = _getHeightOfAdsorbateLayer(aboveAdsorbates, inpCell)
		heightBotAds = _getHeightOfAdsorbateLayer(belowAdsorbates, inpCell)
		numbWaterInCell = len(aboveAdsorbates) + len(belowAdsorbates)
		nBulkWater = nWaterTotal - numbWaterInCell

		#Get the length of the cell for bulk geom
		cellLen = findLatticeParameterToGetTargetDensityForNWater(inpCell, nBulkWater, targDensity, lattParam="c", massDict=self.massDict)
		cellLen -= (heightTopAds + heightBotAds)
		cellLen += extraHeight
		outCell = copy.deepcopy(inpCell)
		lattParams = outCell.getLattParamsList()
		lattParams[-1] = cellLen
		outCell.setLattParams(lattParams)
		outCell.fractCoords = list()

		return outCell,nBulkWater


	def __call__(self, inpCell, nWater, targDensity=None, extraHeight=None):
		return self.getBoxAndNumberOfWaterToAdd(inpCell, nWater, targDensity=targDensity, extraHeight=extraHeight)


#Wont work when they cross PBCs
def _getHeightOfAdsorbateLayer(adsorbateObjs, inpCell):
	#Get all signed distances from a random reference point
	topPlaneEqn = cartHelp.getPlaneEqnForOuterSurfaceAtoms(inpCell)
	outSignedDists = list()
	for adsObj in adsorbateObjs:
		for atomGeom in adsObj.geom:
			currSignedDist = topPlaneEqn.getSignedDistanceOfPointFromPlane(atomGeom[:3])
			outSignedDists.append(currSignedDist)

	return max(outSignedDists)-min(outSignedDists)


def findLatticeParameterToGetTargetDensityForNWater(inpCell, nWater, targDensity, lattParam="c", massDict=None):
	""" Get the lattice parameter for a given UnitCell which will give the targDensity
	
	Args:
		inpCell: (UnitCell object) Contains starting lattice parameters
		nWater: (int) Number of water molecules to put in the cell
		targDensity: (float) Target density. Units are determined by those in inpCell 
		lattParam: (char) The variable lattice parameter, either "a","b" or "c". Function figures out what this needs setting to for getting the targer density assuming other two are fixed
		lenConvFactor: (float) Conversion factor to apply to targDensity. e.g. use BOHR_TO_CM to convert a density fom  
		massDict: (dict) Keys are elements, values are mass per mole. By default its g/mole and ~the values found in any periodic table (i.e. abundance-weighted atomic masses)

	Returns
		 outLattParam: (float) The value the lattice parameter needs setting to in order to get the target density 
 
	"""
	massDict = uCellHelp.getEleKeyToMassDictStandard() if massDict is None else massDict
	massSingleWater = 2*massDict["H"] + massDict["O"]
	totalMass = nWater*massSingleWater*(1/uConvHelp.AVOGADRO_NUMBER)
	surfArea = inpCell.volume*(1/inpCell.lattParams[lattParam])
	return totalMass/(surfArea*targDensity)


def createOneMoleculeHexagonalCellForWater(paramA, paramC, ohDists, angle, azimuth=60, roll=0, pitch=0):
	""" Get a unitCell object with a hexagonal lattice with one water molecule laying flat with Oxygen at the lattice sites. Orientation of water can be controlled with input args, but should generally be left alone especially if you want all water molecules to be in only a single cell.
	
	Args:
		paramA: (float) Lattice parameter a
		paramC: (float) Lattice parameter c
		ohDists: (len 2 iter) O-H bond distances for water molecule
		angle: (float) H-O-H angle for water molecule
		azimuthA: (float) Rotation along azimuth for water molecule A. Default is chosen to minimize odds of an OH bond crossing a cell boundary
 
	Returns
		 outCell: Primitive hexagonal unit cell containing one water molecules
 
	"""
	#Deal with input args
	waterAdsObj = waterAdsHelp.WaterAdsorbateStandard(ohDists, angle, pitch=pitch, roll=roll, azimuthal=azimuth)
	lattSiteFractCoords = [[0.5,0.5,0.5,"X"]]
	return _getHexagonalCellForNWaterMolecules(paramA, paramC, lattSiteFractCoords, [waterAdsObj])


def createTwoMolecularHexagonalCellForWater(paramA, paramC, ohDistsA, angleA, ohDistsB=None, angleB=None, azimuthA=60, azimuthB=-90,
                                            rollA=0, rollB=0, pitchA=0, pitchB=0):
	""" Get a unitCell object with a hexagonal lattice with two water molecules laying flat with Oxygen at the lattice sites. Orientation of water can be controlled with input args, but should generally be left alone especially if you want all water molecules to be in only a single cell.
	
	Args:
		paramA: (float) Lattice parameter a
		paramC: (float) Lattice parameter c
		ohDistsA: (len 2 iter) O-H bond distances for water molecule A
		angleA: (float) H-O-H angle for water molecule A
		ohDistsB: (Optional len-2 iter) O-H bond distances for water molecule B. Defaults to ohDistsA
		angleB: (float) H-O-H angle for water molecule B
		azimuthA: (float) Rotation along azimuth for water molecule A. Default is chosen to minimize odds of an OH bond crossing a cell boundary
		azimuthB: (float) same as A except for water molecule B
 
	Returns
		 outCell: Primitive hexagonal unit cell containing two water molecules
 
	"""

	#Deal with input args
	ohDistsB = ohDistsA if ohDistsB is None else ohDistsB
	angleB = angleA if angleB is None else angleB
	waterAdsA = waterAdsHelp.WaterAdsorbateStandard(ohDistsA, angleA, pitch=pitchA, roll=rollA, azimuthal=azimuthA)
	waterAdsB = waterAdsHelp.WaterAdsorbateStandard(ohDistsB, angleB, pitch=pitchB, roll=rollB, azimuthal=azimuthB)

	lattSiteFractCoords = [[1/3, 1/6, 0.25,"X"], [2/3, 5/6, 0.75,"X"]]

	return _getHexagonalCellForNWaterMolecules(paramA, paramC, lattSiteFractCoords, [waterAdsA,waterAdsB])

def _getHexagonalCellForNWaterMolecules(paramA, paramC, fractCoords, waterAdsObjs):
	#Create the cell
	outCell = uCellHelp.UnitCell(lattParams=[paramA,paramA,paramC], lattAngles=[90,90,120])
	outCell.fractCoords = fractCoords

	#Sort out the co-ordinates
	assert len(fractCoords)==len(waterAdsObjs)
	lattSiteCartCoords = outCell.cartCoords
	outCartCoords = list()
	for lattSite, adsObj in it.zip_longest(lattSiteCartCoords,waterAdsObjs):
		tVector = lattSite[:3]
		currWaterCoords = _getCoordsWithTranslationVectorAdded(adsObj.geom, tVector)
		outCartCoords.extend(currWaterCoords)

	outCell.cartCoords = outCartCoords
	return outCell
	

def _getCoordsWithTranslationVectorAdded(inpCoords, tVect):
	outCoords = list()
	for coord in inpCoords:
		currCoord, currEle = [x for x in coord[:3]], [x for x in coord[3:]]
		translatedCoord = [x+t for x,t in it.zip_longest(currCoord,tVect)]
		outCoords.append(translatedCoord + currEle)
	return outCoords




