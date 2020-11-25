
import itertools as it
import plato_pylib.shared.ucell_class as uCellHelp

from ..adsorption import water_adsorbate as waterAdsHelp

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


	outCell = uCellHelp.UnitCell(lattParams=[paramA,paramA,paramC], lattAngles=[90,90,120])
	lattSiteFractCoords = [[1/3, 1/6, 0.25,"X"], [2/3, 5/6, 0.75,"X"]]
	outCell.fractCoords = lattSiteFractCoords

	#Sort out the co-ordinates
	lattSiteCartCoords = outCell.cartCoords
	tVectA = lattSiteCartCoords[0][:3]
	tVectB = lattSiteCartCoords[1][:3]
	waterACoords = _getCoordsWithTranslationVectorAdded(waterAdsA.geom, tVectA)
	waterBCoords = _getCoordsWithTranslationVectorAdded(waterAdsB.geom, tVectB)
	outCell.cartCoords = waterACoords + waterBCoords

	return outCell	


	

def _getCoordsWithTranslationVectorAdded(inpCoords, tVect):
	outCoords = list()
	for coord in inpCoords:
		currCoord, currEle = [x for x in coord[:3]], [x for x in coord[3:]]
		translatedCoord = [x+t for x,t in it.zip_longest(currCoord,tVect)]
		outCoords.append(translatedCoord + currEle)
	return outCoords



#TODO
def createWaterAdsorbatesUniformlyDistributedInBox(inpCell, nAds, waterAds, edgeRegionSize=None):
	""" 
	
	Args:
		inpCell: (UnitCell object) Contains lattice parameters and angles for the box
		nAds: (int) Number of water molecules to pack into the box
		waterAds: (Adsorbate object) has a .geom where "O" should generally be at [0,0,0]
		edgeRegionSize:  (iter [a,b,c]) Distance (in same units as lattice params) to keep lattice sites (oxygen atom position) away from 

	Returns
		outCell: (UnitCell object) This will have cartesian coordinates defined such that each water molecule is together (e.g. [01,H1,H2, O2,H3,H4] where O1,H1,H2 are for one water molecule and O2,H3,H4 are for the second)
		outAdsObj: (AdsorbateLayer object)
 
	Raises:
		 Errors
	"""
	pass


