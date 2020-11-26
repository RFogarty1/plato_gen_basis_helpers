
import itertools as it

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.shared.unit_convs as uConvHelp

from ..adsorption import water_adsorbate as waterAdsHelp


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
	totalMass = nWater*massSingleWater
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




