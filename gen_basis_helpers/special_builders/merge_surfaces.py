
import copy
import itertools as it

import plato_pylib.shared.ucell_class as uCellHelp

from ..shared import surfaces as surfHelp
from ..shared import cart_coord_utils as cartHelp

def getGeomForSurfaceMergedWithCells(geomSurface, geomAboveSurface, geomBelowSurface=None, deltaVac=0, centreSurface=True):
	""" Get geometry for surface merged with cells above/below. Original use case is adding bulk water cells to an Mg surface
	
	Args:
		geomSurface: (UnitCell object) Geometry for the surface (e.g. a blank metal) with surface plane defined in a/b
		geomAboveSurface: (UnitCell object) Geometry for cell to add above the surface. a/b must be consistent with geomSurface
		geomBelowSurface:  (UnitCell, Optional) Geometry for cell to add below the surface.
		deltaVac: (Float) The amount of vacuum to add to the system. This value can be thought of as vacuum between geomAboveSurface and either geomBelowSurface or geomSurface
		centreSurface: (Bool, Optional) If True the geomSurface is translated to the centre of the cell. This made it easier to test, but i suspect it could cause problems if vacuum regions are required (vacuum regions likely wouldnt be spotted by certain functions if they wernt at the top/bottom of the cell)

	Returns
		outGeom: (UnitCell object) Geometry for the system with cells merged together
 
	Raises:
		 ValueError: If angles or a/b lattice parameters dont match between geomSurface and geomAboveSurface/geomBelowSurface

	Important notes on methodology:
		a) vacuum in geomSurface is kind of ignored. The cells above/below surface are added to the top/bottom surface plane containing atoms
		b) Atoms in the cells geomAboveSurface and geomBelowSurface are centred, such that the same amount of vacuum lies above/below the top/bottom surface planes. 
		c) The total amount of vacuum can be written v_t = 0.5v_g1 + 0.5v_g2 + deltaVac, where v_g1 and v_g2 are the total heights of vacuum regions in geomAboveSurface and geomBelowSurface respectively
		d) The bottom of the geomAboveSurface is placed above the surface. The TOP of geomBelowSurface is placed directly below the surface; thus its possible you want to turn geomBelowSurface upside-down before adding it

	"""
	_checkLattParamsConsistent(geomSurface, geomAboveSurface, geomBelowSurface)

	#Step 1 = Add enough vacuum such that we should be able to fit cartesian co-ordinates above/below the surface
	outGeom = copy.deepcopy(geomSurface)
	nLayers = 1
	extraHeightAbove = cartHelp.getHeightOfCell_abSurface(geomAboveSurface)
	extraHeightBelow = 0 if geomBelowSurface is None else cartHelp.getHeightOfCell_abSurface(geomBelowSurface)
	extraHeight = extraHeightAbove + extraHeightBelow
	outSurfObj = surfHelp.GenericSurface(outGeom, nLayers, lenVac=extraHeight)
	outGeom = outSurfObj.unitCell

	#Step 2: Get geomAbove/geomBelow with co-ordinates in the centre of the cell (equal vacuum above/below the geometry)
	geomAbove = copy.deepcopy(geomAboveSurface)
	geomBelow = None if geomBelowSurface is None else copy.deepcopy(geomBelowSurface)
	uCellHelp.foldAtomicPositionsIntoCell(geomAbove)
	cartHelp.shiftCoordsToLeaveEqualVacAboveAndBelowSurface(geomAbove)
	if geomBelow is not None:
		uCellHelp.foldAtomicPositionsIntoCell(geomBelow)
		cartHelp.shiftCoordsToLeaveEqualVacAboveAndBelowSurface(geomBelow)

	#Step 3: Calculate the (shifted) co-ordinates of geomAboveSurface
	topSurfPlane = cartHelp.getPlaneEqnForOuterSurfaceAtoms(outGeom)
	topSurfNormal = topSurfPlane.coeffs[:3]
	heightTopSurfPlane = topSurfPlane.getDistanceOfPointFromPlane(outGeom.lattVects[0]) #this latt vector by definition is contained on the bottom surface plane
	transVectToGetToTopSurfPlane = cartHelp.getCVectorCorrespondingToSurfHeightShift(heightTopSurfPlane, topSurfNormal, outGeom.lattVects[-1])
	cartCoordsForAbove = list()
	for coord in geomAbove.cartCoords:
		currCoord = [x+t for x,t in zip(coord[:3],transVectToGetToTopSurfPlane)] + [coord[-1]]
		cartCoordsForAbove.append(currCoord)

	#Step 4: Get co-ordinates for geomBelowSurface
	cartCoordsForBelow = list()
	if geomBelow is not None:
		botSurfPlane = cartHelp.getPlaneEqnForOuterSurfaceAtoms(outGeom, top=False)
		heightBotSurfPlane = botSurfPlane.getDistanceOfPointFromPlane(outGeom.lattVects[0])
		deltaHeightForGeomBelow = heightBotSurfPlane - extraHeightBelow
		transVectToGetToBotSurfPlane = cartHelp.getCVectorCorrespondingToSurfHeightShift(deltaHeightForGeomBelow, topSurfNormal, outGeom.lattVects[-1])
		cartCoordsForBelow = list()
		for coord in geomBelow.cartCoords:
			currCoord = [x+t for x,t in zip(coord[:3],transVectToGetToBotSurfPlane)] + [coord[-1]]
			cartCoordsForBelow.append(currCoord)

	#Step 5: Add all the relevant co-ords together into outGeom
	outCartCoords = outGeom.cartCoords
	for coord in cartCoordsForAbove:
		outCartCoords.append(coord)
	for coord in cartCoordsForBelow:
		outCartCoords.append(coord)
	outGeom.cartCoords = outCartCoords


	#Step 6: Figure out the amount of vacuum to set (and set it)
	nLayerStub = 1
	vacAbove = surfHelp.GenericSurface(geomAbove, nLayerStub ,lenVac=0).lenAbsoluteVacuum
	vacBelow = 0 if geomBelow is None else surfHelp.GenericSurface(geomBelow, nLayerStub ,lenVac=0).lenAbsoluteVacuum
	totalVacLength = 0.5*vacAbove + 0.5*vacBelow + deltaVac
	outSurfObj = surfHelp.GenericSurface(outGeom, nLayerStub ,lenVac=0)
	outSurfObj.lenAbsoluteVacuum = totalVacLength
	outGeom = outSurfObj.unitCell

	#TODO: Figure out if this is a remotely sensible thing to EVER do
	#Step 7: Make sure the original surface coords are in the centre of the final cell. Will probably make this an optional thing (default=True)
	if centreSurface:
		alignCell = copy.deepcopy(outGeom)
		numbSurfAtoms = len(geomSurface.cartCoords)
		alignCell.cartCoords = alignCell.cartCoords[:numbSurfAtoms]	
		topSurfPlane = cartHelp.getPlaneEqnForOuterSurfaceAtoms(alignCell)
		botSurfPlane = cartHelp.getPlaneEqnForOuterSurfaceAtoms(alignCell, top=False)
		pointOnBotSurf = botSurfPlane.getPointClosestToOrigin()
		surfThickness = topSurfPlane.getDistanceOfPointFromPlane(pointOnBotSurf)

		currTopSurfHeight = topSurfPlane.getDistanceOfPointFromPlane( alignCell.lattVects[0] )
		fullCellHeight = cartHelp.getHeightOfCell_abSurface(alignCell)
		targetSurfHeight = (fullCellHeight/2) + (surfThickness/2)
		deltaHeightShift = targetSurfHeight - currTopSurfHeight
		translationVector = cartHelp.getCVectorCorrespondingToSurfHeightShift(deltaHeightShift, topSurfPlane.coeffs[:3], outGeom.lattVects[-1])

		finalCartCoords = list()
		for coord in outGeom.cartCoords:
			currCoords = [x+t for x,t in it.zip_longest(coord[:3], translationVector)] + [coord[-1]]
			finalCartCoords.append(currCoords)

		outGeom.cartCoords = finalCartCoords

	return outGeom


def _checkLattParamsConsistent(geomSurface, geomAboveSurface, geomBelowSurface, angleTol=1e-1, paramTol=1e-2):
	abSurf = geomSurface.getLattParamsList()[:2]
	surfAngles = geomSurface.getLattAnglesList()

	#1) Check the geom above surface
	abAbove = geomAboveSurface.getLattParamsList()[:2]
	anglesAbove = geomAboveSurface.getLattAnglesList()
	
	if not _areListValsWithinErrorTol(abSurf, abAbove, paramTol):
		errMsg = "a/b lattice params need to be equal for geomSurface/geomAbove surface; but values are {} and {}".format(abSurf, abAbove)
		raise ValueError(errMsg)
	
	if not _areListValsWithinErrorTol(surfAngles, anglesAbove, angleTol):
		errMsg = "Lattice angles need to be equal for geomSurface/geomAbove surface; but values are {} and {}".format(surfAngles, anglesAbove)
		raise ValueError(errMsg)


	#2) Check the geom below surface if needed
	if geomBelowSurface is None:
		return None

	abBelow = geomBelowSurface.getLattParamsList()[:2]
	anglesBelow = geomBelowSurface.getLattAnglesList()

	if not _areListValsWithinErrorTol(abSurf, abBelow, paramTol):
		errMsg = "a/b lattice params need to be equal for geomSurface/geomBelow surface; but values are {} and {}".format(abSurf, abBelow)
		raise ValueError(errMsg)

	if not _areListValsWithinErrorTol(surfAngles, anglesBelow, angleTol):
		errMsg = "Lattice angles need to be equal for geomSurface/geomBelow surface; but values are {} and {}".format(surfAngles, anglesBelow)
		raise ValueError(errMsg)


def _areListValsWithinErrorTol(listA, listB, errorTol):
	if len(listA) != len(listB):
		return False

	absDiffs = [abs(x-y) for x,y in zip(listA,listB)]
	if any([x>errorTol for x in absDiffs]):
		return False
	return True

	


