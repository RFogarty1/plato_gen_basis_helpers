
import copy
import itertools as it

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.shared.unit_convs as uConvHelp

from . import get_indices_from_geom_core as getIdxCore
from . import get_neb_lists as nebListHelp
from . import calc_dists as calcDistHelp

from ..shared import cart_coord_utils as cartHelp
from ..shared import plane_equations as planeEqnHelp


def getIndicesGroupedBySurfLayerForFirstNLayers(inpGeom, surfEles, top=True, distTol=1e-1, nLayers=1):
	""" Convenience (slooow implementation) function for getting indices for first n-layers of a surface
	
	Args:
		inpGeom: (plato_pylib UnitCell object)
		surfEles: (iter of str) Contains elements that can be found on the surface
		top: (Bool) If True we look for the top surface (atoms with highest c values), if False we get bottom surface indices (most -ve c values)
		distTol: (float) Distance tolerance for considering an atom to be in a plane. Probably want to increase this A LOT of the time 
		nLayers: (int) Number of surface layers to get the indices for

	Returns
		 outIndices: (len-n iter of iters) Length matches nLayers. First entry is an iter of indices associated with the first layer, second entry is an iter of indices associated with the second layer etc.
 
	"""
	topVal = True if top else False
	botVal = False if top else True

	outIndicesGrouped = list()
	outIndicesChained = list()
	for n in range(1,nLayers+1):
		currFinder = GetSurfaceIndicesFromGeomStandard(surfEles, top=topVal, bottom=botVal, distTol=distTol, nLayers=n)
		currIndices = currFinder.getIndicesFromInpGeom(inpGeom)
		newIndices = [x for x in currIndices if x not in outIndicesChained]
		outIndicesChained.extend(newIndices)
		outIndicesGrouped.append(newIndices)

	return outIndicesGrouped



class GetSurfaceIndicesFromGeomStandard(getIdxCore.GetSpecialIndicesFromInpGeomTemplate):

	def __init__(self, surfEles, top=True, bottom=True, distTol=1e-1, nLayers=1):
		""" Initializer
		
		Args:
			surfEles: (iter of str) The elements that may be present in the surface
			top: (Bool) If true return indices at the top surface
			bottom: (Bool) If true return indices at the bottom surface
			distTol: (float) Distance from outer surface plane an atom can be 
			nLayers: (int) The number of layers to restrict to. e.g. 2 will return indices for the first TWO layers

		"""
		self.surfEles = list(surfEles)
		self.top = top
		self.bottom = bottom
		self.distTol = distTol
		self.nLayers = nLayers

	@property
	def filterFuncts(self):
		eleFilter = getIdxCore.FilterToExcludeElesNotInList(self.surfEles)
		surfFilter = getIdxCore.FilterToOuterSurfaceAtoms(top=self.top, bottom=self.bottom, distTol=self.distTol, nLayers=self.nLayers)
		return [eleFilter, surfFilter]


class GetTopWaterLayerIndices():

	def __init__(self, waterDetector=None, maxHeightLayer=uConvHelp.ANG_TO_BOHR, top=True):
		""" Initializer
		
		Args:
			waterDetector: (GetWaterMoleculeIndicesFromGeomStandard object) Used to get indices related to water
			maxLayerHeight: (float) Difference in OXYGEN atom heights between top/bottom of the layer
			top: (Bool) If True then look for the O at the top of the surface, if False look for the O at the bottom of the surface
				 
		Returns
			outIndices: (iter of len-3 iters) The indices corresponding to the top(or bottom) layer of water molecules
	 
		"""
		self.waterDetector = GetWaterMoleculeIndicesFromGeomStandard() if waterDetector is None else waterDetector
		self.maxHeightLayer = maxHeightLayer
		self.top = top

	def getIndicesFromInpGeom(self, inpGeom):
		waterIndices = self.waterDetector.getIndicesFromInpGeom(inpGeom)
		allFractCoords = inpGeom.fractCoords
		allWaterAtomIndices = [x for x in it.chain(*waterIndices)]
		relFractCoords = [x for idx,x in enumerate(allFractCoords) if idx in allWaterAtomIndices and x[-1].upper()=="O"]

		#TODO: Probably have an "optimisation" option to skip this bit; the purpose of which is mainly to make sure
		#we deal with PBCs properly
		useCell = copy.deepcopy(inpGeom)
		useCell.fractCoords = relFractCoords
		uCellHelp.foldAtomicPositionsIntoCell(useCell)

		cartCoords = useCell.cartCoords

		planeEqn = cartHelp.getPlaneEqnForOuterSurfaceAtoms(useCell, top=self.top)
		distsFromPlane = [ planeEqn.getDistanceOfPointFromPlane(x[:3]) for x in cartCoords]
		outWaterIndices = list()
		for idx,val in enumerate(distsFromPlane):
			if val < self.maxHeightLayer:
				outWaterIndices.append( waterIndices[idx] )

		return outWaterIndices

class GetWaterMoleculeIndicesFromGeomStandard():

	def __init__(self, minOH=0.01, maxOH=1.2*uConvHelp.ANG_TO_BOHR, minAngle=100, maxAngle=110):
		self.minOH = minOH
		self.maxOH = maxOH
		self.minAngle = minAngle
		self.maxAngle = maxAngle

	def getIndicesFromInpGeom(self, inpGeom):
		cartCoords = inpGeom.cartCoords
		allAtomsNebLists = nebListHelp.getNeighbourListsForInpCell_imagesMappedToCentral(inpGeom, self.maxOH)
		oxyIndices = [idx for idx,unused in enumerate(cartCoords) if cartCoords[idx][-1].upper()=="O"]

		outIndices = list()
		for oIdx in oxyIndices:
			currNebList = allAtomsNebLists[oIdx]
			if len(currNebList)==2:
				eleTypes = [cartCoords[x][-1] for x in currNebList]
				if all([x.upper()=="H" for x in eleTypes]):
					oCoord, hCoordA, hCoordB = [ cartCoords[idx][:3] for idx in [oIdx]+currNebList ]
					hDistA = calcDistHelp.calcSingleDistBetweenCoords_minImageConv(inpGeom, oCoord, hCoordA)
					hDistB = calcDistHelp.calcSingleDistBetweenCoords_minImageConv(inpGeom, oCoord, hCoordB)
					angle = calcDistHelp.calcSingleAngleBetweenCoords_minImageConv(inpGeom, hCoordA, oCoord, hCoordB)
					if self._areDistsAndAnglesConsistent([hDistA,hDistB], angle):
						outIndices.append( [oIdx] + currNebList )

		return outIndices

	def _areDistsAndAnglesConsistent(self, dists, angle):
		if any([x<self.minOH for x in dists]):
			return False

		if any([x>self.maxOH for x in dists]):
			return False

		if angle<self.minAngle:
			return False

		if angle>self.maxAngle:
			return False

		return True



#TODO:
#			waterDetector: (Optional, GetWaterMoleculeIndicesFromGeomStandard object) This is used to filter out water molecules before looking for surface atoms, only needed if the surface contains oxygen and/or hydrogen atoms
class GetIndicesForVaryTypesOfSurfAtom_waterBilayerAdsorbedSimple():
	""" Gets indices of surface atoms separated by type when a water bilayer exists above the surface

	Original use is the Mg-water system, where 3 types of surface atom exist; "free" which has no water above it,
	"close" which has an oxygen directly adsorbed (usually the water lies ~flat and parralel) and "far" whereby water
	 sits directly above/below the site but too far away to be "close"

	 Note that these sites are mutually exclusive (an atom can only belong to ONE of them)

	"""

	def __init__(self, surfaceDetector, surfEles, maxWaterPlaneDist, maxCloseDist, maxHozDist, waterOxyLabel="O"):
		""" Initializer
		
		Args:
			surfaceDetector: (GetSurfaceIndicesFromGeomStandard) Gets indices of surface atoms; should be set to either top/bottom surface
			surfEles: (iter of str) Elements that could be in the surface. This CANNOT contain waterOxyLabel
			maxWaterPlaneDist: (float) The maximum distance of water from the surface plane (used to limit to the first bilayer). Note the plane position will be set to the middle of the surface atoms
			maxCloseDist: (float) Maximum Surface-O distance for that surface site to be considered "close"
			maxHozDist: (float) Maximum out-of-plane Surface-O distance for a surface site to be considered as a "far" site
			waterOxyLabel: (str) The label for the oxygen in water
		 
		"""
		self.surfaceDetector = surfaceDetector
		self.surfEles = surfEles
		self.maxWaterPlaneDist = maxWaterPlaneDist
		self.maxCloseDist = maxCloseDist
		self.maxHozDist = maxHozDist
		self.waterOxyLabel = waterOxyLabel


	def getIndicesFromInpGeom(self, inpGeom):
		startCoords = inpGeom.cartCoords
		inpIndices = [x for x in range(len(startCoords))]
		dudGetIndicesInstance = None

		#Step 1 - get indices of the surface atoms
		surfIndices = self.surfaceDetector.getIndicesFromInpGeom(inpGeom)

		#Step 2 - get the average surface plane equation
		planeEqn = cartHelp.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(inpGeom) #Direction shouldnt matter
		dVals = list()
		for idx in surfIndices:
			currXyz = startCoords[idx][:3]
			currDVal = planeEqn.calcDForInpXyz(currXyz)
			dVals.append(currDVal)
		surfDVal = sum(dVals)/len(dVals)
		surfPlaneEqn = planeEqnHelp.ThreeDimPlaneEquation(*(planeEqn.coeffs[:3] + [surfDVal]))

		#Step 3 - get indices of the closest adsorbed bilayer
		currFilterFunct = getIdxCore.FilterToAtomsWithinDistanceOfSurfacePlane(surfPlaneEqn, self.maxWaterPlaneDist)
		bilayerIndices = currFilterFunct(dudGetIndicesInstance, inpGeom, inpIndices)
		bilayerIndices = [x for x in bilayerIndices if x not in surfIndices]
		bilayerIndices = [x for x in bilayerIndices if startCoords[x][-1].upper()==self.waterOxyLabel.upper()]

		#Step 4 - Get the "close" layer based on distance to oxygen atoms
		if self.waterOxyLabel in self.surfEles:
			raise ValueError("{} (waterOxyLabel) found in surfEles; this wont work. Please rename either surface-oxygen atoms or water oxygen atoms".format(self.waterOxyLabel))

		nebPairs = [ [self.waterOxyLabel,x] for x in self.surfEles ]
		closeFilterFunct = getIdxCore.FilterToExcludeIndicesWithoutNebsAmongstRemaning(self.maxCloseDist, restrictToPairs=nebPairs)
		closeSurfIndices = [x for x in closeFilterFunct(dudGetIndicesInstance, inpGeom, surfIndices+bilayerIndices) if x in surfIndices]

		#Step 5 - Get the Hup/Hdown sites based on horizontal distance criteria
		inpCoords = [ startCoords[idx][:3] for idx in bilayerIndices ]
		hozDistFilterFunct = getIdxCore.FilterToExcludeIndicesFurtherOutOfPlaneThanCutoff(self.maxHozDist, surfPlaneEqn, inpCoords)
		farSurfIndices = [idx for idx in hozDistFilterFunct(dudGetIndicesInstance, inpGeom, surfIndices) if idx not in closeSurfIndices]

		#Step 6 - The remainig surface atoms must be "free"
		freeSurfIndices = [idx for idx in surfIndices if (idx not in closeSurfIndices) and (idx not in farSurfIndices)]

		return {"free":freeSurfIndices, "far":farSurfIndices, "close":closeSurfIndices}

