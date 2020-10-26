
import copy
import gen_basis_helpers.shared.simple_vector_maths as vectHelp

class BaseSurfaceToSites():
	""" Callable class that takes a surface object and returns positions of adsorbate sites. Also has a function that returns a vector pointing out from the surface (which allows us to place adsorbates at varying distances from the adsorption sites

	"""

	def getOutwardsSurfaceVectorFromSurface(self, inpSurface):
		raise NotImplementedError("")

	def getSurfaceSitesFromInpSurface(self, inpSurface):
		raise NotImplementedError("")

	def __call__(self, inpSurface):
		return self.getSurfaceSitesFromInpSurface(inpSurface)

def getWrappedSurfaceToSitesToReorderBasedOnABCentrality(inpFunct):
	""" Takes a function f(inpSurface)->sitePositions and returns a new version whereby sitePositions are reordered such that the sites closest to the centre of the inpSurface are first in the list (closest to [0.5,0.5,x] fract positions, surface should be defined in AB plane for this)
	
	Args:
		inpFunct: f(inpSurface)->sitePositions
			 
	Returns
		outFunct: f(inpSurface)->reOrderedSitePositions, the output is now ordered based on distances of sites from centre of the cell
 
	"""
	def outFunct(inpSurf):
		startPositions = inpFunct(inpSurf)
		cellCentre = _getCentralCartCoordsForCell(inpSurf.unitCell) 
		sortFunct = lambda x: vectHelp.getDistTwoVectors(x,cellCentre)
		return sorted(startPositions,key=sortFunct)	
	return outFunct

def _getCentralCartCoordsForCell(inpCell):
	tempCell = copy.deepcopy(inpCell)
	fractCoords = tempCell.fractCoords
	fractCoords.append([0.5,0.5,0.5,"X"])
	tempCell.fractCoords = fractCoords
	return tempCell.cartCoords[-1][:3]
	


