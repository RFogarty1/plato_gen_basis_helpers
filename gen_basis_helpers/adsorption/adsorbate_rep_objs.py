
from ..shared import base_surface as baseSurface
from ..shared import simple_vector_maths as vectHelp
import itertools as it

""" Objects used to represent adsorbate layers """

class SurfaceWithAdsorbatesStandard(baseSurface.BaseSurface):

	def __init__(self, cleanSurface, adsorbateLayers, fixAddedVac=True):
		""" Initializer
		
		Args:
			cleanSurface: (BaseSurface) Represents the main surface (e.g. a metal)
			adsorbateLayers: (iter of AdsorbateLayer objects) Each represents a layer of adsorbates; likely will usually only be one
			fixAddedVac: (Optional, bool) If True(Default) the the unit cell be stay the same size as for cleanSurface, which means theasbolute vacuum length between images will (almost always) decrease. False will lead to lenAbsoluteVacuum being fixed once implemented 
 
		"""
		self.cleanSurface = cleanSurface
		self.adsorbateLayers = list(adsorbateLayers)
		self.fixAddedVac = fixAddedVac
		assert self.fixAddedVac is True, "fixAddedVac needs to be True for now (Not implemented the False version)"

	#Test by making simple AdsorbateLayer mocks; they only really need the .cartCoords
	@property
	def unitCell(self):
		outCell = self.cleanSurface.unitCell
		outCoords = outCell.cartCoords
		for x in self.adsorbateLayers:
			currCoords = x.cartCoords
			outCoords += currCoords
		outCell.cartCoords = outCoords
		return outCell

	@property
	def lengthVacuum(self):
		raise NotImplementedError("")

	#Very Tricky to implement, GenericSurface version uses the _singleLayerCell case + lenVac
	@property
	def lenAbsoluteVacuum(self):
		raise NotImplementedError("")



class AdsorbateLayer():
	"""Class holds all info needed to define an adsorption layer

	Attributes:
		sitePositions: (iter of len-3 iters) Each element contains x,y,z co-ords for one site
		adsorbates: (iter of Adsorbate objs/None) Each index maps to sitePositions, values of None mean an unoccupied site
		distances: (iter of floats) Each index maps to adsorbates, specifies the distance of the adsorbate from the site
		surfVector: (len-3 iter) Vector pointing out from the surface. Used to get adsorbate positions for distance > 0
		eqTol: (float, Optional) The maximum distance between two atoms or site positions for them to compare equal
	"""

	def __init__(self, sitePositions, adsorbates, distances, surfVector, eqTol=1e-5):
		self._eqTol = eqTol
		self.sitePositions = list(sitePositions)
		self.adsorbates = list(adsorbates)
		self.distances = list(distances)
		self.surfVector = list(surfVector)

	@property
	def cartCoords(self):
		normdSurfVect = vectHelp.getUnitVectorFromInpVector(self.surfVector)
		outCoords = list()
		for adsA,distA,sitePosA in it.zip_longest(self.adsorbates,self.distances, self.sitePositions):
			if adsA is not None:
				dispVectorDist = [distA*x for x in normdSurfVect]
				dispVectorTotal = [a+b for a,b in it.zip_longest(dispVectorDist,sitePosA)]
				for currCoord in adsA.geom:
					currXyz, currSymbol = currCoord[:3], currCoord[-1]
					outXyz = [a+b for a,b in it.zip_longest(currXyz,dispVectorTotal)]
					outCoords.append(outXyz + [currSymbol])
		return outCoords

	@property
	def cartCoordsNoAtomSymbols(self):
		outCoords = list()
		for x in self.cartCoords:
			outCoords.append(x[:3])
		return outCoords

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)

		#check site positions
		if len(self.sitePositions) != len(other.sitePositions):
			return False

		for posA, posB in it.zip_longest(self.sitePositions, other.sitePositions):
			diffs = [abs(x-y) for x,y in it.zip_longest(posA,posB)]
			if any([x>eqTol for x in diffs]):
				return False

		#check adsorbates are equal
		if len(self.adsorbates) != len(other.adsorbates):
			return False

		for adsA, adsB in it.zip_longest(self.adsorbates,other.adsorbates):
			if adsorbatesSameWithinError(adsA, adsB) is False:
				return False

		#check distances equal
		if len(self.distances) != len(other.distances):
			return False

		for distA, distB in it.zip_longest(self.distances, other.distances):
			if abs(distB-distA) > eqTol:
				return False

		#check surf vectors equal
		if len(self.surfVector) != len(other.surfVector):
			return False
		diffs = [abs(x-y) for x,y in it.zip_longest(self.surfVector,other.surfVector)]
		if any([x > eqTol for x in diffs]):
			return False

		return True

#Probably dont actually inherit from this, just used to define an interface
class Adsorbate():

	@property
	def geom(self):
		""" Returns iter of [x,y,z,Atom] corresponding to the adsorbate geometry. Origin effectively corresponds to the position of the adsorption site
		
		"""
		raise NotImplementedError("")


def adsorbatesSameWithinError(adsA, adsB, errorTol=1e-5):
	""" Test whether two Adsorbate objects are equal, within a tolerance for co-ordinates 
	
	Args:
		adsA: (Adsorbate obj) Mainly contains geometry for an adsorbate
		adsB: (Adsorbate obj)
		errorTol: (float, optional) The maximum error allowed between two co-ordinates for comparing equal
			 
	Returns
		isEqual: (Bool) True if adsA==adsB, else False
 
	"""

	#Firstly check the element identities
	eleValsA = [x[-1] for x in adsA.geom]
	eleValsB = [x[-1] for x in adsB.geom]
	for eleA,eleB in it.zip_longest(eleValsA, eleValsB):
		if eleA!=eleB:
			return False

	#Secondly check the coords
	coordsA = [x[:3] for x in adsA.geom]
	coordsB = [x[:3] for x in adsB.geom]
	for cA, cB in it.zip_longest(coordsA, coordsB):
		diffs = [abs(a-b) for a,b in it.zip_longest(cA,cB)]
		if any([x>errorTol for x in diffs]):
			return False

	return True

