
from . import parse_from_geoms as parseFromGeomsHelp

class DetectSurfaceBasedOnElementsPresent(parseFromGeomsHelp.SurfaceAtomsFromInpGeom):

	def __init__(self, elesToInclude):
		""" Initializer
		
		Args:
			elesToInclude: (iter of case-insensitive strs) e.g. ["Mg","O"] 
				 
		"""
		self.elesToInclude = list( elesToInclude )

	def getSurfaceAtomsFromInpGeom(self, inpGeom):
		cartCoords = inpGeom.cartCoords
		outCoords = list()
		for coord in cartCoords:
			if any( [x.upper()==coord[-1].upper() for x in self.elesToInclude] ):
				outCoords.append(coord)
		return outCoords

