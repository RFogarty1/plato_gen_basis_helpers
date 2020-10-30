

class SurfaceAtomsFromInpGeom():
	""" Callable class for getting a surface object from unitCell. Main job is to filter the adsorbates out from the surface. Interface should be f(inpGeom)->surfAtoms; see self.getSurfaceObjFromInpGeom

	"""
	def getSurfaceAtomsFromInpGeom(self, inpGeom):
		""" Get co-ordinates for surface atoms 
		
		Args:
			inpGeom: (UnitCell object) Contains xyz co-ords
				 
		Returns
			surfCoords: (iter of len-4 iters) Cartesian co-ordinates ([x,y,z,symbol]) for surface atoms
	 
		"""
		raise NotImplementedError("")

	def __call__(self, inpGeom):
		return self.getSurfaceAtomsFromInpGeom(inpGeom)


class AdsorbatesFromInpGeom():
	""" Callable class for getting an iter of adsorbate objects from unitCell. Interface should be f(inpGeom)->adsorbateObjs. see self.getAdsorbateObjsFromInpGeom. 

	"""

	def getAdsorbateObjsFromInpGeom(self, inpGeom):
		""" Get an iter of adsorbate objects
		
		Args:
			inpGeom: (UnitCell object) Contains xyz co-ords
				 
		Returns
			outObjs: (iter of Adsorbate objects) 

		NOTE: This will return adsorbates where ANY atom is in the central cell, thus the number of returned values COULD be larger than the number of adsorbates per cell (and can fluctuate over a closed-system simulation)

		"""
		raise NotImplementedError("")

	def __call__(self, inpGeom, **kwargs):
		return self.getAdsorbateObjsFromInpGeom(inpGeom)


