

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



