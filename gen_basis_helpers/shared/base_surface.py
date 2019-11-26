
''' Module contains base objects for representing surfaces and handling their errors '''



class InvalidSurfaceError(ValueError):
	""" Intended to be used when attempting to create a surface from an invalid bulk-cell. e.g. creating a hcp surface
	    from a bulk fcc structure """

class BaseSurface():

	@property
	def surfaceArea(self):
		""" Surface area (of 1 surface) in same units that lattParams/vectors are specified in
		"""
		raise NotImplementedError()


	@property
	def unitCell(self):
		""" plato_pylib UnitCell object representing the unit-cell formed from the surface including the vacuum region. Note: setting properties on this unitCell will not affect the surface object (a new unitCell should be generated each time this property is called)
		"""
		raise NotImplementedError()

	@property
	def lengthVacuum(self):
		""" Length of the vacuum region above/below the surface. Units are the same as lattParams/vectors
		"""
		raise NotImplementedError()
