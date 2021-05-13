

class FixedIndicesAdsSiteBase():
	"""The summary line for a class docstring should fit on one line.

	Attributes:
		atomIndices (iter of ints): List of the atom indices used to define the position of this site
		siteName: (str) e.g. "top". Probable use is to insert pseudo-atoms into geometries

	"""

	def positionFromGeom(inpGeom, inpCartCoords=None, **kwargs):
		""" When given an input geometry (with same atomic indexing as used to generate this object) this function will return the position of the adsorbate site
		
		Args:
			inpGeom: (plato_pylib UnitCell instance)
			inpCartCoords: (iter of len-4 iters) [x,y,z,pos]. Passing this may speed up the process compared to getting this info from inpGeom
		 
		Returns
			sitePos: (len-3 iter of xyz)
	 
		"""
		raise NotImplementedError("")


	def __eq__(self, other):
		if sorted(self.atomIndices) != sorted(other.atomIndices):
			return False

		if self.siteName != other.siteName:
			return False

		return True



