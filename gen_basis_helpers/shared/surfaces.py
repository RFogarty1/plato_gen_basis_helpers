
import copy
import math

import plato_pylib.utils.supercell as supCell

import gen_basis_helpers.shared.base_surface as baseSurface


class Hcp0001Surface(baseSurface.BaseSurface):

	def __init__(self, bulkUCell, nLayers, lenVac):
		""" Initialiser for Hcp001Surface object
		
		Args:
			bulkUCell: A plato_pylib UnitCell object representing the bulk structure. This should be whats needed
			           to reperesent a single layer, so supercells with [n,m,1] dimensions should be used
			nLayers: The number of surface layers to use
			lenvac: The amount of vacuum required between surface images
		
		Raises:
			InvalidSurfaceError(ValueError): If the angles are wrong for the hcp surface 
		"""
		self._bulkCell = bulkUCell
		self._nLayers = nLayers
		self._lenVac = lenVac

		self._checkLattAnglesConsistent()

	def _checkLattAnglesConsistent(self):
		errorTol = 0.1
		expLattAngles = {"alpha":90.0, "beta":90.0, "gamma":120.0}
		actLattAngles = self._bulkCell.lattAngles
		for key in expLattAngles.keys():
			if abs(expLattAngles[key] - actLattAngles[key]) > errorTol:
				raise baseSurface.InvalidSurfaceError("Hcp surfaces need to be built from unit-cells ith angles of 90/90/120, but input has angles of {}".format(self._bulkCell.lattAngles))


	@property
	def unitCell(self):
		startCell = copy.deepcopy( self._bulkCell )
		outBulkCell = supCell.superCellFromUCell(self._bulkCell,[1,1,self._nLayers]) 

		#Need to add the vacuum, half at top and half at bottom
		lattParams = copy.deepcopy( outBulkCell.lattParams )
		lattParams["c"] += self._lenVac
		cartCoords = copy.deepcopy(outBulkCell.cartCoords)
	
		for row in cartCoords:
			row[2] += 0.5*self._lenVac


		outBulkCell.lattParams = lattParams #NOTE: We MUST change lattParams first, since this changes the cart-coords of atoms
		outBulkCell.cartCoords = cartCoords
	
		return outBulkCell

	@property
	def surfaceArea(self):
		methA = self._calcSurfaceAreaUsingVolumeDivByC()
		methB = self._calcSurfaceAreaUsingTrig()
		allowedDiff = 0.01*max([methA,methB])
		assert abs(methB-methA) < allowedDiff, "Varying estimated for surface area, methA={}, methB={}".format(methA,methB)
		return self._calcSurfaceAreaUsingTrig()

	def _calcSurfaceAreaUsingVolumeDivByC(self):
		return (self._bulkCell.volume / self._bulkCell.lattParams["c"])

	def _calcSurfaceAreaUsingTrig(self):
		#latt params a/b define two sides of a triangle, with an angle of 120 between them. Hence cosine rule can be used to get
		#other info needed to get us a surface area
		lattParams = self._bulkCell.lattParams
		baseSqr = (lattParams["a"]**2) + (lattParams["b"]**2) - (2*lattParams["a"]*lattParams["b"]*math.cos(math.radians(120)))
		base = math.sqrt(baseSqr)
		#Use herons formula to get the area of a triangle from 3 sides
		p = (lattParams["a"] + lattParams["b"] + base) / 2
		areaSqr = p*(p-lattParams["a"])*(p-lattParams["b"])*(p-base)
		area = math.sqrt(areaSqr)
		return area*2 #The cell surface is made of two of these trigangles (its a parallelogram in general)



