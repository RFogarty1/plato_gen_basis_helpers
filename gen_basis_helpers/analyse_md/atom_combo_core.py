

""" Backend code for getting combined atomic distributions of various kinds """

import itertools as it

import numpy as np

from . import calc_dists as calcDistsHelp


#Real core/ abstract go here
class _SparseMatrixCalculatorStandard():

	def __init__(self, populators):
		""" Initializer
		
		Args:
			populators: (_SparseMatrixPopulator) Objects which handle calculating specific matrices
		
		NOTES:
			I havent directly tested the multi-level cases at the time of writing
	 
		"""
		self.populators = populators

	def calcMatricesForGeom(self, inpGeom):
		self.outDict = dict() #Need to reset it between calls
		maxLevel = max([x.maxLevel for x in self.populators])
		for level in range(maxLevel+1):
			for populator in self.populators:
				populator.populateMatrices(inpGeom, self.outDict, level)

	def __eq__(self, other):
		#Probably not ACTUALLY needed; due to the zip_longest below
		if len(self.populators) != len(other.populators):
			return False

		for popA, popB in it.zip_longest(self.populators, other.populators):
			if popA != popB:
				return False

		return True


class _SparseMatrixPopulator():
	""" Generic class to populate a single type of sparse matrix (e.g. distance matrix)

	Attributes:
		maxLevel: (int) The maximum "level" for a matrix. level 0 require only a geometry to calculate, level 1 require geom + matrices calculated at level 0 etc. (e.g. setting something as level 1 may allow us to filter the number of values required, such as when looking at both distance and angular criteria for h-bonding
	

	"""

	def populateMatrices(self, inpGeom, outDict, level):
		""" Function to populate (create if missing) relevant parts of matrices in outDict
		
		Args:
			inpGeom: (plato_pylib UnitCell instance)
			outDict: (dict) Values are sparse matrices which we need to populate
			level: (int) Denotes which level matrices we are populating; level 0 means just the inpGeom is needed, level 1 means matrices calculated at level 0 are also needed
				 
		Returns
			Nothing; modifies outDict entries in place
	 
		"""
		raise NotImplementedError("")


class _GetMultiDimValsToBinFromSparseMatrices():
	""" Purpose is in the function getValsToBin() and gets a bunch of values from matrices. One of these should be initialised per iter of calcOptions objects """

	def __init__(self, oneDimObjs):
		self.oneDimObjs = oneDimObjs


	def getValsToBin(self, sparseMatrixCalc):
		""" Takes a bunch of matrices (e.g. distance matrix, horizontal matrix etc.) and returns an iter of values ready to bin
		
		Args:
			sparseMatrixCalc: (_SparseMatrixCalculator object) This stores all matrices that we calculate (e.g. distance matrices)
				 
		Returns
			valsToBin: (iter of len-n iters) 
	 
		"""
		oneDimBinVals = [ x.getValsToBin(sparseMatrixCalc) for x in self.oneDimObjs ]
		multiBinVals = list()

		for binVals in it.zip_longest(*oneDimBinVals):
			multiBinVals.append( binVals )

		return multiBinVals 

	def __eq__(self, other):
		#Probably not ACTUALLY needed; due to the zip_longest below
		if len(self.oneDimObjs) != len(other.oneDimObjs):
			return False

		for objA, objB in it.zip_longest(self.oneDimObjs, other.oneDimObjs):
			if objA != objB:
				return False

		return True

#This needs subclassing for EACH type of options obj
class _GetOneDimValsToBinFromSparseMatricesBase():
	""" Purpose of this is to get a one-dimensional set of values to bin from the _SparseMatrixCalculator object. One of these should correspond to a single calcOptions object """

	def __init__(self):
		raise NotImplementedError("")

	def getValsToBin(self, sparseMatrixCalculator):
		""" Takes a set of basic matrices (e.g. dist matrices) and returns a 1-dimensional set of values to bin
		
		Args:
			sparseMatrixCalc: (_SparseMatrixCalculator object) This stores all matrices that we calculate (e.g. distance matrices)
				 
		Returns
			valsToBin: (iter of floats) The values to bin
	 
		"""
		raise NotImplementedError("")

	@classmethod
	def fromOptsObj(cls, optsObj):
		raise NotImplementedError("")


#Below are largely implementations of _SparseMatrixPopulator and 
class _DistMatrixPopulator(_SparseMatrixPopulator):
	""" Populator meant for calculating distances between sets of indices """

	def __init__(self, fromIndices, toIndices, level=0):
		""" Initializer
		
		Args: (to/from is probably an arbitrary distinction here)
			fromIndices: (iter of ints) Indices of atoms to calculate distances from
			toIndices: (iter of ints) Indices of atoms to calculate distances to
				 
		"""
		self.fromIndices = fromIndices
		self.toIndices = toIndices
		self.level = 0

	@property
	def maxLevel(self):
		return self.level

	def populateMatrices(self, inpGeom, outDict, level):
		if level!=self.level:
			pass
		else:
			self._populateMatrices(inpGeom, outDict)

	def _populateMatrices(self, inpGeom, outDict):
		try:
			useMatrix = outDict["distMatrix"]
		except KeyError:
			currKwargs = {"indicesA":self.fromIndices, "indicesB":self.toIndices, "sparseMatrix":True}
			outDict["distMatrix"] = calcDistsHelp.calcDistanceMatrixForCell_minImageConv(inpGeom, **currKwargs)
		else:
			self._populatePartiallyPopulatedMatrix(inpGeom, outDict)

	def _populatePartiallyPopulatedMatrix(self, inpGeom, outDict):
		useMatrix = outDict["distMatrix"]

		#Figure out indices we need to calculate
		nanIndices = np.argwhere( np.isnan(useMatrix) )
		relIndicesForward = set( (idxA,idxB) for idxA,idxB in it.product(self.fromIndices,self.toIndices) )

		allNanIndices = set([tuple(x) for x in nanIndices.tolist()])
		forwardIndicesIntersection = relIndicesForward.intersection(allNanIndices)
		forwardFromIndices   = sorted(list( set((x[0] for x in forwardIndicesIntersection)) )) 
		forwardToIndices = sorted(list( set((x[1] for x in forwardIndicesIntersection)) ))


		#Calculate only for the relevant indices
		currKwargs = {"indicesA":forwardFromIndices, "indicesB":forwardToIndices, "sparseMatrix":True}
		newSparseMatrix = calcDistsHelp.calcDistanceMatrixForCell_minImageConv(inpGeom, **currKwargs)

		#Populate the relevant parts of the original matrix 
		outDict["distMatrix"] = np.where( np.isnan(useMatrix), newSparseMatrix, useMatrix)


	def __eq__(self, other):
		if type(self) is not type(other):
			return False

		directCmpAttrs = ["fromIndices", "toIndices", "level"]
		for attr in directCmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if valA != valB:
				return False

		return True


class _PlanarDistMatrixPopulator(_SparseMatrixPopulator):
	""" Populator meant for calculating distances between a plane and a set of indices """

	def __init__(self, indices, planeEqn, level=0):
		""" Intializer
		
		Args:
			indices: (iter of ints)
			planeEqn: (ThreeDimPlaneEquation) The plane equation to calculate distance distribution from

		"""
		self.indices = indices
		self.planeEqn = planeEqn
		self.level = level

	@property
	def maxLevel(self):
		return self.level

	def populateMatrices(self, inpGeom, outDict, level):
		if level!=self.level:
			pass
		else:
			self._populateMatrices(inpGeom, outDict)

	def _populateMatrices(self, inpGeom, outDict):
		try:
			unused = outDict["uniquePlaneEquations"]
		except KeyError:
			self._populateMatricesWhenNonePresent(inpGeom, outDict)
		else:
			self._populateMatricesWhenSomePresent(inpGeom, outDict)

	def _populateMatricesWhenNonePresent(self, inpGeom, outDict):
		outDict["uniquePlaneEquations"], outDict["planarDists"] = [self.planeEqn], list()

		currKwargs = {"indices":self.indices, "planeEqn":self.planeEqn, "sparseMatrix":True}
		outDict["planarDists"].append( calcDistsHelp.calcDistancesFromSurfPlaneForCell(inpGeom, **currKwargs) )

	def _populateMatricesWhenSomePresent(self, inpGeom, outDict):
		#Find the relevant index based on the plane equation
		planeEqnIdx = -1
		for idx,planeEqn in enumerate(outDict["uniquePlaneEquations"]):
			if planeEqn == self.planeEqn:
				planeEqnIdx = idx

		#Decide whether we need to modify an (already present) matrix or create a new one
		if planeEqnIdx == -1:
			self._populateSingleMatrixWhenNotPresent(inpGeom, outDict)
		else:
			self._populateSingleMatrixWhenPresent(inpGeom, outDict, planeEqnIdx)

	def _populateSingleMatrixWhenNotPresent(self, inpGeom, outDict):
		currKwargs = {"indices":self.indices, "planeEqn":self.planeEqn, "sparseMatrix":True}
		outMatrix = calcDistsHelp.calcDistancesFromSurfPlaneForCell(inpGeom, **currKwargs)
		outDict["uniquePlaneEquations"].append(self.planeEqn)
		outDict["planarDists"].append(outMatrix)

	def _populateSingleMatrixWhenPresent(self, inpGeom, outDict, matrixIdx):
		#Figure out which indices are needed
		useMatrix = outDict["planarDists"][matrixIdx]
		nanIndices = np.argwhere( np.isnan(useMatrix) )
		neededIndices = np.intersect1d( np.array(self.indices), nanIndices )

		#Calculate
		currKwargs = {"indices":neededIndices, "planeEqn":self.planeEqn, "sparseMatrix":True}
		extraTermsMatrix = calcDistsHelp.calcDistancesFromSurfPlaneForCell(inpGeom, **currKwargs)

		#Update
		outMatrix = np.where( np.isnan(useMatrix), extraTermsMatrix, useMatrix)
		outDict["planarDists"][matrixIdx] = outMatrix



	def __eq__(self, other):
		if type(self) is not type(other):
			return False

		directCmpAttrs = ["indices", "planeEqn", "level"]
		for attr in directCmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if valA != valB:
				return False

		return True





