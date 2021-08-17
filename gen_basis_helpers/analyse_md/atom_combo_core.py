

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





class _SparseMatrixPopulatorComposite(_SparseMatrixPopulator):


	def __init__(self, populators):
		""" Initializer
		
		Args:
			populators: (iter of _SparseMatrixPopulator instances) Note that these should be able to include other _SparseMatrixPopulatorComposite objects (not explicitly tested though...)

		"""
		self.populators = populators

	def populateMatrices(self, inpGeom, outDict, level):
		for populator in self.populators:
			populator.populateMatrices(inpGeom, outDict, level)

	@property
	def maxLevel(self):
		return max([x.maxLevel for x in self.populators])



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




