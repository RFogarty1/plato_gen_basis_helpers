
import numpy as np

from . import atom_combo_core as atomComboCoreHelp


class _PlanarDistsGetOneDimValsToBin(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, planeEqn, planeDistIndices):
		""" Initializer
		
		Args:
			planeEqn: (ThreeDimPlaneEquation object) The plane equation we want to calculate distances from. Need this to match up to one of the uniquePlaneEquations on the sparseMatrixCalculator we pass to getValsToBin
			planeDistIndices: (iter of ints) The ints of atoms to calculate planar-atom distances
				 
		"""
		self.planeEqn = planeEqn
		self.planeDistIndices = planeDistIndices

	def getValsToBin(self, sparseMatrixCalculator):
		#Get the index of the relevant plane-equation 
		planeEqns = sparseMatrixCalculator.outDict["uniquePlaneEquations"]
		boolVals = [self.planeEqn==x for x in planeEqns]
		assert len([x for x in boolVals if x is True])==1 , "Couldnt find exactly one match for this plane equation"
		planeEqnIdx = boolVals.index(True)

		#
		relPlanarDists = sparseMatrixCalculator.outDict["planarDists"][planeEqnIdx]
		outVals = [ relPlanarDists[idx] for idx in self.planeDistIndices ]
		return outVals

	def __eq__(self, other):
		if type(other) is not type(self):
			return False

		directCmpAttrs = [ "planeEqn", "planeDistIndices" ]
		for attr in directCmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if valA != valB:
				return False
		return True


class _MinDistsGetOneDimValsToBin(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, fromIndices, toIndices):
		""" Initializer
		
		Args:
			fromIndices: (iter of ints) The indices of atoms we calculate distances from. Our output bin values will be len(fromIndices)
			toIndices: (iter of ints) The indices of atoms we calculate distances to
 
		"""
		self.fromIndices = fromIndices
		self.toIndices = toIndices

	def getValsToBin(self, sparseMatrixCalculator):
		relevantMatrix = sparseMatrixCalculator.outDict["distMatrix"]
		outVals = list()
		for idx in self.fromIndices:
			currDists = relevantMatrix[idx][:]
			outVals.append( np.nanmin(currDists[self.toIndices]) )

		return outVals

	def __eq__(self, other):
		if type(other) is not type(self):
			return False

		directCmpAttrs = [ "fromIndices", "toIndices" ]
		for attr in directCmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if valA != valB:
				return False
		return True


