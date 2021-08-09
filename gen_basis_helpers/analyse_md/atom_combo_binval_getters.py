
import itertools as it
import numpy as np

from . import atom_combo_core as atomComboCoreHelp
from . import atom_combo_populators as atomComboPopulatorHelp

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

class _WaterPlanarDistBinValGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, oxyIndices, hyIndices, planeEqn, primaryIdxType="O"):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			planeEqn: (ThreeDimPlaneEquation) The plane equation to calculate distance distribution from
			primaryIdxType: (str) The element of the primary index. "O", "Ha" and "Hb" are the standard options
				 
		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.planeEqn = planeEqn
		self.primaryIdxType = primaryIdxType

	def getValsToBin(self, sparseMatrixCalculator):	
		currBinner = self._getBinner()
		return currBinner.getValsToBin(sparseMatrixCalculator)

	def _getBinner(self):
		indicesToBin = self._getIndicesToBin()
		currBinner = _PlanarDistsGetOneDimValsToBin(self.planeEqn, indicesToBin)
		return currBinner

	def _getIndicesToBin(self):
		if self.primaryIdxType.upper() == "O":
			return self.oxyIndices
		elif self.primaryIdxType.upper() == "HA":
			return [x[0] for x in self.hyIndices]
		elif self.primaryIdxType.upper() == "HB":
			return [x[1] for x in self.hyIndices]
		else:
			raise ValueError("")

class _WaterPlanarMinDistBinValGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, oxyIndices, hyIndices, planeEqn, minDistType="all"):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			planeEqn: (ThreeDimPlaneEquation) The plane equation to calculate distance distribution from
			minDistType: (str) Controls which atoms to get the minimum distance from. Current options are "all","o", and "h" (case insensitive)

		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.planeEqn = planeEqn
		self.minDistType = minDistType

	def getValsToBin(self, sparseMatrixCalculator):
		#Get sets of indices
		oxyIndices = self.oxyIndices
		hyIndicesA = [x[0] for x in self.hyIndices]
		hyIndicesB = [x[1] for x in self.hyIndices]

		#Merge into relevant ones
		if self.minDistType.upper() == "ALL":
			relIndices = [oxyIndices, hyIndicesA, hyIndicesB]
		elif self.minDistType.upper() == "O":
			relIndices = [oxyIndices]
		elif self.minDistType.upper() == "H":
			relIndices = [hyIndicesA, hyIndicesB]
		else:
			raise ValueError("{} is an invalid value for minDistType".format(self.minDistType))

		#Get planar distances for each
		planarBinValGetters = [_PlanarDistsGetOneDimValsToBin(self.planeEqn, indices) for indices in relIndices]
		planarBinVals = [x.getValsToBin(sparseMatrixCalculator) for x in planarBinValGetters]

		#Find the minimum values
		outVals = list()
		for vals in it.zip_longest(*planarBinVals):
			outVals.append(min(vals))

		return outVals


class _DiscHBondCounterBetweenGroupsWithOxyDistFilterOneDimValGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, oxyIndices, hyIndices, distFilterIndices, distFilterVals, acceptor=True, donor=True, maxOO=3.5, maxAngle=35):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			distFilterIndices: (iter of ints) These are the indices we look at X-O distances for. They are used to divide water molecules into two groups
			distFilterVals: (iter) min/max distances to be in the two groups. [ [0,3], [3,5] ] would mean groupA are 0<=x<3 from distFilterIndices while groupB are 3<=x<5 from them
			acceptor: (Bool) If True we calculate dists/angles required to count number of groupA acceptors from groupB. Note changing the order of distFilterValues will give the reverse info (groupB acceptors from groupA)
			donor: (Bool) If True we calculate dists/angles required to count number of groupA donors to groupB. Note changing the order of distFilterValues will give the reverse info (groupB donors to groupA)
			maxOO: (float) The maximum O-O distance between two hydrogen-bonded water. Angles are only calculated when this criterion is fulfilled
			maxAngle: (float) The maximum OA-OD-HD angle for a hydrogen bond; OA = acceptor oxygen, OD=Donor oxygen, HD=donor hydrogen

		Notes:
			a) donor/acceptor: I think setting both will probably end up summing the number of donor and acceptors

		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.distFilterIndices = distFilterIndices
		self.distFilterVals = distFilterVals
		self.acceptor = acceptor
		self.donor = donor
		self.maxOO = maxOO
		self.maxAngle = maxAngle

	def getValsToBin(self, sparseMatrixCalculator):
		distMatrix = sparseMatrixCalculator.outDict["distMatrix"]
		angleMatrix = sparseMatrixCalculator.outDict["angleMatrix"]

		#1) Divide into groups 
		currArgs = [distMatrix, self.oxyIndices, self.distFilterIndices, self.distFilterVals]
		groupAOxyIndices, groupBOxyIndices = atomComboPopulatorHelp._groupOxyByDistMatrix(*currArgs)

		#2) Count the number of hydrogen bonds from groupA TO groupB + output in a 1-dim list
		outVals = list()
		for lIdx,oxyIdx in enumerate(self.oxyIndices):
			currVal = self._getValForOneOxyIdx( lIdx, groupAOxyIndices, groupBOxyIndices, distMatrix, angleMatrix )
			outVals.append(currVal)

		return outVals

	def _getValForOneOxyIdx(self, oxyListIdx, groupAIndices, groupBIndices, distMatrix, angleMatrix):
		oxyIdx = self.oxyIndices[oxyListIdx]

		#We only calc groupA -> groupB; thus the index has to be in groupA to be non-zero here
		if oxyIdx not in groupAIndices:
			return 0

		outVal = 0
		for oxyIdxB in groupBIndices:
			listIdxB = self.oxyIndices.index(oxyIdxB)
			currDist = distMatrix[oxyIdx][oxyIdxB]
			if currDist < self.maxOO and oxyIdxB!=oxyIdx:

				if self.acceptor:
					hyIdxA, hyIdxB = self.hyIndices[listIdxB]
					angleA, angleB = angleMatrix[oxyIdx][oxyIdxB][hyIdxA], angleMatrix[oxyIdx][oxyIdxB][hyIdxB]
					if angleA < self.maxAngle:
						outVal += 1
					if angleB < self.maxAngle:
						outVal += 1

				if self.donor:
					hyIdxA, hyIdxB = self.hyIndices[oxyListIdx]
					angleA, angleB = angleMatrix[oxyIdxB][oxyIdx][hyIdxA], angleMatrix[oxyIdxB][oxyIdx][hyIdxB]
					if angleA < self.maxAngle:
						outVal += 1
					if angleB < self.maxAngle:
						outVal += 1

		return outVal


