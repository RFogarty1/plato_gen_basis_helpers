
import collections
import copy
import itertools as it

import numpy as np


from . import atom_combo_core as atomComboCoreHelp
from . import calc_dists as calcDistsHelp
from . import water_rotations as waterRotHelp

from ..shared import simple_vector_maths as vectHelp



class _DistMatrixPopulator(atomComboCoreHelp._SparseMatrixPopulator):
	""" Populator meant for calculating distances between sets of indices """

	def __init__(self, fromIndices, toIndices, level=0):
		""" Initializer
		
		Args: (to/from is probably an arbitrary distinction here)
			fromIndices: (iter of ints) Indices of atoms to calculate distances from
			toIndices: (iter of ints) Indices of atoms to calculate distances to
				 
		"""
		self.fromIndices = fromIndices
		self.toIndices = toIndices
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
			useMatrix = outDict["distMatrix"]
		except KeyError:
			currKwargs = {"indicesA":self.fromIndices, "indicesB":self.toIndices, "sparseMatrix":True}
			outDict["distMatrix"] = calcDistsHelp._calcDistanceMatrixForCell_minImageConv_memoized(inpGeom, **currKwargs)
		else:
			self._populatePartiallyPopulatedMatrix(inpGeom, outDict)

	#Note: Figuring out which elements were already populated, and avoiding recalculating, was waaaaay too slow
	#Possibly faster way is here, but i couldnt get it to work https://stackoverflow.com/questions/8317022/get-intersecting-rows-across-two-2d-numpy-arrays
	def _populatePartiallyPopulatedMatrix(self, inpGeom, outDict):
		useMatrix = outDict["distMatrix"]

		#Get relevant indices (see note above)
		forwardFromIndices = self.fromIndices
		forwardToIndices = self.toIndices

		#Calculate only for the relevant indices
		currKwargs = {"indicesA":forwardFromIndices, "indicesB":forwardToIndices, "sparseMatrix":True}
		newSparseMatrix = calcDistsHelp._calcDistanceMatrixForCell_minImageConv_memoized(inpGeom, **currKwargs)

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


class _HozDistMatrixPopulator(atomComboCoreHelp._SparseMatrixPopulator):
	""" Populator meant for calculating in-surface-plane distances between groups of atoms; like a lower dimensional (2-d vs 3-d) rdf """

	def __init__(self, fromIndices, toIndices):
		""" Initializer
		
		Args: (to/from is probably an arbitrary distinction here)
			fromIndices: (iter of ints) Indices of atoms to calculate distances from
			toIndices: (iter of ints) Indices of atoms to calculate distances to
				 
		"""
		self.fromIndices= fromIndices
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
			outDict["hozDistMatrix"]
		except KeyError:
			currKwargs = {"indicesA":self.fromIndices, "indicesB":self.toIndices, "sparseMatrix":True}
			outDict["hozDistMatrix"] = calcDistsHelp._calcHozDistMatrixForCell_minImageConv_memoized_debuggable(inpGeom, **currKwargs)
		else:
			self._populatePartiallyPopulatedMatrix(inpGeom, outDict)

	def _populatePartiallyPopulatedMatrix(self, inpGeom, outDict):
		useMatrix = outDict["hozDistMatrix"]

		#Calculate for all indices; too slow to filter down to those needed
		currKwargs = {"indicesA":self.fromIndices, "indicesB":self.toIndices, "sparseMatrix":True}
		newSparseMatrix = calcDistsHelp._calcHozDistMatrixForCell_minImageConv_memoized_debuggable(inpGeom, **currKwargs)

		outDict["hozDistMatrix"] = np.where( np.isnan(useMatrix), newSparseMatrix, useMatrix )


class _PlanarDistMatrixPopulator(atomComboCoreHelp._SparseMatrixPopulator):
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


class _DiatomAngleWithVectorPopulator(atomComboCoreHelp._SparseMatrixPopulator):
	""" Populator meant for calculating angles for diatoms (e.g. hydroxyl O->H vector) with arbitrary vectors

	"""

	def __init__(self, diatomIndices, inpVector, level=0):
		""" Initializer
		
		Args:
			diatomIndices: (iter of len-2 int-iters) e.g. [ [0,1], [4,5] ] 
			inpVector: (len-3 float iter) Vector we want the angle with

		"""
		self.diatomIndices = diatomIndices
		self.inpVector = inpVector
		self.level = level

	@property
	def maxLevel(self):
		return self.level

	def populateMatrices(self, inpGeom, outDict, level):
		if level != self.level:
			pass
		else:
			self._populateMatrices(inpGeom, outDict)

	def _populateMatrices(self, inpGeom, outDict):
		try:
			matrices = outDict["diatom_arbitrary_angle_matrices"]
		except KeyError:
			self._populateForNoMatricesExistCase(inpGeom, outDict)
		else:
			self._populateForMatricesExistCase(inpGeom, outDict)

	def _populateForNoMatricesExistCase(self, inpGeom, outDict):
		outMatrices = list()
		currDict = {"inpVector":self.inpVector}

		currArgs = [inpGeom, self.diatomIndices, self.inpVector]
		currDict["matrix"] = calcDistsHelp.calcSparseDiatomAngleWithArbVectorMatrix(*currArgs)
		outDict["diatom_arbitrary_angle_matrices"] = [currDict]

	def _populateForMatricesExistCase(self, inpGeom, outDict):
		#Check if matrix is present
		relIdx = self._getRelevantMatrixIdx(outDict)
		if relIdx is None:
			self._populateIfRelevantMatrixMissing(inpGeom, outDict)
		else:
			self._populateIfRelevantMatrixExists(inpGeom, outDict)

	#NOTE: Will be efficient ONLY if populating matrix is faster than checking for overlapping indices, which is unlikely here
	#[it works for other cases since pre-compiled C-code is used to build the matrices]
	def _populateIfRelevantMatrixExists(self, inpGeom, outDict):
		newMatrix = calcDistsHelp.calcSparseDiatomAngleWithArbVectorMatrix(inpGeom, self.diatomIndices, self.inpVector)
		relIdx = self._getRelevantMatrixIdx(outDict)

		useMatrix = outDict["diatom_arbitrary_angle_matrices"][relIdx]["matrix"]
		outMatrix = np.where( np.isnan(useMatrix), newMatrix, useMatrix )
		outDict["diatom_arbitrary_angle_matrices"][relIdx]["matrix"] = outMatrix


	def _populateIfRelevantMatrixMissing(self, inpGeom, outDict):
		raise NotImplementedError("")
 
	def _getRelevantMatrixIdx(self, outDict):
		matrixDicts = outDict["diatom_arbitrary_angle_matrices"]
		relIdx = None
		for idx,currDict in enumerate(matrixDicts):
			if np.allclose( np.array(self.inpVector), np.array(currDict["inpVector"]) ):
				relIdx = idx
				break
		return relIdx

#Populate roll/pitch/azimuthal matrices simultaneously
#water_orientations_roll, water_orientations_pitch, water_orientations_azimuthal
class _WaterOrientationPopulator(atomComboCoreHelp._SparseMatrixPopulator):


	def __init__(self, oxyIndices, hyIndices):
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices

	def populateMatrices(self, inpGeom, outDict, level):
		if level == 0:
			self._populateLevelZeroMatrices(inpGeom, outDict)
		elif level == 1:
			self._populateLevelOneMatrices(inpGeom, outDict)

	def _populateLevelZeroMatrices(self, inpGeom, outDict):
		try:
			unused = outDict["pos_vector_matrix"]
		except KeyError:
			self._populateLevelZeroMatrices_nonePresent(inpGeom, outDict)
		else:
			self._populateLevelZeroMatrices_partiallyPresent(inpGeom, outDict)

	def _populateLevelZeroMatrices_nonePresent(self, inpGeom, outDict):
		idxPairs = self._getOxyAndHyIdxPairs()
		posVectorMatrix = calcDistsHelp.getNearestImageVectorsForIdxPairs(inpGeom, idxPairs, sparseMatrix=True)
		outDict["pos_vector_matrix"] = posVectorMatrix

	def _populateLevelZeroMatrices_partiallyPresent(self, inpGeom, outDict):
		outMatrix = outDict["pos_vector_matrix"]

		#Filter out the pairs that have already been calculated
		updateIdxPairs = list()
		for idxPair in self._getOxyAndHyIdxPairs():
			if np.isnan(outMatrix[tuple(idxPair)][0]):
				updateIdxPairs.append(idxPair)

		#Get the values for the updated and put in the matrix
		if len(updateIdxPairs)==0:
			return 0

		newVals = calcDistsHelp.getNearestImageVectorsForIdxPairs(inpGeom, updateIdxPairs, sparseMatrix=False)
		for idx,val in it.zip_longest(updateIdxPairs, newVals):
			outMatrix[tuple(idx)] = val


	def _getOxyAndHyIdxPairs(self):
		idxPairs = [ [ [oxyIdx,hyIndices[0]], [oxyIdx,hyIndices[1]] ] for oxyIdx, hyIndices in it.zip_longest(self.oxyIndices,self.hyIndices) ]
		idxPairs = [ idxPair for idxPair in it.chain(*idxPairs) ]
		return sorted(idxPairs)

	def _populateLevelOneMatrices(self, inpGeom, outDict):
		try:
			unused = outDict["water_rotations_roll_matrix"]
		except:
			self._populateLevelOneMatrices_nonePresent(inpGeom, outDict)
		else:
			self._populateLevelOneMatrices_partiallyPresent(inpGeom, outDict)

	def _populateLevelOneMatrices_nonePresent(self, inpGeom, outDict):
		rollMatrix, pitchMatrix, azimuthalMatrix = self._getOrientationMatricesFromOxyIndices(inpGeom, self.oxyIndices, outDict)
		outDict["water_rotations_roll_matrix"] = rollMatrix
		outDict["water_rotations_pitch_matrix"] = pitchMatrix
		outDict["water_rotations_azimuthal_matrix"] = azimuthalMatrix

	def _populateLevelOneMatrices_partiallyPresent(self, inpGeom, outDict):
		#Get matrices
		inpRollMatrix, inpPitchMatrix = outDict["water_rotations_roll_matrix"], outDict["water_rotations_pitch_matrix"]
		inpAziMatrix = outDict["water_rotations_azimuthal_matrix"]

		#Get OXYGEN indices of is nan
		rollNan, pitchNan, aziNan = [ np.argwhere( np.isnan(matrix) ) for matrix in [inpRollMatrix,inpPitchMatrix,inpAziMatrix]]
		allNan = np.union1d(rollNan, pitchNan)
		allNan = np.union1d(allNan,aziNan)
		neededIndices = np.intersect1d( np.array(self.oxyIndices), allNan )

		if len(neededIndices)==0:
			return 0

		#Calculate relevant values then update all the previous matrices
		updateRoll, updatePitch, updateAzi = self._getOrientationMatricesFromOxyIndices(inpGeom, neededIndices, outDict)

		outRollMatrix = np.where( np.isnan(inpRollMatrix), updateRoll, inpRollMatrix )
		outPitchMatrix = np.where( np.isnan(inpPitchMatrix), updatePitch, inpPitchMatrix )
		outAziMatrix = np.where( np.isnan(inpAziMatrix), updateAzi, inpAziMatrix )

		outDict["water_rotations_roll_matrix"] = outRollMatrix
		outDict["water_rotations_pitch_matrix"] = outPitchMatrix
		outDict["water_rotataions_azimuthal_matrix"] = outAziMatrix


	def _getOrientationMatricesFromOxyIndices(self, inpGeom, oxyIndices, outDict):
		#Initialize output matrices (really just vectors)
		outDim = len(inpGeom.cartCoords)
		rollMatrix, pitchMatrix, azimuthalMatrix = np.empty((outDim)), np.empty((outDim)), np.empty((outDim))
		rollMatrix[:], pitchMatrix[:], azimuthalMatrix[:] = np.nan, np.nan, np.nan


		#
		posVectMatrix = outDict["pos_vector_matrix"]
		outRotMatrices = list()
		for oxyIdx in oxyIndices:
			idxInList = self.oxyIndices.index(oxyIdx)
			hyIndices = self.hyIndices[idxInList]
			ohVectA, ohVectB = posVectMatrix[oxyIdx][hyIndices[0]], posVectMatrix[oxyIdx][hyIndices[1]]
			ohVectA, ohVectB = [vectHelp.getUnitVectorFromInpVector(x) for x in [ohVectA,ohVectB]]
			currRotMatrix = waterRotHelp._getStandardRotationMatrixFromTwoOHVectors(ohVectA,ohVectB)
			outRotMatrices.append(currRotMatrix)

		#Get the rotation angles in order of the oxygen indices
		allRotAngles = waterRotHelp._getWaterStandardRotationCoordsFromMatrices(outRotMatrices)
		for rotIdx,oxyIdx in enumerate(oxyIndices):
			rollMatrix[oxyIdx] = allRotAngles[rotIdx][0]
			pitchMatrix[oxyIdx] = allRotAngles[rotIdx][1]
			azimuthalMatrix[oxyIdx] = allRotAngles[rotIdx][2] 

		return rollMatrix, pitchMatrix, azimuthalMatrix

	@property
	def maxLevel(self):
		return 1

class _WaterMinDist_plusMinDistFilter_populator(atomComboCoreHelp._SparseMatrixPopulator):

	def __init__(self, oxyIndices, hyIndices, toIndices, filterToIndices, minDistType):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			toIndices: (iter of ints) The indices of atoms we calculate the minimum distance TO
			filterToIndices: (iter of ints) 
			minDistType: (str) Controls which atoms to get the minimum distance from. Current options are "all","o", and "h" (case insensitive)

		"""
		self.level = 0
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.toIndices = toIndices
		self.filterToIndices = filterToIndices
		self.minDistType = minDistType

	@property
	def maxLevel(self):
		return self.level

	def populateMatrices(self, inpGeom, outDict, level):
		#1) Sort out the water-toIndices bit
		populatorA = _WaterMinDistPopulator(self.oxyIndices, self.hyIndices, self.toIndices, self.minDistType)
		populatorA.populateMatrices(inpGeom, outDict, level)

		#2) Sort out getting all distances we need to filter to "toIndices" down
		populatorB = _DistMatrixPopulator(self.toIndices, self.filterToIndices, level=self.level)
		populatorB.populateMatrices(inpGeom, outDict, level)



class _WaterMinDistPopulator(atomComboCoreHelp._SparseMatrixPopulator):

	def __init__(self, oxyIndices, hyIndices, toIndices, minDistType):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			toIndices: (iter of ints) The indices of atoms we calculate the minimum distance TO
			minDistType: (str) Controls which atoms to get the minimum distance from. Current options are "all","o", and "h" (case insensitive)

		"""
		self.level = 0
		self.oxyIndices = oxyIndices 
		self.hyIndices = hyIndices
		self.toIndices = toIndices
		self.minDistType = minDistType

	@property
	def maxLevel(self):
		return self.level

	def populateMatrices(self, inpGeom, outDict, level):
		fromIndices = self._getFromIndices()
		populator = _DistMatrixPopulator(fromIndices, self.toIndices, level=self.level)
		populator.populateMatrices(inpGeom,outDict,level)

	def _getFromIndices(self):
		if self.minDistType.upper()=="ALL":
			outIndices = self.oxyIndices + [x for x in it.chain(*self.hyIndices)]
		elif self.minDistType.upper()=="O":
			outIndices = self.oxyIndices
		elif self.minDistType.upper()=="H":
			outIndices = [x for x in it.chain(*self.hyIndices)]
		else:
			raise ValueError("{} is an invalid value for self.minDistType".format(self.minDistType))

		return outIndices

class _WaterPlanarDistPopulator(atomComboCoreHelp._SparseMatrixPopulator):


	def __init__(self, oxyIndices, hyIndices, planeEqn, primaryIdxType="O", primaryOnly=True):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			planeEqn: (ThreeDimPlaneEquation) The plane equation to calculate distance distribution from
			primaryIdxType: (str) The element of the primary index. "O", "Ha" and "Hb" are the standard options
			primaryOnly: (Bool) If True populate only primaryIdxType; else populate distances for ALL water atoms 

		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.planeEqn = planeEqn
		self.primaryIdxType = primaryIdxType
		self.primaryOnly = primaryOnly
		self.level = 0

	@property
	def maxLevel(self):
		return self.level

	@property
	def primaryIndices(self):
		if self.primaryIdxType.upper() == "O":
			return self.oxyIndices
		elif self.primaryIdxType.upper() == "HA":
			return [x[0] for x in self.hyIndices]
		elif self.primaryIdxType.upper() == "HB":
			return [x[1] for x in self.hyIndices]
		else:
			raise ValueError("primaryIdxType = {} is an invalid value".format(self.primaryIdxType))

	def populateMatrices(self, inpGeom, outDict, level):
		if self.primaryOnly:
			populator = _PlanarDistMatrixPopulator(self.primaryIndices, self.planeEqn, level=0)
		else:
			outIndices = self.oxyIndices + [x for x in it.chain(*self.hyIndices)]
			populator = _PlanarDistMatrixPopulator(outIndices, self.planeEqn, level=0)

		populator.populateMatrices(inpGeom, outDict, 0)

class _CountHBondsBetweenGenericGroupsPopulator(atomComboCoreHelp._SparseMatrixPopulator):

	def __init__(self, fromNonHyIndices, fromHyIndices, toNonHyIndices, toHyIndices, acceptor=True, donor=True, maxOO=3.5):
		""" Initializer
		
		Args:
			fromNonHyIndices: (iter of iter of ints) Each entry corresponds to an iter of indices for non-hydrogen atoms (h-bond acceptor atoms) on each molecule
			fromHyIndices: (iter of iter of ints) Each entry corresponds to an iter of indices for hydrogen atoms on each molecule
			toNonHyIndices: (iter of iter of ints) Same as fromNonHyIndices, except for molecules of the second group
			toHyIndices: (iter of iter of ints) Same as fromHyIndices, except for molecules of the second group
			acceptor: (Bool) If True we calculate dists/angles required to count number of groupA acceptors from groupB
			donor: (Bool) If True we calculate dists/angles required to count number of groupA donors to groupB
			maxOO: (float) The maximum X-X distance between two hydrogen-bonded water. For water X are the oxygen atoms; hence the variable name. Angles are only calculated when this criterion is fulfilled
			maxAngle: (float) The maximum XA-XD-XD angle for a hydrogen bond; OA = acceptor oxygen, OD=Donor oxygen, HD=donor hydrogen
	 
		NOTE:
			Don't have multiple NonHyIndices in one entry unless there are no hyIndices. For example, it would be fine to use both oxygen in CO2, but not for (HO)2-CO since theres no way to know which hydrogen is connected to each oxygen

		"""
		self.fromNonHyIndices = fromNonHyIndices
		self.fromHyIndices = fromHyIndices
		self.toNonHyIndices = toNonHyIndices
		self.toHyIndices = toHyIndices
		self.acceptor = acceptor
		self.donor = donor
		self.maxOO = maxOO
		self.nanMatrix = False

	@property
	def maxLevel(self):
		return 1 #level 0 for dist matrices, level 1 for angle matrix

	def populateMatrices(self, inpGeom, outDict, level):
		if level==0:
			self._populateDistMatrices(inpGeom, outDict)
		elif level==1:
			self._populateAngleMatrices(inpGeom, outDict)
		else:
			pass

	def _populateDistMatrices(self, inpGeom, outDict):
		level = 0
		populatorFromIndices = [x for x in it.chain(*self.fromNonHyIndices)]
		populatorToIndices = [x for x in it.chain(*self.toNonHyIndices)]
 
		distPopulator = _DistMatrixPopulator(populatorFromIndices, populatorToIndices)
		distPopulator.populateMatrices(inpGeom, outDict, level)

	def _populateAngleMatrices(self, inpGeom, outDict):
		try:
			unused = outDict["angleMatrix"]
		except KeyError:
			self._populateAngleMatrixWhenNonePresent(inpGeom, outDict)
		else:
			self._populateAngleMatrixWhenSomePresent(inpGeom, outDict)

	def _populateAngleMatrixWhenSomePresent(self, inpGeom, outDict):
		useMatrix = outDict["angleMatrix"]
		outAngleIndices = self._getFullOutAngleIndicesRequired(inpGeom, outDict)

		#Populate output matrix
		allNewAngles = calcDistsHelp.getInterAtomicAnglesForInpGeom(inpGeom, outAngleIndices)

		for currIdx,currAngle in it.zip_longest(outAngleIndices, allNewAngles):
			useMatrix[tuple(currIdx)] = currAngle

	def _populateAngleMatrixWhenNonePresent(self, inpGeom, outDict):
		#Initialise the matrix
		cartCoords = inpGeom.cartCoords
		outMatrix = np.empty( (len(cartCoords), len(cartCoords), len(cartCoords)) )
		if self.nanMatrix:
			outMatrix[:] = np.nan #Almost the FULL runtime for large systems (since the matrix gets SO large + its so sparsely populated)

		#Get the angle indices (fast)
		outAngleIndices = self._getFullOutAngleIndicesRequired(inpGeom, outDict)

		#Get the angles we need (seems to take ~none of the runtime)
		outAngles = calcDistsHelp.getInterAtomicAnglesForInpGeom(inpGeom, outAngleIndices)

		#Generally quite quick
		for currIdx, currAngle in it.zip_longest(outAngleIndices, outAngles):
			outMatrix[tuple(currIdx)] = currAngle

		#Populate
		outDict["angleMatrix"] = outMatrix


	def _getFullOutAngleIndicesRequired(self, inpGeom, outDict):
		distMatrix = outDict["distMatrix"]

		#Get all possible angle indices
		currArgs = [self.fromNonHyIndices, self.toNonHyIndices, self.fromHyIndices, self.toHyIndices, self.maxOO, self.acceptor, self.donor, distMatrix]
		outAngleIndices = _getFullOutAngleIndicesRequired_GENERIC(*currArgs)
		return outAngleIndices



class _DiscHBondCounterBetweenGroupsWithOxyDistFilterPopulator(atomComboCoreHelp._SparseMatrixPopulator):
	""" Populates matrices to allow number of hydrogen bonds between specific groups of water indices to be counted """

	def __init__(self, oxyIndices, hyIndices, distFilterIndices, distFilterVals, acceptor=True, donor=True, maxOO=3.5):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			distFilterIndices: (iter of ints) These are the indices we look at X-O distances for. They are used to divide water molecules into two groups
			distFilterVals: (iter) min/max distances to be in the two groups. [ [0,3], [3,5] ] would mean groupA are 0<=x<3 from distFilterIndices while groupB are 3<=x<5 from them
			acceptor: (Bool) If True we calculate dists/angles required to count number of groupA acceptors from groupB. Note changing the order of distFilterValues will give the reverse info (groupB acceptors from groupA)
			donor: (Bool) If True we calculate dists/angles required to count number of groupA donors to groupB. Note changing the order of distFilterValues will give the reverse info (groupB donors to groupA)
			maxOO: (float) The maximum O-O distance between two hydrogen-bonded water. Angles are only calculated when this criterion is fulfilled

		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.distFilterIndices = distFilterIndices
		self.distFilterVals = distFilterVals
		self.acceptor = acceptor
		self.donor = donor
		self.maxOO = maxOO
		#Optimisation flag that I doubt i'll ever access. Safer to set to True but so much slower
		#If True we initialise our sparse matrix to be full of NaN; else we just leave it as "empty" (i.e. random/undefined) values. 
		self.nanMatrix = False 

	@property
	def maxLevel(self):
		return 1 #level 0 for dist matrices, level 1 for angle matrix
	

	def populateMatrices(self, inpGeom, outDict, level):
		if level==0:
			self._populateDistMatrices(inpGeom, outDict)
		elif level==1:
			self._populateAngleMatrices(inpGeom, outDict)
		else:
			pass

	def _populateDistMatrices(self, inpGeom, outDict):
		level = 0 
		distPopulatorA = _DistMatrixPopulator(self.oxyIndices, self.oxyIndices)
		distPopulatorB = _DistMatrixPopulator(self.oxyIndices, self.distFilterIndices)
		distPopulatorA.populateMatrices(inpGeom, outDict, level)
		if self.distFilterIndices is not None:
			distPopulatorB.populateMatrices(inpGeom, outDict, level)


	def _populateAngleMatrices(self, inpGeom, outDict):
		try:
			unused = outDict["angleMatrix"]
		except:
			self._populateAngleMatrixWhenNonePresent(inpGeom, outDict)
		else:
			self._populateAngleMatrixWhenSomePresent(inpGeom, outDict)


	def _populateAngleMatrixWhenSomePresent(self, inpGeom, outDict):

		#Get all possible angle indices
		outAngleIndices = self._getFullOutAngleIndicesRequired(inpGeom, outDict)

		#Filter based on which are already present
		useMatrix = outDict["angleMatrix"]

		#Note: The filter thing will break if self.nanMatrix is False; this caused a very annoying bug
		#Trying to filter previously calculated angles is generally unlikely to help speed-wise regardless
		filteredOutIndices = outAngleIndices

		#populate the output matrix
		allNewAngles = calcDistsHelp.getInterAtomicAnglesForInpGeom(inpGeom, filteredOutIndices)

		for currIdx, currAngle in it.zip_longest(filteredOutIndices, allNewAngles):
			useMatrix[tuple(currIdx)] = currAngle



	def _populateAngleMatrixWhenNonePresent(self, inpGeom, outDict):
		#Iniitialise the matrix
		cartCoords = inpGeom.cartCoords
		outMatrix = np.empty( (len(cartCoords), len(cartCoords), len(cartCoords)) )
		if self.nanMatrix:
			outMatrix[:] = np.nan #Almost the FULL runtime for large systems (since the matrix gets SO large + its so sparsely populated)

		#Get the angle indices (fast)
		outAngleIndices = self._getFullOutAngleIndicesRequired(inpGeom, outDict)

		#Get the angles we need (seems to take ~none of the runtime)
		outAngles = calcDistsHelp.getInterAtomicAnglesForInpGeom(inpGeom, outAngleIndices)

		#Generally quite quick
		for currIdx, currAngle in it.zip_longest(outAngleIndices, outAngles):
			outMatrix[tuple(currIdx)] = currAngle

		#Populate
		outDict["angleMatrix"] = outMatrix


	def _getFullOutAngleIndicesRequired(self, inpGeom, outDict):
		distMatrix = outDict["distMatrix"]

		#Figure out which oxy indices are in each group
		assert len(self.distFilterVals)==2
		groupAOxyIndices, groupBOxyIndices = _groupOxyByDistMatrix(distMatrix, self.oxyIndices, self.distFilterIndices, self.distFilterVals)

		#Filter which oxygens to calculate angles between (based on their separation)
		outOxyPairs = list()
		for idxA in groupAOxyIndices:
			for idxB in groupBOxyIndices:
				if distMatrix[idxA][idxB] < self.maxOO:
					outOxyPairs.append( [idxA,idxB] )

		#Get the relevant hydrogen indices
		groupAHyIndices, groupBHyIndices = list(),list()
		for idxA in groupAOxyIndices:
			currIdx = self.oxyIndices.index(idxA)
			groupAHyIndices.append( self.hyIndices[currIdx] )

		for idxB in groupBOxyIndices:
			currIdx = self.oxyIndices.index(idxB)
			groupBHyIndices.append( self.hyIndices[currIdx] )

		currArgs = [ [[x] for x in groupAOxyIndices], [[x] for x in groupBOxyIndices], groupAHyIndices, groupBHyIndices, self.maxOO, self.acceptor, self.donor, distMatrix]
		outAngleIndices = _getFullOutAngleIndicesRequired_GENERIC(*currArgs)

		return outAngleIndices


#The GENERIC in caps was originally to differentiate it from a version that worked only for water.
#Though i removed that so i should probably change name of this one at some point
def _getFullOutAngleIndicesRequired_GENERIC(groupANonHyIndices, groupBNonHyIndices, hyIndicesA, hyIndicesB, maxOO, countAcceptor, countDonor, distMatrix):


	def _getAngleIndicesForSingleGroupABPair(groupANonHy, groupBNonHy, groupAHy, groupBHy):
		outAngleIndices = list()

		#0) Special case - small optimisation. 
		if (len(groupANonHy)==1) and (len(groupBNonHy)==1):
			currDist = distMatrix[groupANonHy[0]][groupBNonHy[0]]
			if currDist >= maxOO:
				return outAngleIndices

		#1) Do the error checks
		_checkGenericNonHyAndHyIndicesHaveValidValues(groupANonHy, groupAHy)
		_checkGenericNonHyAndHyIndicesHaveValidValues(groupBNonHy, groupBHy)

		#2) Figure out what angles we need
		if countDonor:
			donorNonHy, donorHy = groupANonHy, groupAHy
			donorPairs = [ [nonHy,hy] for nonHy,hy in it.product(donorNonHy,donorHy) ]
			acceptorNonHy = groupBNonHy
			for donorPair, acceptor in it.product(donorPairs, acceptorNonHy):
				currDist = distMatrix[donorPair[0]][acceptor]
				if currDist < maxOO:
					outAngleIndices.append( [acceptor] + donorPair ) 

		#Basically same code as directly above...
		if countAcceptor:
			donorNonHy, donorHy = groupBNonHy, groupBHy
			donorPairs = [ [nonHy,hy] for nonHy,hy in it.product(donorNonHy,donorHy) ]
			acceptorNonHy = groupANonHy
			for donorPair,acceptor in it.product(donorPairs, acceptorNonHy):
				currDist = distMatrix[donorPair[0]][acceptor]
				if currDist < maxOO:
					outAngleIndices.append( [acceptor] + donorPair )

		return outAngleIndices


	#Filter which non-hy to calculate angles between based on their separation
	outAngleIndices = list()
	for listIdxA, nonHyIndicesA in enumerate(groupANonHyIndices):
		for listIdxB, nonHyIndicesB in enumerate(groupBNonHyIndices):
			currArgs = [nonHyIndicesA, nonHyIndicesB, hyIndicesA[listIdxA], hyIndicesB[listIdxB] ]
			currAngleIndices = _getAngleIndicesForSingleGroupABPair(*currArgs)
			outAngleIndices.extend(currAngleIndices)

	return outAngleIndices

#NOTE: Not explicitly unit tested at time of writing
#Also theres more calls to copy than required
#class _MemoizationWrapperForGettingOutAngleIndices():
#
#	def __init__(self, maxEntries=20):
#		self.argCache = collections.deque(list(), maxlen=maxEntries)
#		self.funct = _getFullOutAngleIndicesRequired_GENERIC
#
#	def __call__(self, groupANonHyIndices, groupBNonHyIndices, hyIndicesA, hyIndicesB, maxOO, countAcceptor, countDonor, distMatrix):
#		#1) Check if args already in cache
#		inpArgs = [groupANonHyIndices, groupBNonHyIndices, hyIndicesA, hyIndicesB, maxOO, countAcceptor, countDonor, distMatrix]
#		directCmpIndices = [0,1,2,3,4,5,6]
#		npCmpIndices = [7]
#
#		for prevArgs,output in self.argCache:
#			match = True
#			#simple to compare args
#			for idx,val in enumerate(prevArgs):
#				if idx in directCmpIndices:
#					if val!=inpArgs[idx]:
#						match = False
#						break
#				elif idx in npCmpIndices:
#					if not np.allclose(val,inpArgs[idx]):
#						match = False
#						break
#				else:
#					raise ValueError("Should never hit this point")
#
#			if match:
#				return copy.deepcopy(output)
#
#		#2) If not, run the function and add to cache
#		output = self.funct(*inpArgs)
#		toAppend = [inpArgs, copy.deepcopy(output)]
#		self.argCache.appendleft(toAppend)
#
#		return copy.deepcopy(output)


#Jus a single entry; this covers most cases where a speedup is possible + comparing distance matrices is SLOW
#_memoized_getFullOutAngleIndicesRequired_GENERIC = _MemoizationWrapperForGettingOutAngleIndices( maxEntries=1 )


def _checkGenericNonHyAndHyIndicesHaveValidValues(nonHyIndices, hyIndices):
	""" Make sure we dont have hyIndices when we have multiple nonHyIndices """
	lenGroupNonHy, lenGroupHy = len(nonHyIndices), len(hyIndices)
	if lenGroupNonHy>1 and lenGroupHy>0:
		raise ValueError("nonHyIndices = {}, hyIndices = {}; cant have ANY hyIndices for >1 nonHyIndices".format(nonHyIndices,hyIndices))


def _groupOxyByDistMatrix(distMatrix, oxyIndices, distFilterIndices, distFilterVals):
	""" Groups oxygen atoms based on MINIMUM distance from distFilterIndices
	
	Args:
		distMatrix: (NxN matrix) Contains distances between atoms. distMatrix[idxA][idxB] contains distance between idxA and idxB atoms
		oxyIndices: (iter of ints) The oxygen indices for each water molecule
		distFilterIndices: (iter of ints) These are the indices we look at X-O distances for. They are used to divide water molecules into two groups
		distFilterVals: (iter) min/max distances to be in the two groups. [ [0,3], [3,5] ] would mean groupA are 0<=x<3 from distFilterIndices while groupB are 3<=x<5 from them

	Returns
		groupIndices: (len-n iter of int iters) The oxy indices for each group. groupIndices[groupIdx] contains in the indices of the groupIdx values. The length of this is the same as distFilterIndices

	NOTES:
		a) I've only tested this for the two-group case so far
		b) Its possible for an index to be in multiple groups if they overlap. This is useful for looking at h-bonding within a group
		c) If you pass None for distFilterIndices then all oxyIndices will be placed in all groups; number of groups is still determined by distFilterVals though
 
	"""
	#Deal with the case where we dont filter anything
	if distFilterIndices is None:
		groupIndices = [ [idx for idx in oxyIndices] for vals in distFilterVals ]
		return groupIndices

	#Figure out the minimum distances
	minDists = list()
	for idx in oxyIndices:
		currDists = distMatrix[idx][:]
		minDists.append( np.nanmin(currDists[distFilterIndices]) )

	#Figure out which oxy indices are in which group
	groupIndices = [list() for x in distFilterVals]

	for lIdx,oxyIdx in enumerate(oxyIndices):
		currMinDist = minDists[lIdx]
		for filterIdx, filterVals in enumerate(distFilterVals):
			if (currMinDist>=filterVals[0]) and (currMinDist<filterVals[1]):
				groupIndices[filterIdx].append(oxyIdx)

	return groupIndices

