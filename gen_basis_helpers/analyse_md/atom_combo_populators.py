


import itertools as it

import numpy as np


from . import atom_combo_core as atomComboCoreHelp
from . import calc_dists as calcDistsHelp


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

	#Note: Figuring out which elements were already populated, and avoiding recalculating, was waaaaay too slow
	#Possibly faster way is here, but i couldnt get it to work https://stackoverflow.com/questions/8317022/get-intersecting-rows-across-two-2d-numpy-arrays
	def _populatePartiallyPopulatedMatrix(self, inpGeom, outDict):
		useMatrix = outDict["distMatrix"]

		#Get relevant indices (see note above)
		forwardFromIndices = self.fromIndices
		forwardToIndices = self.toIndices

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
#		filteredOutIndices = list()
#		for idx in outAngleIndices:
#			if np.isnan( useMatrix[tuple(idx)] ):
#				filteredOutIndices.append( idx )

#		#
#		if len(filteredOutIndices)==0:
#			return 0

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

#		#Get the minimum O-X distances
#		minDists = list()
#		for idx in self.oxyIndices:
#			currDists = distMatrix[idx][:]
#			minDists.append( np.nanmin(currDists[self.distFilterIndices]) )
#
#		#Figure out which oxy indices are in which group
#		groupAOxyIndices, groupBOxyIndices = list(), list()
#		assert len(self.distFilterVals)==2
#		for lIdx, oxyIdx in enumerate(self.oxyIndices):
#			currMinDist = minDists[lIdx]
#			if (currMinDist>=self.distFilterVals[0][0]) and (currMinDist<self.distFilterVals[0][1]):
#				groupAOxyIndices.append( oxyIdx )
#			elif (currMinDist>=self.distFilterVals[1][0]) and (currMinDist<self.distFilterVals[1][1]):
#				groupBOxyIndices.append( oxyIdx )

		#Figure out which oxy indices are in each group
		assert len(self.distFilterVals)==2
		groupAOxyIndices, groupBOxyIndices = _groupOxyByDistMatrix(distMatrix, self.oxyIndices, self.distFilterIndices, self.distFilterVals)


		#Filter which oxygens to calculate angles between (based on their separation)
		outOxyPairs = list()
		for idxA in groupAOxyIndices:
			for idxB in groupBOxyIndices:
				if distMatrix[idxA][idxB] < self.maxOO:
					outOxyPairs.append( [idxA,idxB] )

		
		#Figure out the angles based on donor/acceptors
		#Note the order of angles is [O_acceptor, O_donor, H_donor]
		outAngleIndices = list()
		for oxyPair in outOxyPairs:
			
			#groupA (the left side one) is the acceptor
			if self.acceptor:
				oxyIdxD, oxyIdxA = oxyPair[1], oxyPair[0]
				donorListIdx = self.oxyIndices.index(oxyIdxD)
				hA, hB = self.hyIndices[donorListIdx][0], self.hyIndices[donorListIdx][1]
				angleA, angleB = [oxyIdxA, oxyIdxD, hA], [oxyIdxA, oxyIdxD, hB]
				outAngleIndices.append(angleA)
				outAngleIndices.append(angleB)

			if self.donor:
				oxyIdxD, oxyIdxA = oxyPair[0], oxyPair[1]
				donorListIdx = self.oxyIndices.index(oxyIdxD)
				hA, hB = self.hyIndices[donorListIdx][0], self.hyIndices[donorListIdx][1]
				angleA, angleB = [oxyIdxA, oxyIdxD, hA], [oxyIdxA, oxyIdxD, hB]
				outAngleIndices.append(angleA)
				outAngleIndices.append(angleB)

		return outAngleIndices



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

