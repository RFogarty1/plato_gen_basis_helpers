
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


#TODO: We will have a version that filters indicesA and indicesB
#Maybe need to combine with the classification methods...(though be careful of circular imports)
class _RadialDistsGetValsToBin(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, indicesA, indicesB, minVal=None, maxVal=None):
		""" Initializer
		
		Args:
			indicesA: (iter of ints) Indices to get rdf FROM
			indicesB: (iter of ints) Indices to get rdf TO
			minVal: (float) Minimum value allowed for binning (vals must be >= to)
			maxVal: (float) Maximum value allowes for binning (vals must be < than)

		"""
		self.indicesA = indicesA
		self.indicesB = indicesB
		self.minVal = minVal
		self.maxVal = maxVal

	def getValsToBin(self, sparseMatrixCalculator):
		relevantMatrix = sparseMatrixCalculator.outDict["distMatrix"]
		return _getRdfValsToBinFromMatrix(relevantMatrix, self.indicesA, self.indicesB, minVal=self.minVal, maxVal=self.maxVal)


class _HozDistsGetValsToBin(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, fromIndices, toIndices):
		""" Initializer
		
		Args:
			fromIndices: (iter of ints) Indices to get hoz-dists FROM
			toIndices: (iter of ints) Indices to get hoz-dists TO

		"""
		self.fromIndices = fromIndices
		self.toIndices = toIndices

	def getValsToBin(self, sparseMatrixCalculator):
		relevantMatrix = sparseMatrixCalculator.outDict["hozDistMatrix"]
		return _getRdfValsToBinFromMatrix(relevantMatrix, self.fromIndices, self.toIndices)


def _getRdfValsToBinFromMatrix(relevantMatrix, indicesA, indicesB, minVal=None, maxVal=None):
	#Worried this may be a bit slow
	outVals = np.empty( (len(indicesA), len(indicesB) ) )

	for outRowIdx,inpRowIdx in enumerate(indicesA):
		for outColIdx,inpColIdx in enumerate(indicesB):
			outVals[outRowIdx][outColIdx] = relevantMatrix[inpRowIdx][inpColIdx]

	if (minVal is None) and (maxVal is None):
		return outVals

	#Filter out any values that wont fit into bins anyway
	#Its quickest to do here if we're likely to do a multi-dimensional rdf (number of terms is N instead of dims^{N})
	#TODO: Need to explicitly test these.
	else:
		filteredOutVals = list()
		for rIdx,row in enumerate(outVals):
			if (minVal is not None) and (maxVal is not None):
				currRow = row[ (row<maxVal) & (row>=minVal) ]
			elif minVal is not None:
				currRow = row[ (row>=minVal) ]
			elif maxVal is not None:
				currRow = row[ (row<maxVal) ]
			else:
				raise ValueError("Should be impossible to reach here")
			filteredOutVals.append(currRow)

	return filteredOutVals


class _MinDistsGetOneDimValsToBin(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, fromIndices, toIndices, minVal=-0.01):
		""" Initializer
		
		Args:
			fromIndices: (iter of ints) The indices of atoms we calculate distances from. Our output bin values will be len(fromIndices)
			toIndices: (iter of ints) The indices of atoms we calculate distances to
			minVal: (float) The minimum value to take into account; this is to stop getting 0 when comparing two sets of overlapping groups
 
		"""
		self.fromIndices = fromIndices
		self.toIndices = toIndices
		self.minVal = minVal

	def getValsToBin(self, sparseMatrixCalculator):
		relevantMatrix = sparseMatrixCalculator.outDict["distMatrix"]
		return _getMinDistsToBinFromRelevantMatrix(relevantMatrix, self.fromIndices, self.toIndices, self.minVal)

	def __eq__(self, other):
		if type(other) is not type(self):
			return False

		directCmpAttrs = [ "fromIndices", "toIndices" ]
		for attr in directCmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if valA != valB:
				return False
		return True

class _MinHozDistsGetValsToBin(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, fromIndices, toIndices, minVal=-0.01):
		""" Description of function
		
		Args:
			fromIndices: (iter of ints) The indices of atoms we calculate distances from. Our output bin values will be len(fromIndices)
			toIndices: (iter of ints) The indices of atoms we calculate distances to
			minVal: (float) The minimum value to take into account; this is to stop getting 0 when comparing two sets of overlapping groups
				 
		"""
		self.fromIndices = fromIndices
		self.toIndices = toIndices
		self.minVal = minVal

	def getValsToBin(self, sparseMatrixCalculator):
		#Edge case: We want to return 0, regardless of self.minVal, if we only have the same index for to/from
		if len(self.fromIndices)==1 and len(self.toIndices)==1:
			if self.toIndices == self.toIndices:
				return [0]


		relevantMatrix = sparseMatrixCalculator.outDict["hozDistMatrix"]
		outVals = _getMinDistsToBinFromRelevantMatrix(relevantMatrix, self.fromIndices, self.toIndices, self.minVal)
		return outVals

def _getMinDistsToBinFromRelevantMatrix(relevantMatrix, fromIndices, toIndices, minVal):
	outVals = list()
	for idx in fromIndices:
		currDists = relevantMatrix[idx][:]
		relRow = currDists[toIndices]
		moddedRelRow = np.where(relRow>minVal, relRow, np.nan)
		outVals.append( np.nanmin(moddedRelRow) )

	return outVals


class _WaterOrientationBinValGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, oxyIndices, angleType):
		self.oxyIndices = oxyIndices
		self.angleType = angleType

	def getValsToBin(self, sparseMatrixCalculator):
		#Get relevant matrix (...vector really)
		angleTypeToKey = {"abs_roll":"water_rotations_roll_matrix", "roll":"water_rotations_roll_matrix", "pitch":"water_rotations_pitch_matrix",
		                  "azimuth":"water_rotations_azimuthal_matrix"}
		relAngleKey = angleTypeToKey[self.angleType]
		relAngles = sparseMatrixCalculator.outDict[relAngleKey]

		#Get bin-values from this matrix
		outVals = [relAngles[idx] for idx in self.oxyIndices]

		if self.angleType=="abs_roll":
			outVals = [abs(x) for x in outVals]

		return outVals

class _WaterMinDist_plusMinDistFilter_binValGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, oxyIndices, hyIndices, toIndices, filterToIndices, filterDists, minDistType):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			toIndices: (iter of ints) The indices of atoms we calculate the minimum distance TO
			filterIndices: (iter of ints) Indices of atoms we calculate minDist(toIndices[idxA]) from
			filterDists: (len-2 iter) [minDist, maxDist] for us to consider 

		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.toIndices = toIndices
		self.filterToIndices = filterToIndices
		self.filterDists = filterDists
		self.minDistType = minDistType

	def getValsToBin(self, sparseMatrixCalculator):
		#1) Filter toIndices accordingly
		toIndicesBinValGetter = _MinDistsGetOneDimValsToBin(self.toIndices, self.filterToIndices)
		toIndicesMinDists = toIndicesBinValGetter.getValsToBin(sparseMatrixCalculator)

		filteredIndices = list()
		useFilterDists = sorted(self.filterDists)
		for idx, minDist in it.zip_longest(self.toIndices, toIndicesMinDists):
			if (minDist>=useFilterDists[0]) and (minDist<useFilterDists[1]):
				filteredIndices.append( idx )

		#1.5) Check we have SOME indices; else raise an error
		if len(filteredIndices)==0:
			raise NotImplementedError("Currently no sensible behaviour when indices are filtered down to an empty list")

		#2) Simply use the non-filtering minDistBinValGetter with the filtered indices
		outBinValGetter = _WaterMinDistBinValGetter(self.oxyIndices, self.hyIndices, filteredIndices, self.minDistType)
		return outBinValGetter.getValsToBin(sparseMatrixCalculator)


class _WaterMinDistBinValGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, oxyIndices, hyIndices, toIndices, minDistType, minVal=-0.01):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			toIndices: (iter of ints) The indices of atoms we calculate the minimum distance TO
			minDistType: (str) Controls which atoms to get the minimum distance from. Current options are "all","o", and "h" (case insensitive)		
			minVal: (float) The minimum value to take into account; this is to stop getting 0 when comparing two sets of overlapping groups
 
		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.toIndices = toIndices
		self.minDistType = minDistType
		self.minVal = minVal

	def getValsToBin(self, sparseMatrixCalculator):
		#Min values on a PER WATER basis (so up to 3-minimum values per water)	
		currBinners = self._getBinners()
		outValsAll = [x.getValsToBin(sparseMatrixCalculator) for x in currBinners]

		#Map to a single value per water; in the case of minDistAll this means min(  min(O-X), min(Ha-X), min(Hb-X) )
		outVals = list()
		for vals in it.zip_longest(*outValsAll):
			outVals.append( min(vals) )

		return outVals

	def _getBinners(self):
		outBinners = list()
		relFromIndices = self._getFromIndices()

		#Get min dist for each relevant atom in each water
		for indices in relFromIndices:
			currBinner = _MinDistsGetOneDimValsToBin(indices, self.toIndices, minVal=self.minVal)
			outBinners.append(currBinner)

		return outBinners


	def _getFromIndices(self):
		oxyIndices = self.oxyIndices
		hyIndicesA = [x[0] for x in self.hyIndices]
		hyIndicesB = [x[1] for x in self.hyIndices]

		if self.minDistType.upper()=="ALL":
			outIndices = [oxyIndices, hyIndicesA, hyIndicesB]
		elif self.minDistType.upper()=="O":
			outIndices = [oxyIndices]
		elif self.minDistType.upper()=="H":
			outIndices = [hyIndicesA, hyIndicesB]
		else:
			raise ValueError("{} is an invalid value for self.minDistType".format(self.minDistType))

		return outIndices



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





class _CountHBondsBetweenGenericGroupsBinValGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, fromNonHyIndices, fromHyIndices, toNonHyIndices, toHyIndices, acceptor=True, donor=True, maxOO=3.5, maxAngle=35):
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

		"""
		self.fromNonHyIndices = fromNonHyIndices
		self.fromHyIndices = fromHyIndices
		self.toNonHyIndices = toNonHyIndices
		self.toHyIndices = toHyIndices
		self.acceptor = acceptor
		self.donor = donor
		self.maxOO = maxOO
		self.maxAngle = maxAngle

	def getValsToBin(self, sparseMatrixCalculator):
		distMatrix = sparseMatrixCalculator.outDict["distMatrix"]
		angleMatrix = sparseMatrixCalculator.outDict["angleMatrix"]

		outVals = list()
		sharedKwargs = {"acceptor":self.acceptor, "donor":self.donor, "maxOO":self.maxOO, "maxAngle":self.maxAngle}
		for fromNonHy,fromHy in it.zip_longest(self.fromNonHyIndices, self.fromHyIndices):
			currArgs = [fromNonHy, fromHy, self.toNonHyIndices, self.toHyIndices, distMatrix, angleMatrix]
			currVal = _getNumberHBondsForOneGenericFromGroup(*currArgs, **sharedKwargs)
			outVals.append( currVal )

		return outVals


class _HBondsPropertiesBetweenGenericGroupsBinValGetter(_CountHBondsBetweenGenericGroupsBinValGetter):


	def getValsToBin(self, sparseMatrixCalculator):
		distMatrix = sparseMatrixCalculator.outDict["distMatrix"]
		angleMatrix = sparseMatrixCalculator.outDict["angleMatrix"]

		outVals = list()
		sharedKwargs = {"acceptor":self.acceptor, "donor":self.donor, "maxOO":self.maxOO, "maxAngle":self.maxAngle}
		for fromNonHy, fromHy in it.zip_longest(self.fromNonHyIndices, self.fromHyIndices):
			currArgs = [fromNonHy, fromHy, self.toNonHyIndices, self.toHyIndices, distMatrix, angleMatrix]
			currIndices = _getHBondIndicesForOneGenericFromGroup(*currArgs, **sharedKwargs)
			currVals = self._getValsFromIndices( distMatrix, angleMatrix, currIndices )
			outVals.extend(currVals)

		return outVals

	def _getValsFromIndices(self, distMatrix, angleMatrix, indices):
		raise NotImplementedError("")




class _GetNonHyDistsForHBondsBinValGetter(_HBondsPropertiesBetweenGenericGroupsBinValGetter):

	def _getValsFromIndices(self, distMatrix, angleMatrix, indices):
		outDists = [ distMatrix[currIndices[0]][currIndices[1]] for currIndices in indices ]
		return outDists


class _GetOOHAnglesForHBondsBinValGetter(_HBondsPropertiesBetweenGenericGroupsBinValGetter):
	
	def _getValsFromIndices(self, distMatrix, angleMatrix, indices):
		outAngles = [ angleMatrix[currIndices[0]][currIndices[1]][currIndices[2]] for currIndices in indices ]
		return outAngles


class _GetNonHyToHyDistsForHBondsBinValGetter(_HBondsPropertiesBetweenGenericGroupsBinValGetter):

	def _getValsFromIndices(self, distMatrix, angleMatrix, indices):
#		pass
		outVals = [ distMatrix[currIndices[0]][currIndices[2]] for currIndices in indices ] 
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

		#1.5) Figure out where EACH of the groupA/groupB indices are within the self.oxyIndices
		groupBOxyListIndices = [self.oxyIndices.index(idx) for idx in groupBOxyIndices]
		groupBHyIndices = [self.hyIndices[idx] for idx in groupBOxyListIndices]

		#2) Count the number of hydrogen bonds from groupA TO groupB + output in a 1-dim list
		outVals = list()
		for lIdx,oxyIdx in enumerate(self.oxyIndices):
			oxyIdx = self.oxyIndices[lIdx]
			hyIndices = self.hyIndices[lIdx]
			if oxyIdx in groupAOxyIndices:
				currVal = self._getValForOneOxyIdx( oxyIdx, hyIndices, groupBOxyIndices, groupBHyIndices, distMatrix, angleMatrix)
			else:
				currVal = 0 #We only calc groupA -> groupB; thus the index has to be in groupA to be non-zero here
			outVals.append(currVal)

		return outVals

	def _getValForOneOxyIdx(self, oxyIdx, hyIndices, toOxyIndices, toHyIndices, distMatrix, angleMatrix):
		currArgs = [oxyIdx, hyIndices, toOxyIndices, toHyIndices, distMatrix, angleMatrix]
		currKwargs = {"acceptor":self.acceptor, "donor":self.donor, "maxOO":self.maxOO, "maxAngle":self.maxAngle}
		return _getNumberHBondsForOneOxyIdx(*currArgs, **currKwargs)


#Tons of overlap with _getFullOutAngleIndicesRequired_GENERIC in the populators file
def _getNumberHBondsForOneGenericFromGroup(fromNonHyIndices, fromHyIndices, allToNonHyIndices, allToHyIndices, distMatrix, angleMatrix, acceptor=True, donor=True, maxOO=3.5, maxAngle=35):

	return len( _getHBondIndicesForOneGenericFromGroup(fromNonHyIndices, fromHyIndices, allToNonHyIndices, allToHyIndices, distMatrix, angleMatrix, acceptor=acceptor, donor=donor, maxOO=maxOO, maxAngle=maxAngle) )


#Gives values [O_A, O_D, H_D]
def _getHBondIndicesForOneGenericFromGroup(fromNonHyIndices, fromHyIndices, allToNonHyIndices, allToHyIndices, distMatrix, angleMatrix, acceptor=True, donor=True, maxOO=3.5, maxAngle=35):
	outIndices = list()

	atomComboPopulatorHelp._checkGenericNonHyAndHyIndicesHaveValidValues(fromNonHyIndices, fromHyIndices)
	for toNonHyIndices, toHyIndices in it.zip_longest(allToNonHyIndices, allToHyIndices):
		atomComboPopulatorHelp._checkGenericNonHyAndHyIndicesHaveValidValues(toNonHyIndices, toHyIndices)

		if sorted(fromNonHyIndices) != sorted(toNonHyIndices):
			if donor:
				donorPairs = [ [nonHy,hy] for nonHy, hy in it.product(fromNonHyIndices, fromHyIndices) ]
				for donorPair, acceptorIdx in it.product(donorPairs, toNonHyIndices):
					currDist = distMatrix[donorPair[0]][acceptorIdx]
					if currDist < maxOO:
						currAngle = angleMatrix[acceptorIdx][donorPair[0]][donorPair[1]]
						if currAngle < maxAngle:
							outIndices.append( [acceptorIdx,donorPair[0],donorPair[1]] )
		
			if acceptor:
				donorPairs = [ [nonHy,hy] for nonHy, hy in it.product(toNonHyIndices, toHyIndices) ]
				for donorPair, acceptorIdx in it.product(donorPairs, fromNonHyIndices):
					currDist = distMatrix[donorPair[0]][acceptorIdx]
					if currDist < maxOO:
						currAngle = angleMatrix[acceptorIdx][donorPair[0]][donorPair[1]]
						if currAngle < maxAngle:
							outIndices.append( [acceptorIdx,donorPair[0],donorPair[1]] )

	return outIndices

def _getNumberHBondsForOneOxyIdx(fromOxyIdx, fromHyIdxPair, toOxyIndices, toHyIndices, distMatrix, angleMatrix, acceptor=True, donor=True, maxOO=3.5, maxAngle=35):
	
	outVal = 0
	for toOxyIdx, toHyIdxPair in it.zip_longest(toOxyIndices, toHyIndices):
		currDist = distMatrix[fromOxyIdx][toOxyIdx]
		if currDist<maxOO and toOxyIdx!=fromOxyIdx:
			#Note we want angle [O_A, O_D, H_D]
			if acceptor:
				hyIdxA, hyIdxB = toHyIdxPair
				angleA, angleB = angleMatrix[fromOxyIdx][toOxyIdx][hyIdxA], angleMatrix[fromOxyIdx][toOxyIdx][hyIdxB]
				if angleA < maxAngle:
					outVal += 1
				if angleB < maxAngle:
					outVal += 1

			if donor:
				hyIdxA, hyIdxB = fromHyIdxPair
				angleA, angleB = angleMatrix[toOxyIdx][fromOxyIdx][hyIdxA], angleMatrix[toOxyIdx][fromOxyIdx][hyIdxB]
				if angleA < maxAngle:
					outVal += 1
				if angleB < maxAngle:
					outVal += 1

	return outVal



class _GetDiatomDistsBinvalGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, diatomIndices):
		self.diatomIndices = diatomIndices


	def getValsToBin(self, sparseMatrixCalculator):
		relevantMatrix = sparseMatrixCalculator.outDict["distMatrix"]
		outVals = list()
		for currIndices in self.diatomIndices:
			fromIdx,toIdx = currIndices
			outVals.append( relevantMatrix[fromIdx][toIdx] )

		return outVals


class _GetDiatomAngleWithVectorBinvalGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, diatomIndices, inpVector, leftToRight):
		""" Initializer
		
		Args:
			diatomIndices: (iter of len-2 int-iters) Contains indices of diatoms we want angles for
			inpVector: (len-3 iter) The vector we want the angle with
			leftToRight: (Bool) If True then for indices [a,b] we do angle for (b-a); else we do angle for (a-b)

		"""
		self.diatomIndices = diatomIndices
		self.inpVector = inpVector
		self.leftToRight = leftToRight

	def getValsToBin(self, sparseMatrixCalculator):
		#1) Get the relevant matrix
		relevantDicts = sparseMatrixCalculator.outDict["diatom_arbitrary_angle_matrices"]
		relevantIdx = None
		for idx,currDict in enumerate(relevantDicts):
			if np.allclose( np.array(self.inpVector), np.array(currDict["inpVector"]) ):
				relevantIdx = idx
				break

		#1.5) Raise if we didnt find the correct matrix
		if relevantIdx is None:
			raise ValueError("Couldnt find the relevant matrix")

		#2) Get the values
		outVals = list()
		relevantMatrix = sparseMatrixCalculator.outDict["diatom_arbitrary_angle_matrices"][relevantIdx]["matrix"]
		for currIndices in self.diatomIndices:
			if self.leftToRight:
				useIndices = currIndices
			else:
				useIndices = reversed(currIndices)

			outVals.append( relevantMatrix[tuple(useIndices)] )  

		return outVals







