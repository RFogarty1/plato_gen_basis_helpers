
""" Specific implementations for calculating combined distributions """

import itertools as it
import numpy as np

from . import binned_res as binResHelp
from . import calc_distrib_core as calcDistribCoreHelp
from . import water_rotations as waterRotHelp
from . import calc_dists as calcDistHelp


def getMultipleWaterComboDistribBinsFromOptObjs(inpTraj, optsObjsGroups):
	""" Gets populated NDimensionalBinnedResults objects from a trajectory and options objects
	
	Args:
		inpTraj: (TrajectoryInMemory object)
		optsObjsGroups: (iter of iter of calcOptions objects). Currently these can have "CalcStandardWaterOrientationDistribOptions" or "CalcWaterPlanarDistribOptions_fromOxy" option objs (a mixture is fine). 
			 
	Returns
		outRes: (iter of NDimensionalBinnedResults) One of these per element in optsObjs
 
	Raises:
		ValueError: Should be raised if .waterIndices are not all the same within one element of optsObjsGroups (e.g. all objects in optsObjsGroups[1] must be the same; but they can differ between optsObjsGroups[0] and optsObjsGroups[1]
	"""


	#Figure out if we need to calculate angular AND planar parameters or just one
	#Note if we need angles/planar for ANY ENTRY in optsObjsGroups we calculate for all
	calcAngles, calcPlanar = True, True
	if all( [isinstance(optObj, waterRotHelp.CalcStandardWaterOrientationDistribOptions) for optObj in it.chain(*optsObjsGroups)] ):
		calcPlanar = False
	elif all( [isinstance(optObj, CalcWaterPlanarDistribOptions_fromOxy) for optObj in it.chain(*optsObjsGroups)] ):
		calcAngles = False

	#TODO: All options objs have to have the same planeEqn for now; enforce this here
	if calcPlanar:
		allOptsObjs = [x for x in it.chain(*optsObjsGroups)]
		allPlaneEqns = [x.planeEqn for x in allOptsObjs if isinstance(x,CalcWaterPlanarDistribOptions_fromOxy)]
		assert all([x==allPlaneEqns[0] for x in allPlaneEqns])
		planeEqn = allPlaneEqns[0]

	#Get bin objects (which also checks that each element in optsObjs uses a single set of waterIndices)
	outBinObjs = [_getOutBinObjForCombinedWaterRotations(optsObjGroup) for optsObjGroup in optsObjsGroups]
	groupIndices = [group[0].waterIndices for group in optsObjsGroups]

#	#Figure out ALL the water indices needed [duplicated code from waterRotHelp]
	allIndices = [idxList for idxList in it.chain(*[indices for indices in groupIndices])]
	uniqueIndices = [ list(a) for a in set([tuple(x) for x in allIndices]) ]	
	indicesToUniqueIndicesMap = waterRotHelp._getMapFromIterToUniqueVals(allIndices, uniqueIndices)

	#Bin Values
	for trajStep in inpTraj:
		#1) Calculate anything required
		if calcAngles:
			allAngles = waterRotHelp.getWaterStandardRotationAnglesForInpCell(trajStep.unitCell, uniqueIndices)
		else:
			allAngles = None

		if calcPlanar:
			allPlanarDists = _getPlanarDistsForWaterIndices_usingOxyPos(trajStep.unitCell, uniqueIndices, planeEqn)
		else:
			allPlanarDists = None

		totalIdx = 0 #
		#Bin each 
		for groupIdx in range(len(optsObjsGroups)):
			nWater = len(groupIndices[groupIdx])
			nDims = len(optsObjsGroups[groupIdx])
			currOptsObjs = optsObjsGroups[groupIdx]

			#0) Initialize empty binVals
			currBinVals = [  [None for x in range(nDims)] for water in range(nWater)]
			relevantIndices = [ indicesToUniqueIndicesMap[idx] for idx in range(totalIdx, totalIdx+nWater) ]

			#1) Populate with required angles/planar distances
			_populateBinValsFromAngles(currOptsObjs, allAngles, relevantIndices, currBinVals)
			_populateBinValsFromPlanarDists(currOptsObjs, allPlanarDists, relevantIndices, currBinVals)

			#1) Get ALL angles and planar distances required

			totalIdx += nWater
			outBinObjs[groupIdx].addBinValuesToCounts(currBinVals)

	#Get the normalised counts
	nSteps = len(inpTraj.trajSteps)
	for binObj in outBinObjs:
		_attachCountsNormalisedByNStepsToBin(binObj, nSteps)

	return outBinObjs


def _getPlanarDistsForWaterIndices_usingOxyPos(inpGeom, waterIndices, planeEqn):
	#1) Get oxygen indices from the water indices
	useOxyIndices = list()
	fractCoords = inpGeom.fractCoords
	
	for waterIdxList in waterIndices:
		currOxyIndices = [idx for idx in waterIdxList if fractCoords[idx][-1].upper()=="O"]
		hIndices =   [idx for idx in waterIdxList if fractCoords[idx][-1].upper()=="H"]
		assert len(currOxyIndices)==1
		useOxyIndices.append(currOxyIndices[0])

	allDists = calcDistHelp.calcDistancesFromSurfPlaneForCell(inpGeom, indices=useOxyIndices, planeEqn=planeEqn)
	return allDists

def _populateBinValsFromAngles(optsObjs, allAngles, relevantIndices, inpBinValArray):
	if allAngles is None:
		return None


	angleTypeToIdx = {"roll":0, "pitch":1, "azimuth":2}

	for outIdx,relIdx in enumerate(relevantIndices):
		for optIdx,optObj in enumerate(optsObjs):
			if isinstance(optObj, waterRotHelp.CalcStandardWaterOrientationDistribOptions):
				angleIdx = angleTypeToIdx[optObj.angleType]
				inpBinValArray[outIdx][optIdx] = allAngles[relIdx][angleIdx]


def _populateBinValsFromPlanarDists(optsObjs, allDists, relevantIndices, inpBinValArray):
	if allDists is None:
		return None

	#Very similar to above
	for outIdx,relIdx in enumerate(relevantIndices):
		for optIdx, optObj in enumerate(optsObjs):
			if isinstance(optObj, CalcWaterPlanarDistribOptions_fromOxy):
				inpBinValArray[outIdx][optIdx] = allDists[relIdx]


#Code dealing with pure rotational distributions
#TODO: This function WILL be deprecated very soon
def getMultipleCombinedWaterRotationDistribBinsFromOptObjs(inpTraj ,optsObjsGroups):
	""" DEPRECATED. Similar to getCombinedWaterRotationDistribBinsFromOptsObjs but works on multiple iters of optsObj simultaneously
	
	Args:
		inpTraj: (TrajectoryInMemory object)
		optsObjsGroups: (iter of iter of CalcStandardWaterOrientationDistribOptions objects) Note: Within each iterable each object is expected to define a dimension; Generally expecting only two or three at most (3 likely only useful as a filter function; e.g. to restrict values within a 2-dimensional distribution)

	Returns
		outRes: (iter of NDimensionalBinnedResults) One of these per element in optsObjs
 
	"""
	#Define angle indices for each group
	angleTypeToIdx = {"roll":0, "pitch":1, "azimuth":2}
	groupAngleIndices = list()
	for group in optsObjsGroups:
		angleIndices = [angleTypeToIdx[x.angleType] for x in group]
		groupAngleIndices.append(angleIndices)

	#Get bin objects (which also checks that each element in optsObjs uses a single set of waterIndices)
	outBinObjs = [_getOutBinObjForCombinedWaterRotations(optsObjGroup) for optsObjGroup in optsObjsGroups]
	groupIndices = [group[0].waterIndices for group in optsObjsGroups]

#	#Figure out ALL the water indices needed [duplicated code from waterRotHelp]
	allIndices = [idxList for idxList in it.chain(*[indices for indices in groupIndices])]
	uniqueIndices = [ list(a) for a in set([tuple(x) for x in allIndices]) ]	
	indicesToUniqueIndicesMap = waterRotHelp._getMapFromIterToUniqueVals(allIndices, uniqueIndices)

	#Bin Values
	for trajStep in inpTraj:
		allAngles = waterRotHelp.getWaterStandardRotationAnglesForInpCell(trajStep.unitCell, uniqueIndices)
		totalIdx = 0 #
		#Bin each 
		for groupIdx in range(len(optsObjsGroups)):
			nAngles = len(groupIndices[groupIdx])
#			import pdb
#			pdb.set_trace()
			relevantAngles = [ allAngles[indicesToUniqueIndicesMap[idx]] for idx in range(totalIdx,totalIdx+nAngles) ]
			totalIdx += nAngles
			currAngleIndices = groupAngleIndices[groupIdx]
			valsToBin = _getWaterAnglesToBinFromAllRelevantAngles(relevantAngles, currAngleIndices)
			outBinObjs[groupIdx].addBinValuesToCounts(valsToBin)

	#Get the normalised counts
	nSteps = len(inpTraj.trajSteps)
	for binObj in outBinObjs:
		_attachCountsNormalisedByNStepsToBin(binObj, nSteps)

	return outBinObjs

def getCombinedWaterRotationDistribBinsFromOptsObjs(inpTraj, optsObjs):
	""" Gets a bin object containing the distribution of combined distribs in optsObjs
	
	Args:
		inpTraj: (TrajectoryInMemory object)
		optsObjs: (iter of CalcStandardWaterOrientationDistribOptions objects) Note: Each object is expected to define a dimension; Generally expecting only two or three at most (3 likely only useful as a filter function; e.g. to restrict values within a 2-dimensional distribution)
			 
	Returns
		outRes: (NDimensionalBinnedResults) .binVals["counts"] contains the counts for each in an array/tensor, while .binCentresArray contains the centres of each bin in the same ordering. Note we also have a .binVals["normalised_counts"] that divides by the number of steps (so gives average count per step)

	Raises:
		ValueError: If waterIndices varies between optsObjs.
 
	"""

	#Get bins + water indices (_getOutBinObjForCombinedWaterRotations checks all indices are the same)
	outBinObj = _getOutBinObjForCombinedWaterRotations(optsObjs)
	uniqueWaterIndices = optsObjs[0].waterIndices
	angleTypeToIdx = {"roll":0, "pitch":1, "azimuth":2}
	angleIndices = [angleTypeToIdx[x.angleType] for x in optsObjs]

	#Start binning
	for trajStep in inpTraj:
		allAngles = waterRotHelp.getWaterStandardRotationAnglesForInpCell(trajStep.unitCell, uniqueWaterIndices)
		valsToBin = _getWaterAnglesToBinFromAllRelevantAngles(allAngles, angleIndices)
		outBinObj.addBinValuesToCounts(valsToBin)

	#Get the normalised counts
	nSteps = len(inpTraj.trajSteps)
	_attachCountsNormalisedByNStepsToBin(outBinObj, nSteps)

	return outBinObj

class CalcWaterPlanarDistribOptions_fromOxy(calcDistribCoreHelp.CalcDistribOptionsBase):
	""" Contains options for calculating water-planar distribution using only the oxygen atom position """

	def __init__(self, binResObj, waterIndices, planeEqn=None, volume=None):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			waterIndices: (iter of len-3 ints) Contains the iter of indices for water molecules
			planeEqn: (None or ThreeDimPlaneEquation) The plane equation to calculate rdf from
			volume: (None or float) The volume to use for calculating the rdf. None generally means use the full unit cell volume (may not be sensible for slabs/multi-phase cells)
		"""
		self.distribKey = "rdf"
		self.binResObj = binResObj
		self.waterIndices = waterIndices
		self.planeEqn = planeEqn
		self.volume = volume


def _getOutBinObjForCombinedWaterRotations(optsObjs):
	#Check all water indices are the same
	sortedWaterIndices = [sorted(x.waterIndices) for x in optsObjs]
	if not all([x==sortedWaterIndices[0] for x in sortedWaterIndices]):
		raise ValueError("Water indices must all be the same")

	#Get bins + check that all waterIndices are the same
	binEdges = [x.binResObj.binEdges for x in optsObjs]
	waterIndices = optsObjs[0].waterIndices
	outBinObj = binResHelp.NDimensionalBinnedResults(binEdges)
	outBinObj.initialiseCountsMatrix()
	return outBinObj

def _getWaterAnglesToBinFromAllRelevantAngles(relevantAngles, angleIndices):
	""" Gets a list of values to bin when given a list of relevant angles and the indices to take for each bin
	
	Args:
		relevantAngles: (iter of len-3 iters) Each contains [roll,pitch,azimuth]
		angleIndices: (len-N iter of ints) N is the number of dimensions. Each index says where to take one bin value from, so usually (always?) N should be <=3

	Returns
		outVals: (len-M iter of len-N floats) Each of the M elements corresponds to one entry in relvant angles. Each of the N entries then corresponds to an angle in angle indices
 
	"""
	outVals = list()
	for angles in relevantAngles:
		currVals = [angles[angleIdx] for angleIdx in angleIndices ]
		outVals.append(currVals)
	return outVals

def _attachCountsNormalisedByNStepsToBin(binObj, nSteps):
	binObj.initialiseCountsMatrix(countKey="normalised_counts")
	np.copyto(binObj.binVals["normalised_counts"], binObj.binVals["counts"]) #(dst,src)
	binObj.binVals["normalised_counts"] /= nSteps




