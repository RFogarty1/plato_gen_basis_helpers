
""" Specific implementations for calculating combined distributions """

import itertools as it
import numpy as np

from . import binned_res as binResHelp
from . import water_rotations as waterRotHelp


def getMultipleCombinedWaterRotationDistribBinsFromOptObjs(inpTraj ,optsObjsGroups):
	""" Similar to getCombinedWaterRotationDistribBinsFromOptsObjs but works on multiple iters of optsObj simultaneously
	
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


def _getOutBinObjForCombinedWaterRotations(optsObjs):
	#Check all water indices are the same
	sortedWaterIndices = [sorted(x.waterIndices) for x in optsObjs]
	if not all([x==sortedWaterIndices[0] for x in sortedWaterIndices]):
		raise ValueError("Water indices must all be the same")

	#Get bins + check that all waterIndices are the same
	binEdges = [x.binResObj.binEdges for x in optsObjs]
	waterIndices = optsObjs[0].waterIndices
	angleTypeToIdx = {"roll":0, "pitch":1, "azimuth":2}
	angleIndices = [angleTypeToIdx[x.angleType] for x in optsObjs]
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


