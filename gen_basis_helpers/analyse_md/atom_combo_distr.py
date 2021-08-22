
import itertools as it

import numpy as np


from . import atom_combo_core as atomComboCoreHelp
from . import atom_combo_opts_obj_maps as optsObjsMapHelp
from . import binned_res as binResHelp
from . import calc_distrib_core as calcDistrCoreHelp
from . import calc_dists as calcDistsHelp
from . import calc_radial_distrib_impl as calcRadImpl
from . import water_combo_distrs as waterComboDistrHelp

from ..shared import plane_equations as planeEqnHelp


def getAtomicComboDistrBinsFromOptsObjs(inpTraj, optsObjsGroups):
	""" Gets populated NDimensionalBinnedResults objects from a trajectory and options objects
	
	Args:
		inpTraj: (TrajectoryInMemory object)
		optsObjsGroups: (iter of iter of calcOptions objects). Currently these can have "CalcRdfOptions" (with minDistAToB=True) or "" or "CalcPlanarRdfOptions" option objs (a mixture is fine). 
			 
	Returns
		outRes: (iter of NDimensionalBinnedResults) One of these per element in optsObjs
 
	Raises:
		ValueError: Should be raised if first set of indices are not all the same within one element of optsObjsGroups (e.g. all objects in optsObjsGroups[1] must be the same; but they can differ between optsObjsGroups[0] and optsObjsGroups[1]
	"""

	#1) Check options groups are consistent; each should have the same set of indicesA/indices
	for optObjGroup in optsObjsGroups:
		_checkIndicesConsistentInOptsObjGroup(optObjGroup)

	#2) Setup object to calculate basic matrices (e.g. dist matrix) and one to get bin values from this
	sparseMatrixCalculator =  optsObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter( [optObj for optObj in it.chain(*optsObjsGroups)] )

	binValGetters = list()
	for group in optsObjsGroups:
		currBinValGetter = optsObjsMapHelp.getMultiDimBinValGetterFromOptsObjs(group)
		binValGetters.append( currBinValGetter )

	#3) Setup the bin objects
	outBinObjs = [_getBinObjForOptsObjGroup(group) for group in optsObjsGroups]

	#4) Loop over trajectory + bin values
	for trajStep in inpTraj:
		currGeom = trajStep.unitCell
		sparseMatrixCalculator.calcMatricesForGeom(currGeom)
		for binValGetter,binObj in it.zip_longest(binValGetters,outBinObjs):
			currBinVals = binValGetter.getValsToBin(sparseMatrixCalculator)
			binObj.addBinValuesToCounts(currBinVals)

	#Get the normalised counts
	nSteps = len(inpTraj.trajSteps)
	for binObj in outBinObjs:
		waterComboDistrHelp._attachCountsNormalisedByNStepsToBin(binObj, nSteps)


	return outBinObjs


def _checkIndicesConsistentInOptsObjGroup(optsObjGroup):
	primaryIndices = [_getPrimaryIndicesFromOptObj(optObj) for optObj in optsObjGroup]
	if any([x!=primaryIndices[0] for x in primaryIndices]):
		raise ValueError("Primary indices not consistent within optsObjGroup")


def _getPrimaryIndicesFromOptObj(optObj):
	try:
		outIndices = optObj.indicesA
	except AttributeError:
		try:
			outIndices = optObj.indices
		except AttributeError:
			outIndices = optObj.primaryIndices

	return outIndices


def _getBinObjForOptsObjGroup(optsObjGroup):
	binEdges = list()

	#We want to allow the bin-results from a single options object to be an iterable
	#Original use case was for allowing this interface to work for filtered options objects
	for optsObj in optsObjGroup:
		try:
			currBinObjs = iter(optsObj.binResObj)
		except TypeError:
			binEdges.append( optsObj.binResObj.binEdges )
		else:
			currEdges = [x.binEdges for x in optsObj.binResObj]
			binEdges.extend( currEdges )

	outBinObj = binResHelp.NDimensionalBinnedResults(binEdges)
	outBinObj.initialiseCountsMatrix()
	return outBinObj




