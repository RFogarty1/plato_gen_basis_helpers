
import copy

from .. import binned_res as binResHelp
from .. import classification_distr_opt_objs as classDistrOptObjHelp


def getSolConFarClassifierOpts(hydroxylNonHyIndices, hydroxylHyIndices, waterNonHyIndices, waterHyIndices):
	""" Gets classifier options object for solConFar; these are water which

	a) Are not solAds
	b) Have >=1 hydrogen bond to solAds
	c) Oxygen is <=2 Angstrom horizontal distance from a solAds oxygen
	
	Args:
		hydroxylNonHyIndices: (iter of len-1 int iters) Contains indices for hydroxyl oxygen atoms; e.g. [ [2], [5] ]  
		hydroxylHyIndices: (iter len-1 int iters) Contains indices for hydroxyl hydrogen atoms; e.g. [ [3], [6] ]
		waterNonHyIndices: (iter of len-1 int iters) Contains indices for water oxygen atoms; e.g. [ [7], [8] ] 
		waterHyIndices: (iter of len-2 int iters) Contains indices for water hydrogen atoms; e.g. [ [9,10], [11,12] ] 

	"""
	#1) Get the solCon part
	solConOpts = _getSolConOpts(hydroxylNonHyIndices, hydroxylHyIndices, waterNonHyIndices, waterHyIndices)

	#2) Get the horizontal distance part; this involves being a little bit hacky with ".ClassifyBasedOnHBondingToDynamicGroup"
	hBondsCountBins = binResHelp.BinnedResultsStandard.fromBinEdges([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]) #Likely could pass a dud here really
	solAdsOpts = getSolAdsClassifierOpts(hydroxylNonHyIndices, hydroxylHyIndices, waterNonHyIndices, waterHyIndices)
	filterRanges = [ [0,2] ]
	currArgs = [hBondsCountBins, waterNonHyIndices, waterHyIndices, waterNonHyIndices, waterHyIndices, filterRanges]
	hozDistOpts = classDistrOptObjHelp.ClassifyNonHyAndHyBasedOnMinHozDistsToAtomGroups(*currArgs, minDistVal=0.01)
	hozDistToSolAdsOpts = classDistrOptObjHelp.ClassifyBasedOnHBondingToDynamicGroup(solAdsOpts, hozDistOpts, checkConsistent=False) #Need to disable consistency check due to us using this in a hacky way

	#3) Combine the options together
	currOpts = [ solConOpts, hozDistToSolAdsOpts ]
	solConFarFullClassifierOpts = classDistrOptObjHelp.ClassifyNonHyAndHyChainedAllCommon(currOpts)

	return solConFarFullClassifierOpts


def getSolConCloseClassifierOpts(hydroxylNonHyIndices, hydroxylHyIndices, waterNonHyIndices, waterHyIndices):
	""" Gets classifier options object for solConClose; these are water which

	a) Are not solAds
	b) Have >=1 hydrogen bond to solAds
	c) Oxygen is >2 Angstrom horizontal distance from a solAds oxygen
	
	Args:
		hydroxylNonHyIndices: (iter of len-1 int iters) Contains indices for hydroxyl oxygen atoms; e.g. [ [2], [5] ]  
		hydroxylHyIndices: (iter len-1 int iters) Contains indices for hydroxyl hydrogen atoms; e.g. [ [3], [6] ]
		waterNonHyIndices: (iter of len-1 int iters) Contains indices for water oxygen atoms; e.g. [ [7], [8] ] 
		waterHyIndices: (iter of len-2 int iters) Contains indices for water hydrogen atoms; e.g. [ [9,10], [11,12] ] 

	"""
	#1) Get the solCon part
	solConOpts = _getSolConOpts(hydroxylNonHyIndices, hydroxylHyIndices, waterNonHyIndices, waterHyIndices)

	#2) Get the horizontal distance part; this involves being a little bit hacky with ".ClassifyBasedOnHBondingToDynamicGroup"
	hBondsCountBins = binResHelp.BinnedResultsStandard.fromBinEdges([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]) #Likely could pass a dud here really
	solAdsOpts = getSolAdsClassifierOpts(hydroxylNonHyIndices, hydroxylHyIndices, waterNonHyIndices, waterHyIndices)
	filterRanges = [ [2,10] ]
	currArgs = [hBondsCountBins, waterNonHyIndices, waterHyIndices, waterNonHyIndices, waterHyIndices, filterRanges]
	hozDistOpts = classDistrOptObjHelp.ClassifyNonHyAndHyBasedOnMinHozDistsToAtomGroups(*currArgs, minDistVal=0.01)
	hozDistToSolAdsOpts = classDistrOptObjHelp.ClassifyBasedOnHBondingToDynamicGroup(solAdsOpts, hozDistOpts, checkConsistent=False) #Need to disable consistency check due to us using this in a hacky way

	#3) Combine the options together
	currOpts = [ solConOpts, hozDistToSolAdsOpts ]
	solConCloseFullClassifierOpts = classDistrOptObjHelp.ClassifyNonHyAndHyChainedAllCommon(currOpts)
	
	return solConCloseFullClassifierOpts

def _getSolConOpts(hydroxylNonHyIndices, hydroxylHyIndices, waterNonHyIndices, waterHyIndices):
	#1) Need solAds defined
	hBondsCountBins = binResHelp.BinnedResultsStandard.fromBinEdges([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]) #Likely could pass a dud here really
	solAdsOpts = getSolAdsClassifierOpts(hydroxylNonHyIndices, hydroxylHyIndices, waterNonHyIndices, waterHyIndices)

	#2) Options for h-bonding to solAds
	currArgs = [hBondsCountBins, waterNonHyIndices, waterHyIndices, waterNonHyIndices, waterHyIndices]
	currKwargs = {"nTotalFilterRanges":[ [0.1,1000] ]}
	solConHBondOptPart = classDistrOptObjHelp.ClassifyBasedOnHBondingToGroup_simple(*currArgs, **currKwargs)
	solConHBondToSolAdsPart = classDistrOptObjHelp.ClassifyBasedOnHBondingToDynamicGroup(solAdsOpts, solConHBondOptPart, mutuallyExclusive=True)

	return solConHBondToSolAdsPart


def getSolOtherClassifierOpts(hydroxylNonHyIndices, hydroxylHyIndices, waterNonHyIndices, waterHyIndices):
	""" Gets classifier option for solOther; these are water with <1 hydrogen bond to solAds (and are NOT part of the solAds group)
	
	Args:
		hydroxylNonHyIndices: (iter of len-1 int iters) Contains indices for hydroxyl oxygen atoms; e.g. [ [2], [5] ]  
		hydroxylHyIndices: (iter len-1 int iters) Contains indices for hydroxyl hydrogen atoms; e.g. [ [3], [6] ]
		waterNonHyIndices: (iter of len-1 int iters) Contains indices for water oxygen atoms; e.g. [ [7], [8] ] 
		waterHyIndices: (iter of len-2 int iters) Contains indices for water hydrogen atoms; e.g. [ [9,10], [11,12] ] 
			
	"""
	solAdsOpts = getSolAdsClassifierOpts(hydroxylNonHyIndices, hydroxylHyIndices, waterNonHyIndices, waterHyIndices)

	hBondsCountBins = binResHelp.BinnedResultsStandard.fromBinEdges([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]) #Likely could pass a dud here really
	currArgs = [hBondsCountBins, waterNonHyIndices, waterHyIndices, waterNonHyIndices, waterHyIndices]
	currKwargs = {"nTotalFilterRanges": [ [-0.1,0.1] ]}
	hBondOptsPart = classDistrOptObjHelp.ClassifyBasedOnHBondingToGroup_simple(*currArgs,**currKwargs)
	solOtherToSolAdsOpts = classDistrOptObjHelp.ClassifyBasedOnHBondingToDynamicGroup(solAdsOpts, hBondOptsPart, mutuallyExclusive=True)

	return solOtherToSolAdsOpts


def getSolAdsClassifierOpts(hydroxylNonHyIndices, hydroxylHyIndices, waterNonHyIndices, waterHyIndices):
	""" Gets classifier options object for solAds; these are water which have >= 1 hydrogen bond to hydroxyl groups
	
	Args:
		hydroxylNonHyIndices: (iter of len-1 int iters) Contains indices for hydroxyl oxygen atoms; e.g. [ [2], [5] ]  
		hydroxylHyIndices: (iter len-1 int iters) Contains indices for hydroxyl hydrogen atoms; e.g. [ [3], [6] ]
		waterNonHyIndices: (iter of len-1 int iters) Contains indices for water oxygen atoms; e.g. [ [7], [8] ] 
		waterHyIndices: (iter of len-2 int iters) Contains indices for water hydrogen atoms; e.g. [ [9,10], [11,12] ] 
			
	"""
	hBondsCountBins = binResHelp.BinnedResultsStandard.fromBinEdges([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]) #Likely could pass a dud here really

	#solAds; >= 1 h-bond with hydroxyl
	currArgs = [hBondsCountBins, waterNonHyIndices, waterHyIndices, hydroxylNonHyIndices, hydroxylHyIndices]
	currKwargs = {"nTotalFilterRanges": [ [0.1,1000] ]}
	solAdsOpts = classDistrOptObjHelp.ClassifyBasedOnHBondingToGroup_simple(*currArgs, **currKwargs)

	return solAdsOpts







