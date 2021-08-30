
import copy
import itertools as it

from . import atom_combo_populators as atomComboPopulatorHelp
from . import atom_combo_binval_getters as binValGettersHelp
from . import atom_combo_core as coreComboHelp
from . import calc_distrib_core as calcDistrCoreHelp
from . import calc_radial_distrib_impl as calcRadialDistrImplHelp
from . import distr_opt_objs as distrOptsObjHelp
from . import classification_distr_opt_objs as classDistrOptObjHelp
from . import classification_binval_getters as classBinvalGetterHelp
from . import filtered_atom_combo_opt_objs as filteredAtomComboOptHelp
from . import filtered_atom_combo_binval_getters as filteredAtomBinvalGetterHelp

from ..shared import register_key_decorator as regKeyDecoHelp
from ..shared import plane_equations as planeEqnHelp

#NOTE:Registration imports occur at the bottom of this file
def _importRegistrationModules():
	from . import filtered_atom_combo_obj_maps


#Globals here + Decorators here
_TYPE_TO_POPULATOR_DICT = dict()
_TYPE_TO_BINNER_DICT = dict()

TYPE_TO_POPULATOR_REGISTER_DECO = regKeyDecoHelp.RegisterKeyValDecorator(_TYPE_TO_POPULATOR_DICT)
TYPE_TO_BINNER_REGISTER_DECO = regKeyDecoHelp.RegisterKeyValDecorator(_TYPE_TO_BINNER_DICT)


#Main functions here
def getSparseMatrixCalculatorFromOptsObjIter(optsObjIter):
	""" Function maps an iter of opts objects into a single _SparseMatrixCalculatorStandard object. This should handle calculation of any sparse matrices required in an efficient way
	
	Args:
		optsObjIter: (iter of CalcDistribOptionsBase subclasses)
			 
	Returns
		sparseMatrixCalculator:  (_SparseMatrixCalculatorStandard object)
 
	"""
	populators = [getMatrixPopulatorFromOptsObj(optsObj) for optsObj in optsObjIter]
	return coreComboHelp._SparseMatrixCalculatorStandard(populators)


def getMatrixPopulatorFromOptsObj(optsObj):
	""" Function to get a _SparseMatrixPopulator from an CalcDistribOptionsBase-derived instance. Works based on type-casting, with types->map function being stored in a module-level dictionary
	
	Args:
		optsObj: (subclass of CalcDistribOptionsBase) Contains options for getting some kind of dsitribution function
			 
	Returns
		matrixPopulator: (_SparseMatrixPopulator)
 
	"""
	return _TYPE_TO_POPULATOR_DICT[type(optsObj)](optsObj)


def getMultiDimBinValGetterFromOptsObjs(optsObjIter):
	""" Function maps an iter of opts objects into a single _GetMultiDimValsToBinFromSparseMatrices from a CalcDistribOptionsBase-derived instance. Works based on type-casting, with types->map function being stored in a module-level dictionary
	
	Args:
		optsObjIter: (iter of CalcDistribOptionsBase subclasses)
			 
	Returns
		multiDimBinValsGetter: (_GetMultiDimValsToBinFromSparseMatrices) 
 
	"""
	singleObjs = [getOneDimBinValGetterFromOptsObj(optObj) for optObj in optsObjIter]

	#NOTE: Some options objects wll return iterators of bin-val getters; need to check and handle this
	outObjs = list()
	for obj in singleObjs:
		try:
			iter(obj)
		except TypeError:
			outObjs.append(obj)
		else:
			outObjs.extend(obj)


	return coreComboHelp._GetMultiDimValsToBinFromSparseMatrices(outObjs)


def getOneDimBinValGetterFromOptsObj(optsObj):
	""" Gets a _GetOneDimValsToBinFromSparseMatricesBase object from an CalcDistribOptionsBase-derived instance. This is an object which maps a set of matrices (e.g. a distance matrix) to a one-dimensional set of values to bin (e.g. minimum distance between two atoms for a set of input indices)
	
	Args:
		optsObj: (subclass of CalcDistribOptionsBase) Contains options for getting some kind of dsitribution function
			 
	Returns
		binValGetter: (_GetOneDimValsToBinFromSparseMatricesBase)

	NOTES:
		a) The output can instead be an iterator of binValGetter objects. This can be useful where an options object naturally maps to a multi-dimensional bin structure.
 
	"""
	return _TYPE_TO_BINNER_DICT[type(optsObj)](optsObj)




#Registration of some standard options below
@TYPE_TO_POPULATOR_REGISTER_DECO(distrOptsObjHelp.CalcRdfOptions)
@TYPE_TO_POPULATOR_REGISTER_DECO(calcDistrCoreHelp.CalcRdfOptions)
def _(inpObj):
	return atomComboPopulatorHelp._DistMatrixPopulator(inpObj.indicesA, inpObj.indicesB)

@TYPE_TO_POPULATOR_REGISTER_DECO(distrOptsObjHelp.CalcPlanarDistOptions)
@TYPE_TO_POPULATOR_REGISTER_DECO(calcRadialDistrImplHelp.CalcPlanarRdfOptions)
def _(inpObj):
	planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0) if inpObj.planeEqn is None else inpObj.planeEqn
	return atomComboPopulatorHelp._PlanarDistMatrixPopulator(inpObj.indices, planeEqn)


@TYPE_TO_POPULATOR_REGISTER_DECO(distrOptsObjHelp.DiscHBondCounterWithOxyDistFilterOptions)
def _(inpObj):
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices, inpObj.distFilterIndices, inpObj.distFilterVals]
	currKwargs = {"acceptor":inpObj.acceptor, "donor":inpObj.donor, "maxOO":inpObj.maxOO}
	return atomComboPopulatorHelp._DiscHBondCounterBetweenGroupsWithOxyDistFilterPopulator(*currArgs, **currKwargs)


@TYPE_TO_POPULATOR_REGISTER_DECO(distrOptsObjHelp.CountHBondsBetweenWaterGroupsOptions)
def _(inpObj):
	currArgs = [inpObj.fromOxyIndices, inpObj.fromHyIndices, inpObj.toOxyIndices, inpObj.toHyIndices]
	currKwargs = {"acceptor":inpObj.acceptor, "donor":inpObj.donor, "maxOO":inpObj.maxOO}
	return atomComboPopulatorHelp._CountHBondsBetweenGroupsPopulator(*currArgs, **currKwargs)


@TYPE_TO_POPULATOR_REGISTER_DECO(distrOptsObjHelp.WaterPlanarDistOptions)
def _(inpObj):
	planeEqn =  _getDefaultPlaneEquation() if inpObj.planeEqn is None else inpObj.planeEqn
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices, planeEqn]
	currKwargs = {"primaryIdxType":inpObj.primaryIdxType}
	return atomComboPopulatorHelp._WaterPlanarDistPopulator(*currArgs, **currKwargs)

@TYPE_TO_POPULATOR_REGISTER_DECO(distrOptsObjHelp.WaterMinPlanarDistOptions)
def _(inpObj):
	planeEqn = _getDefaultPlaneEquation() if inpObj.planeEqn is None else inpObj.planeEqn
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices, planeEqn]
	currKwargs = {"primaryOnly":False, "primaryIdxType":inpObj.primaryIdxType} #PrimaryIdxType basically irrelevant here i think...
	return atomComboPopulatorHelp._WaterPlanarDistPopulator(*currArgs, **currKwargs)

@TYPE_TO_POPULATOR_REGISTER_DECO(distrOptsObjHelp.WaterMinDistOptions)
def _(inpObj):
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices, inpObj.toIndices, inpObj.minDistType]
	return atomComboPopulatorHelp._WaterMinDistPopulator(*currArgs)

@TYPE_TO_POPULATOR_REGISTER_DECO(distrOptsObjHelp.WaterMinDistPlusMinDistFilterOptions)
def _(inpObj):
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices, inpObj.toIndices, inpObj.filterToIndices, inpObj.minDistType]
	return atomComboPopulatorHelp._WaterMinDist_plusMinDistFilter_populator(*currArgs)

@TYPE_TO_POPULATOR_REGISTER_DECO(distrOptsObjHelp.WaterOrientationOptions)
def _(inpObj):
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices]
	return atomComboPopulatorHelp._WaterOrientationPopulator(*currArgs)

@TYPE_TO_POPULATOR_REGISTER_DECO(distrOptsObjHelp.CalcHozDistOptions)
def _(inpObj):
	currArgs = [inpObj.indicesFrom, inpObj.indicesTo]
	return atomComboPopulatorHelp._HozDistMatrixPopulator(*currArgs)

@TYPE_TO_POPULATOR_REGISTER_DECO(classDistrOptObjHelp.WaterCountTypesMinDistAndHBondSimpleOpts)
def _(inpObj):
	#Need a DiscHBondCounter (with ridic filter things?) and a distFilter counter i guess
	#a) Create the DiscHBond counter populators
	distFilterVals = [ [0,1000], [0,1000] ] #Populator shouldnt need to know about distFilterVals really; since we count H-bonds between EVERY water
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices, inpObj.distFilterIndices, distFilterVals] 
	currKwargs = {"maxOO":inpObj.maxOOHBond, "donor":True, "acceptor":True}
	return atomComboPopulatorHelp._DiscHBondCounterBetweenGroupsWithOxyDistFilterPopulator(*currArgs, **currKwargs)


@TYPE_TO_POPULATOR_REGISTER_DECO(classDistrOptObjHelp.WaterAdsorbedClassifier_usingMinHozDistsBetweenAdsorptionSitesOptsObj)
def _(inpObj):
	#1) Get the populator relevant to the hydrogen-bonding stuff (and min dist to adsorbate sites)
	distFilterVals = [ [0,1000], [0,1000] ] #Populator shouldnt need to know about distFilterVals really; since we count H-bonds between EVERY water
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices, inpObj.distFilterIndices, distFilterVals] 
	currKwargs = {"maxOO":inpObj.maxOOHBond, "donor":True, "acceptor":True}
	hBondPopulator = atomComboPopulatorHelp._DiscHBondCounterBetweenGroupsWithOxyDistFilterPopulator(*currArgs, **currKwargs)

	#2) Get the populator relevant to finding horizontal distances between adsorption sites
	indicesFrom, indicesTo = inpObj.distFilterIndices, inpObj.distFilterIndices
	hozDistPopulator = atomComboPopulatorHelp._HozDistMatrixPopulator(indicesFrom, indicesTo)

	compositePopulator = coreComboHelp._SparseMatrixPopulatorComposite([hBondPopulator, hozDistPopulator])

	return compositePopulator


@TYPE_TO_POPULATOR_REGISTER_DECO(classDistrOptObjHelp.AtomClassifyBasedOnDistsFromIndicesSimpleOpts)
def _(inpObj):
	fromIndices, toIndices = inpObj.atomIndices, inpObj.distFilterIndices
	return atomComboPopulatorHelp._DistMatrixPopulator(fromIndices, toIndices)



#Registration of standard binners below
@TYPE_TO_BINNER_REGISTER_DECO(distrOptsObjHelp.CalcRdfOptions)
@TYPE_TO_BINNER_REGISTER_DECO(calcDistrCoreHelp.CalcRdfOptions)
def _(inpObj):
	fromIndices, toIndices = inpObj.indicesA, inpObj.indicesB
	if inpObj.minDistAToB is False:
		if inpObj.filterBasedOnBins is True:
			minVal, maxVal = min(inpObj.binResObj.binEdges), max(inpObj.binResObj.binEdges)
			outObj = binValGettersHelp._RadialDistsGetValsToBin(fromIndices, toIndices, minVal=minVal, maxVal=maxVal)
		else:
			outObj = binValGettersHelp._RadialDistsGetValsToBin(fromIndices, toIndices)
	else:
		outObj = binValGettersHelp._MinDistsGetOneDimValsToBin(fromIndices, toIndices)

	return outObj

@TYPE_TO_BINNER_REGISTER_DECO(distrOptsObjHelp.CalcHozDistOptions)
def _(inpObj):
	fromIndices, toIndices = inpObj.indicesFrom, inpObj.indicesTo
	if inpObj.minDistAToB is False:
		outObj = binValGettersHelp._HozDistsGetValsToBin(fromIndices, toIndices)
	else:
		outObj = binValGettersHelp._MinHozDistsGetValsToBin(fromIndices, toIndices, minVal=inpObj.minDistVal)

	return outObj

@TYPE_TO_BINNER_REGISTER_DECO(distrOptsObjHelp.CalcPlanarDistOptions)
@TYPE_TO_BINNER_REGISTER_DECO(calcRadialDistrImplHelp.CalcPlanarRdfOptions)
def _(inpObj):
	planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0) if inpObj.planeEqn is None else inpObj.planeEqn
	return binValGettersHelp._PlanarDistsGetOneDimValsToBin(planeEqn, inpObj.indices)


@TYPE_TO_BINNER_REGISTER_DECO(distrOptsObjHelp.CountHBondsBetweenWaterGroupsOptions)
def _(inpObj):
	currArgs = [inpObj.fromOxyIndices, inpObj.fromHyIndices, inpObj.toOxyIndices, inpObj.toHyIndices]
	currKwargs = {"acceptor":inpObj.acceptor, "donor":inpObj.donor, "maxOO":inpObj.maxOO, "maxAngle":inpObj.maxAngle}
	return binValGettersHelp._CountHBondsBetweenWaterGroupsBinValGetter(*currArgs,**currKwargs)


@TYPE_TO_BINNER_REGISTER_DECO(distrOptsObjHelp.DiscHBondCounterWithOxyDistFilterOptions)
def _(inpObj):
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices, inpObj.distFilterIndices, inpObj.distFilterVals]
	currKwargs = {"acceptor":inpObj.acceptor, "donor":inpObj.donor, "maxOO":inpObj.maxOO, "maxAngle":inpObj.maxAngle}
	return binValGettersHelp._DiscHBondCounterBetweenGroupsWithOxyDistFilterOneDimValGetter(*currArgs, **currKwargs)


@TYPE_TO_BINNER_REGISTER_DECO(distrOptsObjHelp.WaterPlanarDistOptions)
def _(inpObj):
	planeEqn =  _getDefaultPlaneEquation() if inpObj.planeEqn is None else inpObj.planeEqn
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices, planeEqn]
	currKwargs = {"primaryIdxType":inpObj.primaryIdxType}
	return binValGettersHelp._WaterPlanarDistBinValGetter(*currArgs, **currKwargs)

@TYPE_TO_BINNER_REGISTER_DECO(distrOptsObjHelp.WaterMinPlanarDistOptions)
def _(inpObj):
	planeEqn =  _getDefaultPlaneEquation() if inpObj.planeEqn is None else inpObj.planeEqn
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices, planeEqn]
	currKwargs = {"minDistType":inpObj.minDistType}
	return binValGettersHelp._WaterPlanarMinDistBinValGetter(*currArgs, **currKwargs)

@TYPE_TO_BINNER_REGISTER_DECO(distrOptsObjHelp.WaterMinDistOptions)
def _(inpObj):
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices, inpObj.toIndices, inpObj.minDistType]
	currKwargs = {"minVal":inpObj.minVal}
	return binValGettersHelp._WaterMinDistBinValGetter(*currArgs, **currKwargs)

@TYPE_TO_BINNER_REGISTER_DECO(distrOptsObjHelp.WaterMinDistPlusMinDistFilterOptions)
def _(inpObj):
	currArgs = [inpObj.oxyIndices, inpObj.hyIndices, inpObj.toIndices, inpObj.filterToIndices, inpObj.filterDists, inpObj.minDistType]
	return binValGettersHelp._WaterMinDist_plusMinDistFilter_binValGetter(*currArgs)

@TYPE_TO_BINNER_REGISTER_DECO(distrOptsObjHelp.WaterOrientationOptions)
def _(inpObj):
	currArgs = [inpObj.oxyIndices, inpObj.angleType]
	return binValGettersHelp._WaterOrientationBinValGetter(*currArgs)

@TYPE_TO_BINNER_REGISTER_DECO(classDistrOptObjHelp.AtomClassifyBasedOnDistsFromIndicesSimpleOpts)
def _(inpObj):
	outObjs = list()
	for idx,unused in enumerate(inpObj.distFilterRanges):
		currArgs = [inpObj.atomIndices, inpObj.distFilterIndices, inpObj.distFilterRanges[idx]]
		currKwargs = {"minDistVal":inpObj.minDistVal}
		currObj = classBinvalGetterHelp._AtomsWithMinDistRangeCountBinvalGetter(*currArgs, **currKwargs)
		outObjs.append(currObj)
	return outObjs

@TYPE_TO_BINNER_REGISTER_DECO(classDistrOptObjHelp.WaterCountTypesMinDistAndHBondSimpleOpts)
def _(inpObj):
	outObjs = list()
	for idx,unused in enumerate(inpObj.distFilterRanges):
		currArgs = [inpObj.oxyIndices, inpObj.hyIndices, inpObj.distFilterIndices, inpObj.distFilterRanges[idx],
		            inpObj.nDonorFilterRanges[idx], inpObj.nAcceptorFilterRanges[idx], inpObj.nTotalFilterRanges[idx],
		            inpObj.maxOOHBond, inpObj.maxAngleHBond]
		currObj = classBinvalGetterHelp._WaterCountTypeBinvalGetter(*currArgs)
		outObjs.append(currObj)
	return outObjs

@TYPE_TO_BINNER_REGISTER_DECO(classDistrOptObjHelp.WaterAdsorbedClassifier_usingMinHozDistsBetweenAdsorptionSitesOptsObj)
def _(inpObj):
	outObjs = list()
	for idx,unused in enumerate(inpObj.distFilterRanges):
		currArgs = [inpObj.oxyIndices, inpObj.hyIndices, inpObj.distFilterIndices, inpObj.distFilterRanges[idx],
		            inpObj.nDonorFilterRanges[idx], inpObj.nAcceptorFilterRanges[idx], inpObj.nTotalFilterRanges[idx],
		            inpObj.maxOOHBond, inpObj.maxAngleHBond, inpObj.adsSiteMinHozToOtherAdsSiteRanges[idx]]
		currObj = classBinvalGetterHelp._AdsorbedWaterCountTypeWithAdsSiteHozDistsBinvalGetter(*currArgs)
		outObjs.append(currObj)
	return outObjs

#Utility functions
def _getDefaultPlaneEquation():
	return planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0)




#Imports
_importRegistrationModules()
