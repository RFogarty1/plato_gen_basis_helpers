

from . import atom_combo_populators as atomComboPopulatorHelp
from . import atom_combo_binval_getters as binValGettersHelp
from . import atom_combo_core as coreComboHelp
from . import calc_distrib_core as calcDistrCoreHelp
from . import calc_radial_distrib_impl as calcRadialDistrImplHelp
from . import distr_opt_objs as distrOptsObjHelp


from ..shared import register_key_decorator as regKeyDecoHelp
from ..shared import plane_equations as planeEqnHelp


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
	return coreComboHelp._GetMultiDimValsToBinFromSparseMatrices(singleObjs)


def getOneDimBinValGetterFromOptsObj(optsObj):
	""" Gets a _GetOneDimValsToBinFromSparseMatricesBase object from an CalcDistribOptionsBase-derived instance. This is an object which maps a set of matrices (e.g. a distance matrix) to a one-dimensional set of values to bin (e.g. minimum distance between two atoms for a set of input indices)
	
	Args:
		optsObj: (subclass of CalcDistribOptionsBase) Contains options for getting some kind of dsitribution function
			 
	Returns
		binValGetter: (_GetOneDimValsToBinFromSparseMatricesBase)
 
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
	return atomComboPopulatorHelp._WaterPlanarDistPopulator


#Registration of standard binners below
@TYPE_TO_BINNER_REGISTER_DECO(distrOptsObjHelp.CalcRdfOptions)
@TYPE_TO_BINNER_REGISTER_DECO(calcDistrCoreHelp.CalcRdfOptions)
def _(inpObj):
	fromIndices, toIndices = inpObj.indicesA, inpObj.indicesB
	if inpObj.minDistAToB is False:
		raise NotImplementedError("")
	return binValGettersHelp._MinDistsGetOneDimValsToBin(fromIndices, toIndices)


@TYPE_TO_BINNER_REGISTER_DECO(distrOptsObjHelp.CalcPlanarDistOptions)
@TYPE_TO_BINNER_REGISTER_DECO(calcRadialDistrImplHelp.CalcPlanarRdfOptions)
def _(inpObj):
	planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0) if inpObj.planeEqn is None else inpObj.planeEqn
	return binValGettersHelp._PlanarDistsGetOneDimValsToBin(planeEqn, inpObj.indices)

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

#Utility functions
def _getDefaultPlaneEquation():
	return planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0)

