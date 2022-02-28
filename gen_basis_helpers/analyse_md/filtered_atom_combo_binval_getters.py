
import itertools as it

from . import atom_combo_core as atomComboCoreHelp
from . import atom_combo_opts_obj_maps as atomComboOptObjMapHelp
from . import atom_combo_binval_getters as atomComboBinvalGetterHelp

from ..shared import register_key_decorator as regKeyDecoHelp

_BINVAL_GETTER_TYPE_TO_ATOM_ATOM_MOD_DICT = dict()
_BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_DICT = dict()
_BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_DICT = dict()
_BINVAL_GETTER_TYPE_TO_WATER_DERIVATIVE_MOD_DICT = dict()

BINVAL_GETTER_TYPE_TO_ATOM_ATOM_MOD_REGISTER_DECO = regKeyDecoHelp.RegisterKeyValDecorator(_BINVAL_GETTER_TYPE_TO_ATOM_ATOM_MOD_DICT)
BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_REGISTER_DECO = regKeyDecoHelp.RegisterKeyValDecorator(_BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_DICT)
BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO = regKeyDecoHelp.RegisterKeyValDecorator(_BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_DICT)
BINVAL_GETTER_TYPE_TO_WATER_DERIVATIVE_MOD_DICT = regKeyDecoHelp.RegisterKeyValDecorator(_BINVAL_GETTER_TYPE_TO_WATER_DERIVATIVE_MOD_DICT)

def _modBinvalGetterWithFilterIndicesAtomAtomDict(binValGetter, groupIndices, useGroups,**kwargs):
	_BINVAL_GETTER_TYPE_TO_ATOM_ATOM_MOD_DICT[type(binValGetter)](binValGetter, groupIndices, useGroups, **kwargs)

#Note: This is specifically for getting distributions between water of different types
#Different code would be needed to get (for example) water of one type to a fixed set of indices
def _modBinvalGetterWithFilterIndicesWaterWaterDict(binValGetter, groupIndices, useGroups, **kwargs):
	_BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_DICT[type(binValGetter)](binValGetter, groupIndices, useGroups, **kwargs)
	
def _modBinvalGetterWithFilterIndicesNonHyAndHyGenericDict(binvalGetter, groupIndices, useGroups, **kwargs):
	_BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_DICT[type(binvalGetter)](binvalGetter, groupIndices, useGroups, **kwargs)

def _modBinvalGetterWithFilterIndicesForWaterDerivativeDict(binvalGetter, groupIndices, useGroups, **kwargs):
	_BINVAL_GETTER_TYPE_TO_WATER_DERIVATIVE_MOD_DICT[type(binvalGetter)](binvalGetter, groupIndices, useGroups, **kwargs)


class FilteredAtomComboBinvalGetterGeneric(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, classifiers, binValGetter, useGroups):
		""" Initializer
		
		Args:
			classifiers: (iter of _AtomsWithinMinDistRangeClassifier objects) May have a more expanded range of classifiers in future
			binValGetter: (_GetOneDimValsToBinFromSparseMatricesBase object)
			useGroups: (len-x iter) Each element represents a group to use. Generally expecing len-1 or len-2 at most
 
		"""
		self.classifiers = classifiers
		self.binValGetter = binValGetter
		self.useGroups = useGroups

	#Basically same as the WaterWater one (minus the toIdxType)
	def getValsToBin(self, sparseMatrixCalculator):
		#1) Get groups [THIS IS DOING IT WRONG]
		groupIndices = [x.classify(sparseMatrixCalculator) for x in self.classifiers]

		#2) Modify the binval getter with relevant indices
		_modBinvalGetterWithFilterIndicesAtomAtomDict(self.binValGetter, groupIndices, self.useGroups)

		#3) Just return the binval answer now its been modified
		outVals =  self.binValGetter.getValsToBin(sparseMatrixCalculator)

		return outVals

class WaterToWaterFilteredAtomComboBinvalGetterGeneric(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):
	
	def __init__(self, waterClassifiers, binValGetter, useGroups, toIdxType):
		""" Initializer
		
		Args:
			waterClassifiers: (iter of _WaterClassifierBase objects)
			binValGetter: (_GetOneDimValsToBinFromSparseMatricesBase object)
			useGroups: (len-x iter) Each element represents a group to use. Generally expecing len-1 or len-2 at most
			toIdxType: (str) Should be "O","H" or "all"
 
		"""
		self.waterClassifiers = waterClassifiers
		self.binValGetter = binValGetter
		self.useGroups = useGroups
		self.toIdxType = toIdxType

	def getValsToBin(self, sparseMatrixCalculator):
		#1) Get groups
		groupIndices = [x.classify(sparseMatrixCalculator) for x in self.waterClassifiers]

		#2) Modify the binval getter with the relevant group indices
		_modBinvalGetterWithFilterIndicesWaterWaterDict(self.binValGetter, groupIndices, self.useGroups, toIdxType=self.toIdxType)

		#3) Just return the binval answer now its been modified
		outVals =  self.binValGetter.getValsToBin(sparseMatrixCalculator)

		return outVals

#TODO: I should probably make the water-water use the same code in the backend really 
class GenericNonHyAndHyFilteredAtomComboBinvalGetter_simple(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, nonHyAndHyClassifiers, binValGetter, useGroups, useNonHyIdx=True, useIdxEach=0):
		""" Initializer
		
		Args:
			nonHyAndHyClassifiers: (iter of ClassifierBase objects) The classify method should return two lists of indices for each group; with the first containing non-hydrogen indices (usually only one) and the second containing hydrogen indices. Likely will only ever use "_GenericNonHyAndHyClassiferUsingHBondsToGroup_simple"
			binValGetter: (_GetOneDimValsToBinFromSparseMatricesBase object)
			useGroups: (len-x iter) Each element represents a group to use. Generally expecing len-1 or len-2 at most
			useNonHyIdx: (Bool) If True we represent our group with one of the non-hydrogen indices (if false we use a hydrogen index)
			useIdxEach: (int) The index to use in the list of hy/nonHyIndices
 
		"""
		self.nonHyAndHyClassifiers = nonHyAndHyClassifiers
		self.binValGetter = binValGetter
		self.useGroups = useGroups
		self.useNonHyIdx = useNonHyIdx
		self.useIdxEach = useIdxEach

	def getValsToBin(self, sparseMatrixCalculator):
		#1) Get groups
		groupIndices = [x.classify(sparseMatrixCalculator) for x in self.nonHyAndHyClassifiers]

		#2) Modify the binval getter with the relevant group indices
		currKwargs = {"useNonHyIdx":self.useNonHyIdx, "useIdxEach":self.useIdxEach}
		_modBinvalGetterWithFilterIndicesNonHyAndHyGenericDict(self.binValGetter, groupIndices, self.useGroups, useNonHyIdx=self.useNonHyIdx, useIdxEach=self.useIdxEach)

		#3) Get values using the modified getter
		outVals = self.binValGetter.getValsToBin(sparseMatrixCalculator)

		return outVals



#Modification for atom-atom case
@BINVAL_GETTER_TYPE_TO_ATOM_ATOM_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._RadialDistsGetValsToBin)
def _(binValGetter, groupIndices, useGroups):
	fromIndices = groupIndices[useGroups[0]]
	toIndices = groupIndices[useGroups[1]]

	binValGetter.indicesA, binValGetter.indicesB = fromIndices, toIndices


@BINVAL_GETTER_TYPE_TO_ATOM_ATOM_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._MinHozDistsGetValsToBin)
@BINVAL_GETTER_TYPE_TO_ATOM_ATOM_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._HozDistsGetValsToBin)
def _(binValGetter, groupIndices, useGroups):
	fromIndices = groupIndices[useGroups[0]]
	toIndices = groupIndices[useGroups[1]]

	binValGetter.fromIndices, binValGetter.toIndices = fromIndices, toIndices

@BINVAL_GETTER_TYPE_TO_ATOM_ATOM_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._PlanarDistsGetOneDimValsToBin)
def _(binValGetter, groupIndices, useGroups):
	indices = groupIndices[useGroups[0]]
	binValGetter.planeDistIndices = indices


#Modification function for the generic nonhy-hy case
@BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._PlanarDistsGetOneDimValsToBin)
def _(binvalGetter, groupIndices, useGroups, useNonHyIdx, useIdxEach):
	indices = groupIndices[useGroups[0]]
	useIndices = _getAtomicIndicesForNonHyToHyGeneric(indices, useNonHyIdx, useIdxEach)
	binvalGetter.planeDistIndices = useIndices

@BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._WaterOrientationBinValGetter)
def _(binvalGetter, groupIndices, useGroups, **kwargs):
	allIndices = groupIndices[useGroups[0]]
	assert all([len(x)==1 for x in allIndices[0]])
	oxyIndices = [x[0] for x in allIndices[0]]
	binvalGetter.oxyIndices = oxyIndices


@BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._CountNWithinDistancesGetValsToBin)
def _(binvalGetter, groupIndices, useGroups, useNonHyIdx, useIdxEach):
	#Sort out mapping to attributes
	typeToFromAttr = {atomComboBinvalGetterHelp._RadialDistsGetValsToBin:"indicesA",
	                  atomComboBinvalGetterHelp._HozDistsGetValsToBin:"fromIndices"}

	typeToToAttr = {atomComboBinvalGetterHelp._RadialDistsGetValsToBin:"indicesB",
	                  atomComboBinvalGetterHelp._HozDistsGetValsToBin:"toIndices"}

	fromIndicesAttr = typeToFromAttr[type(binvalGetter.getDistsBinValGetter)]
	toIndicesAttr = typeToToAttr[type(binvalGetter.getDistsBinValGetter)]

	#Figure out the indices to use and modify the binval getter
	allFromIndices = groupIndices[useGroups[0]]
	useFromIndices = _getAtomicIndicesForNonHyToHyGeneric(allFromIndices, useNonHyIdx, useIdxEach)
	setattr(binvalGetter.getDistsBinValGetter, fromIndicesAttr, useFromIndices)

	if len(useGroups)==1:
		pass
	else:
		allToIndices = groupIndices[useGroups[1]]
		useToIndices = _getAtomicIndicesForNonHyToHyGeneric(allToIndices, useNonHyIdx, useIdxEach)
		setattr(binvalGetter.getDistsBinValGetter, toIndicesAttr, useToIndices)



#Same as hozdist except attrs are indicesA/indicesB instead of fromIndices/toIndices
@BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._RadialDistsGetValsToBin)
def _(binvalGetter, groupIndices, useGroups, useNonHyIdx, useIdxEach):
	allFromIndices = groupIndices[useGroups[0]] #indicesA
	useFromIndices = _getAtomicIndicesForNonHyToHyGeneric(allFromIndices, useNonHyIdx, useIdxEach)
	binvalGetter.indicesA = useFromIndices

	if len(useGroups)==1:
		pass
	else:
		allToIndices = groupIndices[useGroups[1]]
		useToIndices = _getAtomicIndicesForNonHyToHyGeneric(allToIndices, useNonHyIdx, useIdxEach)
		binvalGetter.indicesB = useToIndices

@BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._HozDistsGetValsToBin)
@BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._MinHozDistsGetValsToBin)
def _(binvalGetter, groupIndices, useGroups, useNonHyIdx, useIdxEach):
	allFromIndices = groupIndices[useGroups[0]]
	useFromIndices = _getAtomicIndicesForNonHyToHyGeneric(allFromIndices, useNonHyIdx, useIdxEach)
	binvalGetter.fromIndices = useFromIndices

	if len(useGroups)==1:
		pass
	else:
		allToIndices = groupIndices[useGroups[1]]
		useToIndices = _getAtomicIndicesForNonHyToHyGeneric(allToIndices, useNonHyIdx, useIdxEach)
		binvalGetter.toIndices = useToIndices


@BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._GetOOHAnglesForHBondsBinValGetter)
@BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._GetNonHyDistsForHBondsBinValGetter)
@BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._CountHBondsBetweenGenericGroupsBinValGetter)
def _(binvalGetter, groupIndices, useGroups, **unused):
	fromNonHyIndices, fromHyIndices = groupIndices[useGroups[0]][0], groupIndices[useGroups[0]][1]
	binvalGetter.fromNonHyIndices, binvalGetter.fromHyIndices = fromNonHyIndices, fromHyIndices

	if len(useGroups)==1:
		pass
	else:
		toNonHyIndices, toHyIndices = groupIndices[useGroups[1]][0], groupIndices[useGroups[1]][1]
		binvalGetter.toNonHyIndices = toNonHyIndices
		binvalGetter.toHyIndices = toHyIndices

def _getAtomicIndicesForNonHyToHyGeneric(inpIndices, useNonHyIdx, useIdxEach):
	typeIdx = 0 if useNonHyIdx else 1
	outIndices = [ inpVals[useIdxEach] for inpVals in inpIndices[typeIdx] ]
	return outIndices

@BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._GetDiatomHozDistsBinvalGetter)
@BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._GetDiatomAngleWithVectorBinvalGetter)
@BINVAL_GETTER_TYPE_TO_NONHY_HY_GENERIC_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._GetDiatomDistsBinvalGetter)
def _(binvalGetter, groupIndices, useGroups, **unused):
	nonHyIndices, hyIndices = groupIndices[useGroups[0]]
	assert all([len(x)==1 for x in nonHyIndices])
	assert all([len(x)==1 for x in hyIndices])
	outDiatomIndices = [ [nonHyIdx[0],hyIdx[0]] for nonHyIdx,hyIdx in it.zip_longest(nonHyIndices,hyIndices) ]
	binvalGetter.diatomIndices = outDiatomIndices


#Modification functions for generic water derivative
@BINVAL_GETTER_TYPE_TO_WATER_DERIVATIVE_MOD_DICT(atomComboBinvalGetterHelp._PlanarDistsGetOneDimValsToBin)
def _(binvalGetter, groupIndices, useGroups, useNonHyIdx, useIdxEach):
	raise NotImplementedError("")



#Modification functions for water-water case
@BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._WaterMinDistBinValGetter)
def _(binValGetter, groupIndices, useGroups, toIdxType=None):
	if toIdxType is None:
		raise ValueError("")

	binValGetter.oxyIndices = groupIndices[useGroups[0]][0]
	binValGetter.hyIndices  = groupIndices[useGroups[0]][1]

	if toIdxType.upper() == "O":
		binValGetter.toIndices = groupIndices[useGroups[1]][0]
	elif toIdxType.upper() == "H":
		binValGetter.toIndices = [x for x in it.chain(*groupIndices[useGroups[1]][1])]
	elif toIdxType.upper() == "ALL":
		binValGetter.toIndices = groupIndices[useGroups[1]][0] + [x for x in it.chain(*groupIndices[useGroups[1]][1])]
	else:
		raise ValueError("{} is an invalid value for self.toIdxType".format(self.toIdxType))

@BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._RadialDistsGetValsToBin)
def _(binValGetter, groupIndices, useGroups, toIdxType=None):
	if toIdxType is None:
		raise ValueError("")

	fromOxyIndices = groupIndices[useGroups[0]][0]
	toOxyIndices, toHyIndices = groupIndices[useGroups[1]][0], [idx for idx in it.chain(*groupIndices[useGroups[1]][1])]

	binValGetter.indicesA = fromOxyIndices
	binValGetter.indicesB = _getToIndicesForWaterWater(toOxyIndices, toHyIndices, toIdxType)


#Only differes from the _RadialDistsGetValsToBin since the attrs are named differently
@BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._MinHozDistsGetValsToBin)
@BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._HozDistsGetValsToBin)
def _(binValGetter, groupIndices, useGroups, toIdxType=None):
	if toIdxType is None:
		raise ValueError("")

	fromOxyIndices = groupIndices[useGroups[0]][0]
	toOxyIndices, toHyIndices = groupIndices[useGroups[1]][0], [idx for idx in it.chain(*groupIndices[useGroups[1]][1])]

	binValGetter.fromIndices = fromOxyIndices
	binValGetter.toIndices = _getToIndicesForWaterWater(toOxyIndices, toHyIndices, toIdxType)


@BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._CountHBondsBetweenGenericGroupsBinValGetter)
def _(binValGetter, groupIndices, useGroups, toIdxType=None):
	fromIndices, toIndices = groupIndices[useGroups[0]], groupIndices[useGroups[1]]
	binValGetter.fromNonHyIndices = [[x] for x in fromIndices[0]]
	binValGetter.fromHyIndices = fromIndices[1]
	binValGetter.toNonHyIndices = [[x] for x in toIndices[0]]
	binValGetter.toHyIndices = toIndices[1]


@BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._WaterOrientationBinValGetter)
def _(binValGetter, groupIndices, useGroups, toIdxType=None):
	oxyIndices = groupIndices[ useGroups[0] ][0] #Should be only one so....
	binValGetter.oxyIndices = oxyIndices


@BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_REGISTER_DECO(atomComboBinvalGetterHelp._PlanarDistsGetOneDimValsToBin)
def _(binValGetter, groupIndices, useGroups, toIdxType=None):
	if toIdxType is None:
		raise ValueError("")

	toOxyIndices, toHyIndices = groupIndices[useGroups[0]][0],  [idx for idx in it.chain(*groupIndices[useGroups[0]][1])]

	toIndices = _getToIndicesForWaterWater(toOxyIndices, toHyIndices, toIdxType)
	binValGetter.planeDistIndices = toIndices


def _getToIndicesForWaterWater(toOxyIndices, toHyIndices, toIdxType):
	if toIdxType.upper() == "O":
		outIndices = toOxyIndices
	elif toIdxType.upper() == "H":
		outIndices = toHyIndices
	elif toIdxType.upper() == "ALL":
		outIndices = toOxyIndices + toHyIndices
	else:
		raise ValueError("{} is an invalid value for self.toIdxType".format(self.toIdxType))

	return outIndices


