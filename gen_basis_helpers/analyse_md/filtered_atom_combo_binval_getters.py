
import itertools as it

from . import atom_combo_core as atomComboCoreHelp
from . import atom_combo_opts_obj_maps as atomComboOptObjMapHelp
from . import atom_combo_binval_getters as atomComboBinvalGetterHelp

from ..shared import register_key_decorator as regKeyDecoHelp

_BINVAL_GETTER_TYPE_TO_ATOM_ATOM_MOD_DICT = dict()
_BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_DICT = dict()

BINVAL_GETTER_TYPE_TO_ATOM_ATOM_MOD_REGISTER_DECO = regKeyDecoHelp.RegisterKeyValDecorator(_BINVAL_GETTER_TYPE_TO_ATOM_ATOM_MOD_DICT)
BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_REGISTER_DECO = regKeyDecoHelp.RegisterKeyValDecorator(_BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_DICT)


def _modBinvalGetterWithFilterIndicesAtomAtomDict(binValGetter, groupIndices, useGroups,**kwargs):
	_BINVAL_GETTER_TYPE_TO_ATOM_ATOM_MOD_DICT[type(binValGetter)](binValGetter, groupIndices, useGroups, **kwargs)

#Note: This is specifically for getting distributions between water of different types
#Different code would be needed to get (for example) water of one type to a fixed set of indices
def _modBinvalGetterWithFilterIndicesWaterWaterDict(binValGetter, groupIndices, useGroups, **kwargs):
	_BINVAL_GETTER_TYPE_TO_WATER_WATER_MOD_DICT[type(binValGetter)](binValGetter, groupIndices, useGroups, **kwargs)
	



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


