
import itertools as it

from . import atom_combo_core as coreComboHelp
from . import atom_combo_opts_obj_maps as atomComboOptObjMaps
from . import atom_combo_populators as atomComboPopulatorHelp
from . import filtered_atom_combo_opt_objs as filteredAtomComboOptHelp
from . import classification_distr_opt_objs as classDistrOptObjHelp
from . import classification_binval_getters as classBinvalGetterHelp
from . import filtered_atom_combo_binval_getters as filteredAtomBinvalGetterHelp

from ..shared import register_key_decorator as regKeyDecoHelp

#May only even be neccesary for water-water case specifically???
_MOD_POPULATOR_BASED_ON_TYPE_DICT = dict()
MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO = regKeyDecoHelp.RegisterKeyValDecorator(_MOD_POPULATOR_BASED_ON_TYPE_DICT)


@atomComboOptObjMaps.TYPE_TO_POPULATOR_REGISTER_DECO(filteredAtomComboOptHelp.FilteredAtomComboOptsObjGeneric)
def _(inpObj):
	distrOptsPopulators = [atomComboOptObjMaps.getMatrixPopulatorFromOptsObj(optObj) for optObj in inpObj.distrOpts]
	classifierPopulators = [atomComboOptObjMaps.getMatrixPopulatorFromOptsObj(inpObj.classificationOpts)] 
	assert inpObj.classificationOpts.atomIndices == inpObj.atomIndices

	#Problem: These populators totally ignore the atomIndices on the overall index object.
	#To be on the safe side I modify them here to get everything needed
	for populator in distrOptsPopulators:
		_MOD_POPULATOR_BASED_ON_TYPE_DICT[type(populator), type(inpObj)] (populator, inpObj)

	#Combine individual populators into a composite
	outPopulator = coreComboHelp._SparseMatrixPopulatorComposite( distrOptsPopulators + classifierPopulators )

	return outPopulator


@atomComboOptObjMaps.TYPE_TO_POPULATOR_REGISTER_DECO(filteredAtomComboOptHelp.WaterToWaterFilteredAtomComboOptsObjGeneric)
def _(inpObj):
	distrOptsPopulators = [atomComboOptObjMaps.getMatrixPopulatorFromOptsObj(optObj) for optObj in inpObj.distrOpts]
	classifierPopulators = [atomComboOptObjMaps.getMatrixPopulatorFromOptsObj(inpObj.classificationOpts)] 

	# Modify certain args on the default populators based on the options object;
	# the initial example involves water-water minDist, where the indices i calculate distances TO can vary based on the "toIdx" argument in the options object
	#I could later build in a try/except so that we can ignore cases where i havent registered a key; but for now 
	#this more cautious approach seems best (I may put a specific overide on the opts obj to skip this step for key errors later)
	for populator, toIdxType in it.zip_longest(distrOptsPopulators, inpObj.toIdxTypes):
		_MOD_POPULATOR_BASED_ON_TYPE_DICT[type(populator), type(inpObj)]( populator,inpObj, toIdxType )


	#Combine individual populators into a composite
	outPopulator = coreComboHelp._SparseMatrixPopulatorComposite( distrOptsPopulators + classifierPopulators )

	return outPopulator


#Note 1: The binval getters are only really partially created at this stage; they are modified each step to take into account which
#group each atom belongs to
#Note 2: In future we may need a proper interface to map to the classifier; but for now this implementation
#should probably cover everything
@atomComboOptObjMaps.TYPE_TO_BINNER_REGISTER_DECO(filteredAtomComboOptHelp.WaterToWaterFilteredAtomComboOptsObjGeneric)
def _(inpObj):
	if type(inpObj.classificationOpts) is not classDistrOptObjHelp.WaterCountTypesMinDistAndHBondSimpleOpts:
		raise NotImplementedError("")

	_checkGroupIndicesConsistent(inpObj)

	#1) Get water classifier objects
	classifiers = list()
	currOptObj = inpObj.classificationOpts
	for idx,unused in enumerate(currOptObj.distFilterRanges):
		currArgs = [currOptObj.oxyIndices, currOptObj.hyIndices, currOptObj.distFilterIndices,
		            currOptObj.distFilterRanges[idx], currOptObj.nDonorFilterRanges[idx],
		            currOptObj.nAcceptorFilterRanges[idx], currOptObj.nTotalFilterRanges[idx],
		            currOptObj.maxOOHBond, currOptObj.maxAngleHBond]
		currObj = classBinvalGetterHelp._WaterClassifierMinDistAndNumberHBonds(*currArgs)
		classifiers.append(currObj)

	#2) Get the binval getters with default arguments (i.e. unfiltered indices or maybe even dud indices)
	binValGetters = [atomComboOptObjMaps.getOneDimBinValGetterFromOptsObj(optObj) for optObj in inpObj.distrOpts]

	#Create the actual objects (need 1 binval getter per property)
	outObjs = list()
	for binValGetter, useGroup, toIdxType in it.zip_longest(binValGetters, inpObj.useGroups, inpObj.toIdxTypes):
		currArgs = [ classifiers, binValGetter, useGroup, toIdxType ]
		outObjs.append( filteredAtomBinvalGetterHelp.WaterToWaterFilteredAtomComboBinvalGetterGeneric(*currArgs) )
	return outObjs


@atomComboOptObjMaps.TYPE_TO_BINNER_REGISTER_DECO(filteredAtomComboOptHelp.FilteredAtomComboOptsObjGeneric)
def _(inpObj):
	if type(inpObj.classificationOpts) is not classDistrOptObjHelp.AtomClassifyBasedOnDistsFromIndicesSimpleOpts:
		raise NotImplementedError("")

	_checkGroupIndicesConsistent(inpObj)

	#1) Get individual classifier objects[TODO: This is something another registration-dictionary COULD handle in future]
	classifiers = list()
	currOptObj = inpObj.classificationOpts
	for idx,unused in enumerate(currOptObj.distFilterRanges):
		currArgs = [currOptObj.atomIndices, currOptObj.distFilterIndices, currOptObj.distFilterRanges[idx]]
		currObj = classBinvalGetterHelp._AtomsWithinMinDistRangeClassifier(*currArgs, minDistVal=currOptObj.minDistVal)
		classifiers.append( currObj )

	#2) Get the binval getters with default arguments
	binValGetters = [atomComboOptObjMaps.getOneDimBinValGetterFromOptsObj(optObj) for optObj in inpObj.distrOpts]

	#3) Create the actual objects (need 1 binval getter per property)
	outObjs = list()
	for binValGetter, useGroup in it.zip_longest(binValGetters, inpObj.useGroups):
		currArgs = [ classifiers, binValGetter, useGroup ]
		outObjs.append( filteredAtomBinvalGetterHelp.FilteredAtomComboBinvalGetterGeneric(*currArgs) )

	return outObjs


def _checkGroupIndicesConsistent(inpObj):
	firstGroupIndices = [x[0] for x in inpObj.useGroups]
	if not all([x==firstGroupIndices[0] for x in firstGroupIndices]):
		raise ValueError("firstGroupIndices (in inpObj.useGroups) need to all be the same; but are {}".format(firstGroupIndices))


#Functions for modifying populators; this needs doing at creation time (since we populate matrices BEFORE we filter atoms into groups for a given geometry)
@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._WaterMinDistPopulator,filteredAtomComboOptHelp.WaterToWaterFilteredAtomComboOptsObjGeneric) )
def _(populator, optsObj, toIdxType):
	oxyIndices, hyIndices = optsObj.oxyIndices, optsObj.hyIndices
	if toIdxType.upper() == "O":
		toIndices = oxyIndices
	elif toIdxType.upper() == "H":
		toIndices = [idx for idx in it.chain(*hyIndices)]
	elif toIdxType.upper() == "ALL":
		toIndices = oxyIndices + [idx for idx in it.chain(*hyIndices)]
	else:
		raise ValueError("{} is an invalid value to toIdxType".format(toIdxType))

	populator.toIndices = toIndices


@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._HozDistMatrixPopulator, filteredAtomComboOptHelp.WaterToWaterFilteredAtomComboOptsObjGeneric) )
@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._DistMatrixPopulator   , filteredAtomComboOptHelp.WaterToWaterFilteredAtomComboOptsObjGeneric) )
def _(populator, optsObj, toIdxType):
	oxyIndices, hyIndices = optsObj.oxyIndices, optsObj.hyIndices
	fromIndices = oxyIndices
	if toIdxType.upper() == "O":
		toIndices = oxyIndices
	elif toIdxType.upper() == "H":
		toIndices = [idx for idx in it.chain(*hyIndices)]
	elif toIdxType.upper() == "ALL":
		toIndices = oxyIndices + [idx for idx in it.chain(*hyIndices)]
	else:
		raise ValueError("{} is an invalid value to toIdxType".format(toIdxType))

	populator.fromIndices = fromIndices
	populator.toIndices = toIndices


@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._HozDistMatrixPopulator, filteredAtomComboOptHelp.FilteredAtomComboOptsObjGeneric) )
@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._DistMatrixPopulator, filteredAtomComboOptHelp.FilteredAtomComboOptsObjGeneric) )
def _(populator, optsObj):
	populator.fromIndices = optsObj.atomIndices
	populator.toIndices = optsObj.atomIndices

