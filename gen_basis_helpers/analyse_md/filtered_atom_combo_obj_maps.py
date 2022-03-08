
import copy
import itertools as it

from . import atom_combo_core as coreComboHelp
from . import atom_combo_opts_obj_maps as atomComboOptObjMaps
from . import atom_combo_populators as atomComboPopulatorHelp
from . import filtered_atom_combo_opt_objs as filteredAtomComboOptHelp
from . import classification_distr_opt_objs as classDistrOptObjHelp
from . import classification_binval_getters as classBinvalGetterHelp
from . import classifier_objs as classifierObjHelp
from . import filtered_atom_combo_binval_getters as filteredAtomBinvalGetterHelp
from . import water_derivative_species_classifiers as waterDerivClassifierHelp

from ..shared import register_key_decorator as regKeyDecoHelp

#May only even be neccesary for water-water case specifically???
_MOD_POPULATOR_BASED_ON_TYPE_DICT = dict()
MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO = regKeyDecoHelp.RegisterKeyValDecorator(_MOD_POPULATOR_BASED_ON_TYPE_DICT)

#
_CLASSIFIER_FROM_OPTS_OBJ_FROM_TYPE_DICT = dict()
TYPE_TO_CLASSIFIER_REGISTER_DECO = regKeyDecoHelp.RegisterKeyValDecorator(_CLASSIFIER_FROM_OPTS_OBJ_FROM_TYPE_DICT)

def getClassifiersFromOptsObj(inpOptsObj):
	return _CLASSIFIER_FROM_OPTS_OBJ_FROM_TYPE_DICT[type(inpOptsObj)](inpOptsObj)


@TYPE_TO_CLASSIFIER_REGISTER_DECO(classDistrOptObjHelp.CompositeClassiferOptsSimple)
def _(inpObj):
	classifiers = list()
	for optsObj in inpObj.optsObjs:
		currClassifiers = getClassifiersFromOptsObj(optsObj)
		classifiers.extend(currClassifiers)
	return classifiers

@TYPE_TO_CLASSIFIER_REGISTER_DECO(classDistrOptObjHelp.WaterCountTypesMinDistAndHBondSimpleOpts)
def _(inpObj):
	classifiers = list()
	for idx,unused in enumerate(inpObj.distFilterRanges):
		currArgs = [inpObj.oxyIndices, inpObj.hyIndices, inpObj.distFilterIndices,
		            inpObj.distFilterRanges[idx], inpObj.nDonorFilterRanges[idx],
		            inpObj.nAcceptorFilterRanges[idx], inpObj.nTotalFilterRanges[idx],
		            inpObj.maxOOHBond, inpObj.maxAngleHBond]
		currObj = classifierObjHelp._WaterClassifierMinDistAndNumberHBonds(*currArgs)
		classifiers.append(currObj)
	return classifiers


@TYPE_TO_CLASSIFIER_REGISTER_DECO(classDistrOptObjHelp.WaterAdsorbedClassifier_usingMinHozDistsBetweenAdsorptionSitesOptsObj)
def _(inpObj):
	classifiers = list()
	for idx,unused in enumerate(inpObj.distFilterRanges):
		currArgs = [inpObj.oxyIndices, inpObj.hyIndices, inpObj.distFilterIndices, inpObj.distFilterRanges[idx],
		            inpObj.nDonorFilterRanges[idx], inpObj.nAcceptorFilterRanges[idx], inpObj.nTotalFilterRanges[idx],
		            inpObj.maxOOHBond, inpObj.maxAngleHBond, inpObj.adsSiteMinHozToOtherAdsSiteRanges[idx]]
		currObj = classifierObjHelp._WaterClassifierMinDistHBondsAndAdsSiteHozDists(*currArgs)
		classifiers.append(currObj)
	return classifiers


@TYPE_TO_CLASSIFIER_REGISTER_DECO(classDistrOptObjHelp.AtomClassifyBasedOnDistsFromIndicesSimpleOpts)
def _(inpObj):
	classifiers = list()
	for idx,unused in enumerate(inpObj.distFilterRanges):
		currArgs = [inpObj.atomIndices, inpObj.distFilterIndices, inpObj.distFilterRanges[idx]]
		currObj = classifierObjHelp._AtomsWithinMinDistRangeClassifier(*currArgs, minDistVal=inpObj.minDistVal)
		classifiers.append( currObj )
	return classifiers


@TYPE_TO_CLASSIFIER_REGISTER_DECO(classDistrOptObjHelp.ClassifyBasedOnHBondingToGroup_simple)
def _(inpObj):
	classifiers = list()
	for idx, unused in enumerate(inpObj.nDonorFilterRanges):
		currArgs = [inpObj.fromNonHyIndices, inpObj.fromHyIndices, inpObj.toNonHyIndices, inpObj.toHyIndices,
		            inpObj.nDonorFilterRanges[idx], inpObj.nAcceptorFilterRanges[idx],
		            inpObj.nTotalFilterRanges[idx], inpObj.maxOOHBond, inpObj.maxAngleHBond]
		currObj = classifierObjHelp._GenericNonHyAndHyClassiferUsingHBondsToGroup_simple(*currArgs)
		classifiers.append( currObj )
	return classifiers

#Note: Not too sure about the inpObj.firstClassifierObjs argument; untested at time of writing
@TYPE_TO_CLASSIFIER_REGISTER_DECO(classDistrOptObjHelp.ClassifyBasedOnHBondingToDynamicGroup)
def _(inpObj):

	firstGroupClassifiers  = getClassifiersFromOptsObj(inpObj.dynGroupOptObj)
	secondGroupClassifiers = getClassifiersFromOptsObj(inpObj.thisGroupOptObj)
	outClassifiers = list()


	for firstClass, secondClass in it.zip_longest(firstGroupClassifiers, secondGroupClassifiers):
		currArgs = [firstClass, secondClass]
		currKwargs = {"mutuallyExclusive":inpObj.mutuallyExclusive}
		currClassifier = classifierObjHelp._ClassiferUsingHBondsToDynamicGroup(*currArgs, **currKwargs)
		outClassifiers.append( currClassifier )

	return outClassifiers


@TYPE_TO_CLASSIFIER_REGISTER_DECO(classDistrOptObjHelp.ClassifyNonHyAndHyBasedOnMinHozDistsToAtomGroups)
def _(inpObj):
	outClassifiers = list()

	sharedKwargs = {"useIndicesFrom":inpObj.useIndicesFrom, "useIndicesTo":inpObj.useIndicesTo, "minDistVal":inpObj.minDistVal}
	for idx,unused in enumerate(inpObj.distFilterRanges):
		currArgs = [inpObj.fromNonHyIndices, inpObj.fromHyIndices, inpObj.toNonHyIndices, inpObj.toHyIndices,
		            inpObj.distFilterRanges[idx]]
		currClassifier = classifierObjHelp._NonHyAndHyClassifierBasedOnMinHozDistToGroup(*currArgs, **sharedKwargs)
		outClassifiers.append( currClassifier )

	return outClassifiers

@TYPE_TO_CLASSIFIER_REGISTER_DECO(classDistrOptObjHelp.ClassifyNonHyAndHyChainedAllCommon)
def _(inpObj):
	useClassifiers = [x for x in it.chain(*[getClassifiersFromOptsObj(x) for x in inpObj.classifierOpts])]
	outClassifier = classifierObjHelp._NonHyAndHyChainedClassifier_allCommon(useClassifiers)
	return [outClassifier]


@TYPE_TO_CLASSIFIER_REGISTER_DECO(classDistrOptObjHelp.WaterDerivativeBasedOnDistanceClassifierOptsObj)
def _(inpObj):
	classifiers = list()
	for idx, unused in enumerate(inpObj.nNebs):
		currArgs = [inpObj.oxyIndices, inpObj.hyIndices]
		currKwargs = {"maxOHDist":inpObj.maxOHDist, "nNebs": inpObj.nNebs[idx]}
		currObj = waterDerivClassifierHelp._WaterDerivativeDistanceOnlyClassifierGeneric(*currArgs, **currKwargs)
		classifiers.append( currObj )
	return classifiers



#
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

#NOTE: Same as for "FilteredAtomComboOptsObjGeneric" initially; put separately still though
@atomComboOptObjMaps.TYPE_TO_POPULATOR_REGISTER_DECO(filteredAtomComboOptHelp.HydroxylDiatomFromNonHyAndHyFilteredOptsObj_simple)
@atomComboOptObjMaps.TYPE_TO_POPULATOR_REGISTER_DECO(filteredAtomComboOptHelp.WaterDerivativeFilteredOptsObj_simple)
@atomComboOptObjMaps.TYPE_TO_POPULATOR_REGISTER_DECO(filteredAtomComboOptHelp.GenericNonHyAndHyFilteredOptsObj_simple)
def _(inpObj):
	distrOptsPopulators = [atomComboOptObjMaps.getMatrixPopulatorFromOptsObj(optObj) for optObj in inpObj.distrOpts]
	classifierPopulators = [atomComboOptObjMaps.getMatrixPopulatorFromOptsObj(inpObj.classificationOpts)]

	#Want to ensure the indices on populators match that of the primaryIndices for what we're filtering
	#Probably a similar problem for the classifier populators too actually but....
	for populator in distrOptsPopulators: 
		_MOD_POPULATOR_BASED_ON_TYPE_DICT[type(populator), type(inpObj)](populator,inpObj)

	#Combine individual into a composite
	outPopulator = coreComboHelp._SparseMatrixPopulatorComposite( distrOptsPopulators + classifierPopulators )

	return outPopulator

#SLIGHTLY different to *_simple since we're only allowed 1 distrOpts here
@atomComboOptObjMaps.TYPE_TO_POPULATOR_REGISTER_DECO(filteredAtomComboOptHelp.GenericNonHyAndHyFilteredOptsObj_getAverageVal)
def _(inpObj):
	distrOptsPopulator = atomComboOptObjMaps.getMatrixPopulatorFromOptsObj(inpObj.distrOpts)
	classifierPopulators = [atomComboOptObjMaps.getMatrixPopulatorFromOptsObj(inpObj.classificationOpts)]

	#Want to ensure the indices on populators match that of the primaryIndices for what we're filtering
	#Probably a similar problem for the classifier populators too actually but....
	_MOD_POPULATOR_BASED_ON_TYPE_DICT[type(distrOptsPopulator), type(inpObj)](distrOptsPopulator, inpObj)

	#Combine individual into a composite
	outPopulator = coreComboHelp._SparseMatrixPopulatorComposite( [distrOptsPopulator] + classifierPopulators )

	return outPopulator

#Note 1: The binval getters are only really partially created at this stage; they are modified each step to take into account which
#group each atom belongs to
#Note 2: In future we may need a proper interface to map to the classifier; but for now this implementation
#should probably cover everything
@atomComboOptObjMaps.TYPE_TO_BINNER_REGISTER_DECO(filteredAtomComboOptHelp.WaterToWaterFilteredAtomComboOptsObjGeneric)
def _(inpObj):

	_checkGroupIndicesConsistent(inpObj)

	#1) Get water classifier objects
	if inpObj.classificationObjs is not None:
		classifiers = inpObj.classificationObjs
	else:
		classifiers = getClassifiersFromOptsObj(inpObj.classificationOpts)

#	classifiers = getClassifiersFromOptsObj(inpObj.classificationOpts)

	#2) Get the binval getters with default arguments (i.e. unfiltered indices or maybe even dud indices)
	binValGetters = [atomComboOptObjMaps.getOneDimBinValGetterFromOptsObj(optObj) for optObj in inpObj.distrOpts]

	#Create the actual objects (need 1 binval getter per property)
	#Note: We only use the actual classifiers for the first case and refernce for the others; since there all the same
	outObjs = list()
	useClassifiers = [classifiers] + [classifierObjHelp.getByReferenceClassifiers(classifiers) for x in range(1,len(binValGetters))]
	for binValGetter, useGroup, toIdxType, currClassifiers in it.zip_longest(binValGetters, inpObj.useGroups, inpObj.toIdxTypes, useClassifiers):
		currArgs = [ currClassifiers, binValGetter, useGroup, toIdxType ]
		outObjs.append( filteredAtomBinvalGetterHelp.WaterToWaterFilteredAtomComboBinvalGetterGeneric(*currArgs) )
	return outObjs



#Only the "GenericNonHyAndHyFilteredAtomComboBinvalGetter_simple" line differs from the atomic case
@atomComboOptObjMaps.TYPE_TO_BINNER_REGISTER_DECO(filteredAtomComboOptHelp.HydroxylDiatomFromNonHyAndHyFilteredOptsObj_simple)
@atomComboOptObjMaps.TYPE_TO_BINNER_REGISTER_DECO(filteredAtomComboOptHelp.WaterDerivativeFilteredOptsObj_simple)
@atomComboOptObjMaps.TYPE_TO_BINNER_REGISTER_DECO(filteredAtomComboOptHelp.GenericNonHyAndHyFilteredOptsObj_simple)
def _tempBinvalGetter_genericNonHyAndHyFiltered(inpObj):

	_checkGroupIndicesConsistent(inpObj)

	#1) Get individual classifier objects
	try:
		classificationObjs = inpObj.classificationObjs
	except AttributeError:
		classifiers = getClassifiersFromOptsObj(inpObj.classificationOpts)
	else:
		if inpObj.classificationObjs is not None:
			classifiers = inpObj.classificationObjs
		else:
			classifiers = getClassifiersFromOptsObj(inpObj.classificationOpts)

	#2) Get the binval getters with default arguments
	binValGetters = [atomComboOptObjMaps.getOneDimBinValGetterFromOptsObj(optObj) for optObj in inpObj.distrOpts]

	#3) Create the actual objects (need 1 binval getter per property)
	#Note: We only use the actual classifiers for the first case and refernce for the others; since there all the same
	outObjs = list()
	useClassifiers = [classifiers] + [classifierObjHelp.getByReferenceClassifiers(classifiers) for x in range(1,len(binValGetters))]

	for binValGetter, useGroup, currClassifiers in it.zip_longest(binValGetters, inpObj.useGroups, useClassifiers):
		currArgs = [ currClassifiers, binValGetter, useGroup ]
		currKwargs = {"useNonHyIdx":inpObj.useNonHyIdx, "useIdxEach":inpObj.useIdxEach}
		outObjs.append( filteredAtomBinvalGetterHelp.GenericNonHyAndHyFilteredAtomComboBinvalGetter_simple(*currArgs,**currKwargs) )

	return outObjs

@atomComboOptObjMaps.TYPE_TO_BINNER_REGISTER_DECO(filteredAtomComboOptHelp.GenericNonHyAndHyFilteredOptsObj_getAverageVal)
def _(inpObj):
	#Temporarily changing distrOpts so it works with old interface; using copy destroys classification obj use
	inpObj.distrOpts = [inpObj.distrOpts]
	basicObjs = _tempBinvalGetter_genericNonHyAndHyFiltered(inpObj)
	inpObj.distrOpts = inpObj.distrOpts[0]


	outObjs = [filteredAtomBinvalGetterHelp.GenericNonHyAndHyFilteredAtomComboBinvalGetter_simple_getAverage(x) for x in basicObjs]
	if len(outObjs)>1:
		raise NotImplementedError("Likely no harm just removing this line; but i couldnt be bothered to think about if thatd cause any issues")
	return outObjs

@atomComboOptObjMaps.TYPE_TO_BINNER_REGISTER_DECO(filteredAtomComboOptHelp.FilteredAtomComboOptsObjGeneric)
def _(inpObj):

	_checkGroupIndicesConsistent(inpObj)

	#1) Get individual classifier objects
	try:
		classificationObjs = inpObj.classificationObjs
	except AttributeError:
		classifiers = getClassifiersFromOptsObj(inpObj.classificationOpts)
	else:
		if inpObj.classificationObjs is not None:
			classifiers = inpObj.classificationObjs
		else:
			classifiers = getClassifiersFromOptsObj(inpObj.classificationOpts)

	#2) Get the binval getters with default arguments
	binValGetters = [atomComboOptObjMaps.getOneDimBinValGetterFromOptsObj(optObj) for optObj in inpObj.distrOpts]

	#3) Create the actual objects (need 1 binval getter per property)
	#Note: We only use the actual classifiers for the first case and refernce for the others; since there all the same
	outObjs = list()
	useClassifiers = [classifiers] + [classifierObjHelp.getByReferenceClassifiers(classifiers) for x in range(1,len(binValGetters))]

	for binValGetter, useGroup, currClassifiers in it.zip_longest(binValGetters, inpObj.useGroups, useClassifiers):
		currArgs = [ currClassifiers, binValGetter, useGroup ]
		outObjs.append( filteredAtomBinvalGetterHelp.FilteredAtomComboBinvalGetterGeneric(*currArgs) )

	return outObjs


def _checkGroupIndicesConsistent(inpObj):
	firstGroupIndices = [x[0] for x in inpObj.useGroups]
	if not all([x==firstGroupIndices[0] for x in firstGroupIndices]):
		raise ValueError("firstGroupIndices (in inpObj.useGroups) need to all be the same; but are {}".format(firstGroupIndices))


#Functions for modifying populators; this needs doing at creation time (since we populate matrices BEFORE we filter atoms into groups for a given geometry)

@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (coreComboHelp._SparseMatrixPopulatorComposite, filteredAtomComboOptHelp.GenericNonHyAndHyFilteredOptsObj_simple) )
def _(populator, optsObj):
	for subPopulator in populator.populators:
		_MOD_POPULATOR_BASED_ON_TYPE_DICT[type(subPopulator), type(optsObj)] (populator, optsObj)



@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._WaterMinDistPopulator,filteredAtomComboOptHelp.WaterToWaterFilteredAtomComboOptsObjGeneric) )
def _(populator, optsObj, toIdxType):
	populator.toIndices = _getToIndicesFromWaterToWaterOptsObjAndToIdxType(optsObj, toIdxType)


@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._HozDistMatrixPopulator, filteredAtomComboOptHelp.WaterToWaterFilteredAtomComboOptsObjGeneric) )
@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._DistMatrixPopulator   , filteredAtomComboOptHelp.WaterToWaterFilteredAtomComboOptsObjGeneric) )
def _(populator, optsObj, toIdxType):
	toIndices = _getToIndicesFromWaterToWaterOptsObjAndToIdxType(optsObj, toIdxType)
	oxyIndices, hyIndices = optsObj.oxyIndices, optsObj.hyIndices
	fromIndices = oxyIndices

	populator.fromIndices = fromIndices
	populator.toIndices = toIndices

@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._PlanarDistMatrixPopulator, filteredAtomComboOptHelp.WaterToWaterFilteredAtomComboOptsObjGeneric) )
def _(populator, optsObj, toIdxType):
	populator.indices = _getToIndicesFromWaterToWaterOptsObjAndToIdxType(optsObj, toIdxType)


def _getToIndicesFromWaterToWaterOptsObjAndToIdxType(optsObj, toIdxType):
	oxyIndices, hyIndices = optsObj.oxyIndices, optsObj.hyIndices
	if toIdxType.upper() == "O":
		toIndices = oxyIndices
	elif toIdxType.upper() == "H":
		toIndices = [idx for idx in it.chain(*hyIndices)]
	elif toIdxType.upper() == "ALL":
		toIndices = oxyIndices + [idx for idx in it.chain(*hyIndices)]
	else:
		raise ValueError("{} is an invalid value to toIdxType".format(toIdxType))
	
	return toIndices


@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._WaterOrientationPopulator, filteredAtomComboOptHelp.WaterToWaterFilteredAtomComboOptsObjGeneric) )
def _(populator, optsObj, toIdxType):
	populator.oxyIndices = optsObj.oxyIndices
	populator.hyIndices = optsObj.hyIndices


@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._CountHBondsBetweenGenericGroupsPopulator, filteredAtomComboOptHelp.WaterToWaterFilteredAtomComboOptsObjGeneric) )
def _(populator, optsObj, toIdxType):
	populator.fromNonHyIndices, populator.FromHyIndices = [[x] for x in optsObj.oxyIndices], optsObj.hyIndices
	populator.toNonHyIndices, populator.toHyIndices = [[x] for x in optsObj.oxyIndices], optsObj.hyIndices



@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._HozDistMatrixPopulator, filteredAtomComboOptHelp.FilteredAtomComboOptsObjGeneric) )
@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._DistMatrixPopulator, filteredAtomComboOptHelp.FilteredAtomComboOptsObjGeneric) )
def _(populator, optsObj):
	populator.fromIndices = optsObj.atomIndices
	populator.toIndices = optsObj.atomIndices

@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._PlanarDistMatrixPopulator, filteredAtomComboOptHelp.GenericNonHyAndHyFilteredOptsObj_getAverageVal) )
@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._PlanarDistMatrixPopulator, filteredAtomComboOptHelp.GenericNonHyAndHyFilteredOptsObj_simple) )
@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._PlanarDistMatrixPopulator, filteredAtomComboOptHelp.FilteredAtomComboOptsObjGeneric) )
def _(populator, optsObj):
	populator.indices = optsObj.atomIndices


# The generic nonhy/hy h-bond one
@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._WaterOrientationPopulator, filteredAtomComboOptHelp.GenericNonHyAndHyFilteredOptsObj_simple) )
def _(populator, optsObj):
	outOxyIndices, outHyIndices = list(), list()

	#Get oxygen indices
	for currIndices in optsObj.fromNonHyIndices:
		assert len(currIndices)==1, "{} non-hydrogen indices detected; This needs to be 1 to get water orientations".format( len(currIndices) )
		outOxyIndices.append( currIndices[0] )

	#Get hydrogen indices
	for currIndices in  optsObj.fromHyIndices:
		assert len(currIndices)==2, "{} hydrogen indices detected; This needs to be 2 to get water orientations".format( len(currIndices) )
		outHyIndices.append( currIndices )

	populator.oxyIndices = outOxyIndices
	populator.hyIndices = outHyIndices


@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._DistMatrixPopulator,    filteredAtomComboOptHelp.GenericNonHyAndHyFilteredOptsObj_simple) )
@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._HozDistMatrixPopulator, filteredAtomComboOptHelp.GenericNonHyAndHyFilteredOptsObj_simple) )
def _(populator, optsObj):
	populator.fromIndices = optsObj.atomIndices

	#We modify the toIndices if they wont get used OR if .fromIndices will be used for them (e.g. for calcualting distributions BETWEEN groups)
	outToIndices = list()
	bothGroupsUsed, singleGroupUsed = _areSingleAndOrTwoGroupsUsed(optsObj.useGroups)

	if bothGroupsUsed:
		outToIndices = optsObj.atomIndices
	if singleGroupUsed:
		outToIndices = list( set(outToIndices + populator.toIndices ) )

	populator.toIndices = outToIndices


@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._HozDistMatrixPopulator, filteredAtomComboOptHelp.HydroxylDiatomFromNonHyAndHyFilteredOptsObj_simple) )
@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._DistMatrixPopulator, filteredAtomComboOptHelp.HydroxylDiatomFromNonHyAndHyFilteredOptsObj_simple) )
def _(populator, optsObj):
	populator.fromIndices = [x[0] for x in optsObj.fromNonHyIndices]
	populator.toIndices = [ x[0] for x in optsObj.fromHyIndices ]

@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._DiatomAngleWithVectorPopulator, filteredAtomComboOptHelp.HydroxylDiatomFromNonHyAndHyFilteredOptsObj_simple) )
def _(populator, optsObj):
	assert all([len(x)==1 for x in optsObj.fromNonHyIndices])
	assert all([len(x)==1 for x in optsObj.fromHyIndices])
	diatomIndices = [ [nonHy[0],hy[0]] for nonHy,hy in it.zip_longest(optsObj.fromNonHyIndices, optsObj.fromHyIndices) ] 
	populator.diatomIndices = diatomIndices

@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._CountHBondsBetweenGenericGroupsPopulator, filteredAtomComboOptHelp.GenericNonHyAndHyFilteredOptsObj_getAverageVal) )
@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._CountHBondsBetweenGenericGroupsPopulator, filteredAtomComboOptHelp.GenericNonHyAndHyFilteredOptsObj_simple) )
def _(populator, optsObj):
	populator.fromNonHyIndices = optsObj.fromNonHyIndices
	populator.fromHyIndices = optsObj.fromHyIndices

	#
	outToNonHyIndices, outToHyIndices = list(), list()
	bothGroupsUsed, singleGroupUsed = _areSingleAndOrTwoGroupsUsed(optsObj.useGroups)
	if bothGroupsUsed:
		outToNonHyIndices, outToHyIndices = optsObj.fromNonHyIndices, optsObj.fromHyIndices
	if singleGroupUsed:
		outToNonHyIndices = list( set([tuple(x) for x in outToNonHyIndices] + [tuple(x) for x in populator.toNonHyIndices]) )
		outToHyIndices = list( set([tuple(x) for x in outToHyIndices] + [tuple(x) for x in populator.toHyIndices]) )

	#
	populator.toNonHyIndices, populator.toHyIndices = outToNonHyIndices, outToHyIndices


def _areSingleAndOrTwoGroupsUsed(useGroups):
	bothGroupsUsed, singleGroupUsed = False, False
	for currGroups in useGroups:
		if len(currGroups)==1:
			singleGroupUsed = True
		else:
			bothGroupsUsed = True

	return bothGroupsUsed, singleGroupUsed


#The water derivative classifier options cases
@MOD_POPULATOR_BASED_ON_TYPE_DICT_REGISTER_DECO( (atomComboPopulatorHelp._PlanarDistMatrixPopulator, filteredAtomComboOptHelp.WaterDerivativeFilteredOptsObj_simple) )
def _(populator, optsObj):
	if optsObj.useNonHyIdx:
		outIndices = optsObj.oxyIndices
	else:
		outIndices = optsObj.hyIndices
	populator.indices = outIndices



