

#Note: This import is to register any functions for creating populators/binval getters
# related to options in this module; rather than to directly use any functions from it
#from . import filtered_atom_combo_obj_maps

#May or may not keep this; depends really
class _FilteredAtomComboOptsObjBase():
	pass



class FilteredAtomComboOptsObjGeneric():
	""" Generic options object used to get distributions between atoms based on how they are classified by the provided "classifier" """

	def __init__(self, atomIndices, classificationOpts, distrOpts, useGroups, classificationObjs=None):
		""" Initializer
		
		Args:
			atomIndices: (iter of ints) The indices of all atoms you want to include
			classificationOpts: (Options object for classification) This contains options for how we classify the indices into N-groups 
			distrOpts: (iter of CalcDistribOptionsBase) Each defines a distribution we want calculated using filtered lists of atomIndices
			useGroups: (iter of int-iters) Groups to calculate between. The group indices are determined by "classificationOpts". E.g. [ [0,1], [0] ] Would indicate to calculate distribution between groups [0,1] for distrOpts[0] and for group 0 for distrOpts[1]
			classificationObjs: (iter of ClassifierBase objects) If present these take priority over classificationOpts. Original purpose was to allow "byReference" classifiers to be used to speed up code when multiple of the same classifiers were used

		Notes:
			a) Even if classificationObjs are passed,  classificationOpts is still needed for figuring out which parts of various matrices to populate (e.g. which part of the distance matrix)

		"""
		self.atomIndices = atomIndices
		self.classificationOpts = classificationOpts
		self.distrOpts = distrOpts
		self.useGroups = useGroups
		self.classificationObjs = classificationObjs

	@property
	def primaryIndices(self):
		return self.atomIndices

	@property
	def binResObj(self):
		return [x.binResObj for x in self.distrOpts]

#May need an option to switch "primaryIdxType". For now its presumably always oxygen but...
#^^Easy to modify with an option later; and easy to set a default since should be the same for all
class WaterToWaterFilteredAtomComboOptsObjGeneric():
	""" Class used to specify options for how to calculate distributions between dynamically assigned (e.g. based on n-hbonds at one step) groups of water molecules """

	def __init__(self, oxyIndices, hyIndices, toIdxTypes, classificationOpts, distrOpts, useGroups, classificationObjs=None):
		""" Initializer 
		
		Args:
			oxyIndices: (iter of ints)
			hyIndices: (iter of int-iters)
			toIdxType: (iter of str) For each distr opt obj can be "O","H" or "all". When calcualting (for example) rdf between water groups this determines whether the second group contains oxygen indices, hydrogen indices or all water indices 
			classificationOpts: (Options object for classification) This contains options for how we classify the indices into N-groups 
			distrOpts: (iter of CalcDistribOptionsBase) Each defines a distribution we want calculated using filtered list of oxyIndices
			useGroups: (iter of int-iters) Groups to calculate between. The group indices are determined by "classificationOpts". E.g. [ [0,1], [0] ] Would indicate to calculate distribution between groups [0,1] for distrOpts[0] and for group 0 for distrOpts[1]
			classificationObjs: (iter of ClassifierBase objects) If present these take priority over classificationOpts. Original purpose was to allow "byReference" classifiers to be used to speed up code when multiple of the same classifiers were used

		NOTES:
			a) useGroups: The first index needs to be the same for all of them
			b) Even if classificationObjs are passed,  classificationOpts is still needed for figuring out which parts of various matrices to populate (e.g. which part of the distance matrix)

		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.toIdxTypes = toIdxTypes
		self.classificationOpts = classificationOpts
		self.distrOpts = distrOpts
		self.useGroups = useGroups
		self.classificationObjs = classificationObjs

	@property
	def primaryIndices(self):
		return self.oxyIndices

	@property
	def binResObj(self):
		return [x.binResObj for x in self.distrOpts]
