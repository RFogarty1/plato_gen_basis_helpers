

#Note: This import is to register any functions for creating populators/binval getters
# related to options in this module; rather than to directly use any functions from it
#from . import filtered_atom_combo_obj_maps

#May or may not keep this; depends really
class _FilteredAtomComboOptsObjBase():
	pass



#May need an option to switch "primaryIdxType". For now its presumably always oxygen but...
#^^Easy to modify with an option later; and easy to set a default since should be the same for all
class WaterToWaterFilteredAtomComboOptsObjGeneric():
	""" Class used to specify options for how to calculate distributions between dynamically assigned (e.g. based on n-hbonds at one step) groups of water molecules """

	def __init__(self, oxyIndices, hyIndices, toIdxTypes, classificationOpts, distrOpts, useGroups):
		""" Initializer 
		
		Args:
			oxyIndices: (iter of ints)
			hyIndices: (iter of int-iters)
			toIdxType: (iter of str) For each distr opt obj can be "O","H" or "all". When calcualting (for example) rdf between water groups this determines whether the second group contains oxygen indices, hydrogen indices or all water indices 
			classificationOpts: (Options object for classification) This contains options for how we classify the indices into N-groups 
			distrOpts: (iter of CalcDistribOptionsBase) Each defines a distribution we want calculated using filtered list of oxyIndices
			useGroups: (iter of int-iters) Groups to calculate between. The group indices are determined by "classificationOpts". E.g. [ [0,1], [0] ] Would indicate to calculate distribution between groups [0,1] for distrOpts[0] and for group 0 for distrOpts[1]

		NOTES:
			useGroups: The first index needs to be the same for all of them
 
		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.toIdxTypes = toIdxTypes
		self.classificationOpts = classificationOpts
		self.distrOpts = distrOpts
		self.useGroups = useGroups


	@property
	def primaryIndices(self):
		return self.oxyIndices

	@property
	def binResObj(self):
		return [x.binResObj for x in self.distrOpts]

