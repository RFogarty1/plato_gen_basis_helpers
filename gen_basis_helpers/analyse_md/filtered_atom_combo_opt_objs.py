

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

#TODO: useGroups = [ [0] ] for something like rdf should probably calculate it to a set of fixed indices, whereas useGroups = [ [0,1] ] should calcualte between classified groups
#TODO: Probably want a separate class for dealing with dynamic "to" indices? Or maybe i 
#TODO: Should likely merge the "WaterToWaterFilteredAtomComboOptsObjGeneric" with this in future; at least make use the same backends
class GenericNonHyAndHyFilteredOptsObj_simple():
	""" Class used to specify options for how to calculate distributions between dynamically assigned groups of molecules, where classsification returns non-hydrogen/hydrogen atom indices separately (generally the case when counting h-bonds or similar)

	"""

	def __init__(self, fromNonHyIndices, fromHyIndices, classificationOpts, distrOpts, useGroups, useNonHyIdx=True, useIdxEach=0, classificationObjs=None):
		""" Initializer
		
		Args:
			fromNonHyIndices: (iter of int-iters) Each is the non-hydrogen indices from one group (e.g. may be [ [0] ] for a single water)
			fromHyIndices: (iter of int-iters) Each is hydrogen indices for one group (e.g. may be [ [1,2] ] for a single water)
			classificationOpts: (Options object for classification) This contains options for how we classify the indices into N-groups. Should lead to something that 
			distrOpts: (iter of CalcDistribOptionsBase) Each defines a distribution we want calculated using filtered list of oxyIndices
			useGroups: (iter of int-iters) Groups to calculate between. The group indices are determined by "classificationOpts". E.g. [ [0,1], [0] ] Would indicate to calculate distribution between groups [0,1] for distrOpts[0] and for group 0 for distrOpts[1]
			useNonHyIdx: (Bool) If True we represent our group with one of the non-hydrogen indices (if false we use a hydrogen index)
			useIdxEach: (Int) The index to use in the list of hy/nonHyIndices.
			classificationObjs: (iter of ClassifierBase objects) If present these take priority over classificationOpts. Original purpose was to allow "byReference" classifiers to be used to speed up code when multiple of the same classifiers were used

		"""
		self.fromNonHyIndices = fromNonHyIndices
		self.fromHyIndices = fromHyIndices
		self.classificationOpts = classificationOpts
		self.distrOpts = distrOpts
		self.useGroups = useGroups
		self.useNonHyIdx = useNonHyIdx
		self.useIdxEach = useIdxEach
		self.classificationObjs = classificationObjs

	@property
	def primaryIndices(self):
		if self.useNonHyIdx:
			outVals = [x[self.useIdxEach]  for x in self.fromNonHyIndices]
		else:
			outVals = [x[self.useIdxEach] for x in self.fromHyIndices]

		return outVals

	#Used to match the interface of FilteredAtomComboOptsObjGeneric as close as possible (lets me resue some other functions)
	@property
	def atomIndices(self):
		return self.primaryIndices

	@property
	def binResObj(self):
		return [x.binResObj for x in self.distrOpts]


class GenericNonHyAndHyFilteredOptsObj_getAverageVal():
	""" Class used to specify options for how to calculate an average value between dynamically assigned groups of molecules, where classsification returns non-hydrogen/hydrogen atom indices separately (generally the case when counting h-bonds or similar)

	Original use case was to count the average number of h-bonds formed for all water in a group (so you can plot time vs <h-bonds> for adsorbed water for example)

	"""
	def __init__(self, fromNonHyIndices, fromHyIndices, classificationOpts, distrOpts, useGroups, useNonHyIdx=True, useIdxEach=0, classificationObjs=None):
		""" Initializer
		
		Args:
			fromNonHyIndices: (iter of int-iters) Each is the non-hydrogen indices from one group (e.g. may be [ [0] ] for a single water)
			fromHyIndices: (iter of int-iters) Each is hydrogen indices for one group (e.g. may be [ [1,2] ] for a single water)
			classificationOpts: (Options object for classification) This contains options for how we classify the indices into N-groups. Should lead to something that 
			distrOpts: (CalcDistribOptionsBase) Each defines a distribution we want calculated using filtered list of indices
			useGroups: (int-iter) Groups to calculate between. The group indices are determined by "classificationOpts". E.g. [0,1] Would indicate to calculate distribution between groups zero and one
			useNonHyIdx: (Bool) If True we represent our group with one of the non-hydrogen indices (if false we use a hydrogen index)
			useIdxEach: (Int) The index to use in the list of hy/nonHyIndices.
			classificationObjs: (iter of ClassifierBase objects) If present these take priority over classificationOpts. Original purpose was to allow "byReference" classifiers to be used to speed up code when multiple of the same classifiers were used

		"""
		self.fromNonHyIndices = fromNonHyIndices
		self.fromHyIndices = fromHyIndices
		self.classificationOpts = classificationOpts
		self.distrOpts = distrOpts
		self.useGroups = useGroups
		self.useNonHyIdx = useNonHyIdx
		self.useIdxEach = useIdxEach
		self.classificationObjs = classificationObjs

	@property
	def primaryIndices(self):
		if self.useNonHyIdx:
			outVals = [x[self.useIdxEach]  for x in self.fromNonHyIndices]
		else:
			outVals = [x[self.useIdxEach] for x in self.fromHyIndices]

		return outVals

	#Used to match the interface of FilteredAtomComboOptsObjGeneric as close as possible (lets me resue some other functions)
	@property
	def atomIndices(self):
		return self.primaryIndices




class HydroxylDiatomFromNonHyAndHyFilteredOptsObj_simple():
	""" Class used to look at hydroxyl O-H properties based on classifications given by NonHy/hy classifiers. Original use case is to look at O-H/surface normal angle based on whether a hydrogen bond is formed or not

	This MIGHT be possible to do with the normal GenericNonHyAndHyFilteredOptsObj_simple aswell if you simply use len-1 indices for hy/nonhy. Likely that this options object is redundant.....

	"""

	def __init__(self, fromNonHyIndices, fromHyIndices, classificationOpts, distrOpts, useGroups, useNonHyIdx=True, classificationObjs=None):
		""" Initializer
		
		Args:
			fromNonHyIndices: (iter of int-iters) Each is the non-hydrogen indices from one group (e.g. may be [ [0] ] for a single water)
			fromHyIndices: (iter of int-iters) Each is hydrogen indices for one group (e.g. may be [ [1,2] ] for a single water)
			classificationOpts: (Options object for classification) This contains options for how we classify the indices into N-groups. Should lead to something that 
			distrOpts: (iter of CalcDistribOptionsBase) Each defines a distribution we want calculated using filtered list; SHOULD INVOLVE DIATOMS
			useGroups: (iter of int-iters) Groups to calculate between. The group indices are determined by "classificationOpts". E.g. [ [0,1], [0] ] Would indicate to calculate distribution between groups [0,1] for distrOpts[0] and for group 0 for distrOpts[1]
			useNonHyIdx: (Bool) If True we represent our group with the non-hydrogen indices (if false we use a hydrogen index)
			classificationObjs: (iter of ClassifierBase objects) If present these take priority over classificationOpts. Original purpose was to allow "byReference" classifiers to be used to speed up code when multiple of the same classifiers were used

		"""
		self.fromNonHyIndices = fromNonHyIndices
		self.fromHyIndices = fromHyIndices
		self.classificationOpts = classificationOpts
		self.distrOpts = distrOpts
		self.useGroups = useGroups
		self.useNonHyIdx = useNonHyIdx
		self.classificationObjs = classificationObjs
		self.useIdxEach = 0 #Needed to reuse some "GenericNonHyAndHyFilteredOptsObj_simple" interfaces

	@property
	def primaryIndices(self):
		if self.useNonHyIdx:
			outVals = [x[0] for x in self.fromNonHyIndices]
		else:
			outVals = [x[0] for x in self.fromHyIndices]

		return outVals

	@property
	def atomIndices(self):
		return self.primaryIndices

	@property
	def binResObj(self):
		return [x.binResObj for x in self.distrOpts]





class WaterDerivativeFilteredOptsObj_simple():
	""" Class used to get properties for water-derivative molecules (e.g. hydroxyl/hydronium). Original use was finding planar positions of products for H2O->OH- + H+

	"""

	def __init__(self, oxyIndices, hyIndices, classificationOpts, distrOpts, useGroups, useNonHyIdx=True, useIdxEach=0):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints)
			hyIndices: (iter of ints)
			classificationOpts: (Options object for classification) This contains options for how we classify the indices into N-groups 
			distrOpts: (iter of CalcDistribOptionsBase) Each defines a distribution we want calculated using filtered list of oxyIndices
			useGroups: (iter of int-iters) Groups to calculate between. The group indices are determined by "classificationOpts". E.g. [ [0,1], [0] ] Would indicate to calculate distribution between groups [0,1] for distrOpts[0] and for group 0 for distrOpts[1]
			useNonHyIdx: (Bool) If True we represent our group with one of the non-hydrogen indices (if false we use a hydrogen index)
			useIdxEach: (Int) The index to use in the list of hy/nonHyIndices.
				 
		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.classificationOpts = classificationOpts
		self.distrOpts = distrOpts
		self.useGroups = useGroups
		self.useNonHyIdx = useNonHyIdx
		self.useIdxEach = useIdxEach

	#Used so it still works with atom-combo distribution stuff.
	@property
	def primaryIndices(self):
		return self.oxyIndices

	@property
	def binResObj(self):
		return [x.binResObj for x in self.distrOpts]

