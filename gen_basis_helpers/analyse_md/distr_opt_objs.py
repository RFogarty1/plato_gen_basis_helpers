

from . import calc_distrib_core as calcDistrCoreHelp
from . import calc_radial_distrib_impl as calcRadImpl




class CalcRdfOptions(calcDistrCoreHelp.CalcRdfOptions):
	pass


class CalcPlanarDistOptions(calcRadImpl.CalcPlanarRdfOptions):
	pass


class DiscHBondCounterWithOxyDistFilterOptions(calcDistrCoreHelp.CalcDistribOptionsBase):
	""" Object contains options for counting discrete number of hydrogen bonds from one group to another """

	def __init__(self, oxyIndices, hyIndices, distFilterIndices=None, distFilterVals=None, acceptor=True, donor=True, maxOO=3.5, maxAngle=35):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			distFilterIndices: (iter of ints) These are the indices we look at X-O distances for. They are used to divide water molecules into two groups. Default=None; meaning ALL water are placed in both groups (so we count hydrogen bonds between all water)
			distFilterVals: (iter) min/max distances to be in the two groups. [ [0,3], [3,5] ] would mean groupA are 0<=x<3 from distFilterIndices while groupB are 3<=x<5 from them. Default vals are [ [0,1000], [0,1000] ]
			acceptor: (Bool) If True we calculate dists/angles required to count number of groupA acceptors from groupB. Note changing the order of distFilterValues will give the reverse info (groupB acceptors from groupA)
			donor: (Bool) If True we calculate dists/angles required to count number of groupA donors to groupB. Note changing the order of distFilterValues will give the reverse info (groupB donors to groupA)
			maxOO: (float) The maximum O-O distance between two hydrogen-bonded water. Angles are only calculated when this criterion is fulfilled
			maxAngle: (float) The maximum OA-OD-HD angle for a hydrogen bond; OA = acceptor oxygen, OD=Donor oxygen, HD=donor hydrogen
 
		"""
		self.distribKey = "nHbonds"
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.distFilterVals = distFilterVals #Default to [[0,1000],[0,1000]
		self.distFilterIndices = distFilterIndices #Default to None; means dont filter anything
		self.donor = donor
		self.acceptor = acceptor
		self.maxOO = maxOO
		self.maxAngle = maxAngle

	@property
	def primaryIndices(self):
		return self.oxyIndices

