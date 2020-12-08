


import MDAnalysis.analysis.rdf as rdfHelp
from . import mdanalysis_interface as mdAnalysisInter


#Note: The MDAnalysis interface is defined here
#https://docs.mdanalysis.org/2.0.0-dev0/documentation_pages/analysis/rdf.html
def getSimpleEleEleRdf(traj, eleA, eleB, distRange, nBins=None):
	""" Gets radial distribution function between two elements
	
	Args:
		traj: (TrajectoryInMemory object) Contains all information on the trajectory
		eleA: (str) Str representation for the first element
		eleB: (str) Str representation for the second element
		distRange: (len-2 iter) [minDist, maxDist] Defines the range of distances to calculate over
		nBins: (int) The number of bins to use. Default is distRange/10
			 
	Returns
		 rdfResults: (RdfBinnedResultsSimple) Contains the radial distribution function between two elements

	Raises:
		 ValueError: if distRange[-1]>L/2 where L is the shortest lattice parameter. Note: You can go to longer range by creating supercells for each trajectory step, but its unlikely to be a good idea generally

	"""
	nBins = int( (distRange[1]-distRange[0])*10 ) if nBins is None else nBins

	#Check that the range doesnt go beyond L/2, else nonsense results will be given
	lattParams = traj.trajSteps[0].unitCell.getLattParamsList()
	if any([x/2<distRange[-1] for x in lattParams]):
		maxLattParam = max(lattParams)
		raise ValueError("Maximum range for rdf is L/2 ({}) but {} requested".format( maxLattParam/2, distRange[-1]) )

	#Use MDAnalysis to do the hard work
	universeObj = mdAnalysisInter.getSimpleAtomicUniverseObjFromTrajObj(traj)
	groupA = universeObj.select_atoms("name {}".format(eleA))
	groupB = universeObj.select_atoms("name {}".format(eleB))
	exclusionBlock = (1,1) if eleA==eleB else None #Stops (for example) atom0-atom0 being counted
	rdfObj = rdfHelp.InterRDF(groupA, groupB, nbins=nBins, range=distRange, exclusion_block=exclusionBlock)
	output = rdfObj.run()

	#Get the results in the format we use
	dists = output.bins
	binEdges = output.edges
	rdfVals = output.rdf
	counts = output.count

	return RdfBinnedResultsSimple(dists, rdfVals, binEdges=binEdges, counts=counts)

class RdfBinnedResultsSimple():
	""" Representation of radial distribution function results

	"""
	def __init__(self, dists, rdf, binEdges=None, counts=None):
		""" Initializer
		
		Args:
			dists: (iter of floats) Distance from atomA. Generally x-axis values and usually will represent centre of histogram bins. MUST BE IN ORDER
			rdf: (iter of floats) Values for the radial distribution function at each dist.
			binEdges: (iter of floats) Edges of the bins. Length will be len(dists)+1. 
			counts: (iter of floats) Raw counts for each bin, i.e. unnormalised values
 
		"""
		self._eqTol = 1e-5
		self.dists = dists
		self.rdf = rdf
		self.binEdges = binEdges
		self.counts = counts


	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)
		
		nonOptNumbAttrs = ["dists","rdf"]
		optNumbAttrs = ["binEdges", "counts"]

		for attr in nonOptNumbAttrs:
			valsA, valsB = getattr(self,attr), getattr(other,attr)
			if len(valsA)!=len(valsB):
				return False
			for vA, vB in zip(valsA, valsB):
				if abs(vA-vB) > eqTol:
					return False

		for attr in optNumbAttrs:
			valsA, valsB = getattr(self,attr), getattr(other,attr)
			if (valsA is None and valsB is None):
				pass
			elif (valsA is None) or (valsB is None):
				return False
			else:
				if len(valsA)!=len(valsB):
					return False
				for vA, vB in zip(valsA, valsB):
					if abs(vA-vB) > eqTol:
						return False

		return True

