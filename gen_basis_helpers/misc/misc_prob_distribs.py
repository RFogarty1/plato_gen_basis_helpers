
import numpy as np
import scipy.special

def getHexagonalNAdsAdjacentProbability(nStartSites, nAdded, nPopulatedAdjacent, nAdjacentSites=6):
	""" Function to calculate the probability of N adsorbates being adjacent to reference site; given that all were added at random. This was originally used for putting the degree of clustering into context for water adsorbed on an Mg surface; unlikely that I'll reuse for anything else
	
	Args:
		nStartSites: (int) The total number of avaiable sites at start of adding adsorbates
		nAdded: (int) Number of adsorbates we add (effectively number of "attempts" to get an adsorbate on an adjacent site)
		nPopulatedAdjacent: (int) We get the probability that the number of occupied adsorption sites will b "nPopulatedAdjacent"
		nAdjacentSites: (int) Number of adjacent sites. Its 6 if we assume hcp(0001) surface with just 1 site occupied at the start
			 
	Returns
		 p: (float) The probability of this occuring
 
	"""
	#1) Check for edge cases (such as nAdded < nAdsOnAdjacentSites)
	if nAdded > nStartSites:
		raise ValueError("Trying to add {} adsorbates to {} sites is impossible".format(nAdded, nStartSites))

	#2) Get the probability of the nAdsOnAdjacentSites all being populated by the initial adsorbates
	probAdjContribs = [(nAdjacentSites-idx)/(nStartSites-idx) for idx in range(nPopulatedAdjacent)]
	probAdjTerm = np.product(probAdjContribs)

	#3) Get the probability of the next lot of added adsorbates occupying non-adjacent sites
	maxIdx = nAdded-nPopulatedAdjacent
	probNonAdjContribs = [ (nStartSites-nAdjacentSites-idx)/(nStartSites-nPopulatedAdjacent-idx) for idx in range(maxIdx)]
	probNonAdjTerm = np.product(probNonAdjContribs)

	#4) Deal with the various combinations of how this could happen
	combinationsTerm = scipy.special.comb(nAdded, nPopulatedAdjacent)

	#5)Combine all and return

	return probAdjTerm*probNonAdjTerm*combinationsTerm

