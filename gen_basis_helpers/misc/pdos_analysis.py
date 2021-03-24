
import copy
import itertools as it
import sim_xps_spectra.mol_spectra.standard_objs as simXpsStdObjsHelp
import sim_xps_spectra.mol_spectra.spectrum_creator as specCreatorHelp
import sim_xps_spectra.broad_functs.create_broaden_functs as broadFunctHelp



def getGaussianBroadenedTotalDosFromFragment(pdosFragment, xVals, fwhm, multByOcc=False, normFactor=1):
	""" Get a total density of states broadened by gaussian functions
	
	Args:
		pdosFragment: (PdosFragmentStandard object)
		xVals: (iter of floats) The x-values to evaluate at
		fhwm: (float) The width of the Gaussian function
		multByOcc: (Bool) Whether to multiply intensities by occupancies. If True, you may want to set normFactor to 1/max(occ)
		normFactor: (float) Area of the broadening function
 
	Returns
		outSpec: (iter of len-2 iters) Basically and nx2 array; first column is xVals while second is tDoS
 
	"""
	#Create a fragment with a single breakdown
	useFrag = copy.deepcopy(pdosFragment)
	useFrag.breakdowns = [ [1] for x in useFrag.eigenValues ] 
	useFrag.breakdownLabels = ["total"]

	#Create the spectrum
	specObj = getGaussianBroadenedPdosFromFragmentsSimple( [useFrag], ["total_dos"], xVals, fwhm, multByOcc=multByOcc, normFactor=normFactor)

	return specObj.totalSpectralContributions



def getGaussianBroadenedPdosFromFragmentsSimple(pdosFragments, fragNames, xVals, fwhm, multByOcc=False, separateShells=False, normFactor=1):
	""" Gets a set of partial density of states broadened with a Gaussian function
	
	Args:
		pdosFragments: (iter of PdosFragmentStandard objects)
		fragNames: (iter of str describing each pdos fragment)
		xVals: (iter of floats) The x-values to evaluate at
		multByOcc: (Bool) Whether to multiply intensities by occupancies. If True, you may want to set normFactor to 1/max(occ)
		fwhm: (float) The width of the Gaussian function
		separateShells: (Bool) Whether to split the pdos fragments into separate shells (s,p,d) or treat them all together
		normFactor: (float) Area under each broadening function

	Returns
		outSpectrum: (sim_xps GenSpectraOutputStandard object) Contains data on each individual pdos contribution and (less useful here) their sum. Individual contribs can be searched by using their fragName and (if separateShells=True) their breakdownHeaders
 
	"""
	#Step 1 = get the full fragment objects
	allFrags = list()
	for currFrag,fragName in it.zip_longest(pdosFragments,fragNames):
		currOutFrags = getSimXpsSpecFragmentsFromPdosFragmentStandard(currFrag, fragName=fragName, multByOcc=multByOcc, separateShells=separateShells)
		allFrags.extend(currOutFrags)

	#Step 2 = get the broadening function
	bFunct = broadFunctHelp.createNormalisedGauFunctFromCentreAndFWHM(0.0,fwhm,area=normFactor)

	#Step 3 = Generate the spectrum from the fragment objects and the correct inputs
	currKwargs = {"spectraFrags":allFrags, "normBFunct":bFunct, "xVals":xVals}
	specCreator = specCreatorHelp.SpectrumCreatorStandard(**currKwargs)

	return specCreatorHelp.createSpectrumFromStandardCreator(specCreator)



def getSimXpsSpecFragmentsFromPdosFragmentStandard(pdosFragment, fragName=None, multByOcc=False, separateShells=False):
	""" Gets sim_xps* code SpectrumFragmentStandard objects from CP2K PdosFragment objects
	
	Args:
		pdosFragment: (PdosFragmentStandard object) Used in cp2k modules but not neccesarily specific to CP2K
		fragName: (str) A unique name for this fragment. Almost always needs setting
		multByOcc: (Bool) Whether to multiply intensities by occupancies
		separateShells: (Bool) If True then return separate objects for individual shells
			 
	Returns
		frags: (iter of SpectrumFragmentStandard objects) These are key objects in sim_xps* code. Usually len will be one, but maybe more if using separateShells
  
	"""
	fragName = "" if fragName is None else fragName
	energies = pdosFragment.eigenValues
	if separateShells is False:
		outLabels = [ simXpsStdObjsHelp.MolFragLabel(fragKey=fragName, eleKey="", aoKey="") ]
		intensities = [ [sum(x) for x in pdosFragment.breakdowns] ] 
	else:
		intensities = list()
		for idx in range(len(pdosFragment.breakdowns[0])):
			currIntensities = [x[idx] for x in pdosFragment.breakdowns]
			intensities.append(currIntensities)
		outLabels = [simXpsStdObjsHelp.MolFragLabel(fragKey=fragName, eleKey="", aoKey=header) for header in pdosFragment.breakdownHeaders]


	if multByOcc:
		occs = pdosFragment.occs

		for idx,iVals in enumerate(intensities):
			intensities[idx] = [occ*i for occ,i in it.zip_longest(occs,iVals)]

	outObjs = list()
	for label,currIntensities in it.zip_longest(outLabels, intensities):
		currObj = specCreatorHelp.SpectrumFragmentStandard(energies, currIntensities, label)
		outObjs.append(currObj)

	return outObjs

