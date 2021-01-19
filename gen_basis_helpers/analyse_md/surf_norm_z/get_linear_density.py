
from .. import mdanalysis_interface as mdAnalysisInter
from .. import binned_res as binResHelp

import plato_pylib.shared.unit_convs as uConvHelp

import MDAnalysis.analysis.lineardensity as linDenHelp

def getLinearDensityZ(traj, binWidth, groupIndices=None, lenUnitConv=1, massDict=None):
	""" Description of function
	
	Args:
		traj: (TrajectoryInMemory object)
		binWidth: (float) Width of each bin
		groupIndices: (iter of int) Indices of atoms to include
		lenUnitConv: (Optional, float) Unit conversion to apply to lengths (e.g. angstrom to cm will be common)
		massDict: (Optional, dict) Keys are element symbols, values are mass to use
 
	Returns
		outVals: (BinnedResultsStandard) Contains bins used and values obtained
 
	"""

	#TODO: Probably do a try/except for the charges and masses. Possibly reset to original if values wernt set (maybe always reset them...)
	#Get the universe object
	universeObj = mdAnalysisInter.getSimpleAtomicUniverseObjFromTrajObj(traj)
	mdAnalysisInter.addMassesToUniverseObj(universeObj, eleToMassDict=massDict)
	universeObj.add_TopologyAttr("charges") #Anoyingly needed in order to get linear density, since it tries to also 

	#Use MDAnalysis to do the hard part
	if groupIndices is None:
		selectedAtoms = universeObj.select_atoms("all")
	else:
		selectedAtoms = mdAnalysisInter.getSelectAtomsObjFromIndices(universeObj,groupIndices)

	linDenRunner = linDenHelp.LinearDensity(selectedAtoms, binsize=binWidth)
	linDenRunner.run()
	results = linDenRunner.results

	#Figure out the unit conversions
	cmToAng = 1/uConvHelp.ANG_TO_CM
	volConv  = (1/(cmToAng**3))
	volConv *= 1/(lenUnitConv**3)

	#Convert results into our preferred format
	width = linDenRunner.binsize
	linDenVals = [x*volConv for x  in results["z"]["pos"]]
	binCentres = [ (idx*width)+(0.5*width) for idx,unused in enumerate(results["z"]["pos"]) ]
	binVals = {"lin_den":linDenVals}
	outObj = binResHelp.BinnedResultsStandard.fromConstantBinWidth(binCentres, width, binVals)

	return outObj
