
import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import MDAnalysis as mdAnalysisLib

#TODO: Maybe uncomment the code, such that masses can easily be added at conversion time
#Very tricky + likely fragile
#def getSimpleAtomicUniverseObjFromTrajObj(trajObj, addMasses=True, eleToMassDict=None):
def getSimpleAtomicUniverseObjFromTrajObj(trajObj):
	""" Gets an MDAnalysis Universe object from a gen_basis_helpers trajectory object. 
	
	Args:
		trajObj: (TrajectoryInMemory instance)
			 
	Returns
		 outObject: (MDAnalysis universe object) This is the central object for the MDAnalysis library
 
	"""

#		addMasses: (Bool) If True then add atomic masses to the output object
#		eleToMassDict: (dict) Keys are element strings, values are the mass to use for that element. Defaults to a sensible dict


	topology = getSimpleAtomicTopologyFromTrajObj(trajObj)
	trajectory = getMDAMemoryReaderTrajFromTrajInMemoryObj(trajObj)
	nAtoms = len(trajObj.trajSteps[0].unitCell.cartCoords)
	outObj = mdAnalysisLib.core.universe.Universe.empty( nAtoms )
	outObj.trajectory = trajectory
	outObj._topology = topology
	mdAnalysisLib.core.universe._generate_from_topology(outObj)
#	if addMasses:
#		addMassesToUniverseObj(outObj, eleToMassDict=eleToMassDict)

	return outObj


def addMassesToUniverseObj(uniObj, eleToMassDict=None):
	""" Adds masses to a universe object using an eleToMassDict
	
	Args:
		eleToMassDict: (dict) Keys are element strings, values are the mass to use for that element. Defaults to a sensible dict
			 
	"""
	eleToMassDict = eleToMassDict if eleToMassDict is not None else uCellHelp.getEleKeyToMassDictStandard()
	useDict = {k.capitalize():val for k,val in eleToMassDict.items()}
	massVals = list()
	for ele in uniObj.atoms.names:
		currMass = useDict[ele.capitalize()]
		massVals.append( currMass )

	uniObj.add_TopologyAttr("masses", values=massVals)


def getSimpleAtomicTopologyFromTrajObj(trajObj):
	""" Gets a Simple MD Analysis topology object from a TrajectoryInMemory object. This contains only atom names
	
	Args:
		trajObj: (TrajectoryInMemory object) Will only actually use the first trajectory step
			 
	Returns
		outObj: (mdAnalysis Topology object) This contains atom names only. 
 
	"""
	eleNames = [x[-1] for x in trajObj.trajSteps[0].unitCell.cartCoords]
	nAtoms = len(eleNames)
	nResidues, nSegments = nAtoms,1
	residueIndices = [x for x in range(nAtoms)]
	segmentIndices = [0 for x in range(nAtoms)]

	namesAttr = mdAnalysisLib.core.topologyattrs.Atomnames(eleNames)

	currKwargs = {"n_res":nResidues, "n_seg":nSegments, "attrs":[namesAttr],
	              "atom_resindex":residueIndices, "residue_segindex":segmentIndices}
	topologyObj = mdAnalysisLib.core.topology.Topology(nAtoms, **currKwargs)
	return topologyObj


def getMDAMemoryReaderTrajFromTrajInMemoryObj(trajObj):
	""" Gets a trajectory reader object (used in MD Analysis code) from a TrajectoryInMemory object
	
	Args:
		trajObj: (TrajectoryInMemory object)
			 
	Returns
		outObj: (mdanalysis MemoryReader object)
 
	"""
	outGeomArrays = list()
	outGeomDims = list()
	for step in trajObj:
		currCart = np.array([x[:3] for x in step.unitCell.cartCoords])
		currDims = step.unitCell.getLattParamsList() + step.unitCell.getLattAnglesList()
		outGeomArrays.append(currCart)
		outGeomDims.append(currDims)

	outObj = mdAnalysisLib.coordinates.memory.MemoryReader( np.array(outGeomArrays), dimensions=outGeomDims )
	return outObj


def getSelectAtomsObjFromIndices(uniObj, atomIndices):
	selectArgs = list()
	for idx in atomIndices:
		currComm = "index {}".format(idx)
		selectArgs.append(currComm)
	return uniObj.select_atoms(*selectArgs)
