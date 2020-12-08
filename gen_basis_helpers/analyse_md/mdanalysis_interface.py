
import numpy as np

import MDAnalysis as mdAnalysisLib

#Very tricky + likely fragile
def getSimpleAtomicUniverseObjFromTrajObj(trajObj):
	topology = getSimpleAtomicTopologyFromTrajObj(trajObj)
	trajectory = getMDAMemoryReaderTrajFromTrajInMemoryObj(trajObj)
	nAtoms = len(trajObj.trajSteps[0].unitCell.cartCoords)
	outObj = mdAnalysisLib.core.universe.Universe.empty( nAtoms )
	outObj.trajectory = trajectory
	outObj._topology = topology
	mdAnalysisLib.core.universe._generate_from_topology(outObj)
	return outObj

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

