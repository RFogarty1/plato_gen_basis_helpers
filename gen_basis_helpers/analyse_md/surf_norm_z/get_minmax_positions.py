



def getLayerThicknessDataFromMinMaxData(minMaxData):
	""" Simple function for getting surface thicknesses from min and max surf positions data
	
	Args:
		minMaxData: (iter of len-2 iter) [min,max]; the minimum and maximum z-val (or any val) bounding a surface
			 
	Returns
		surfThicknesses: (iter of floats) max-min for each value in minMaxData
 
	"""
	outData = list()
	for minVal,maxVal in minMaxData:
		outData.append( maxVal-minVal )
	return outData

def getAverageMinMaxZPositionsForSubsetOfAtomGroups(traj, groupIndices, nMaxZ=None, nMinZ=None):
	""" Gets the average z positions of a group of atom indices BUT only using n-atoms with the lowest/highest z-values (use case is to get an estimate of changing surface height without knowing that the same atoms will always be at the surface) 
	
	Args:
		traj: (TrajectoryBase object)
		groupIndices: (int list) Indices for atoms in this group
		nMaxZ: (int) The maximum number of indices used to calculate the maxZ
		nMinZ: (int) The maximum number of indices used to calculate minZ

	Returns
		outVals: iter of len-2 iters, Each element has [minZ,maxZ]. Id nMaxZ/nMinZ arent set then both values (minZ,maxZ) will be the same
 
	"""
	nMinZ = len(groupIndices) if nMinZ is None else nMinZ
	nMaxZ = len(groupIndices) if nMaxZ is None else nMaxZ

	outVals = list()

	for currStep in traj:
		currCartCoords = currStep.unitCell.cartCoords
		zVals = sorted([currCartCoords[idx][2] for idx in groupIndices])
		revSorted = sorted(zVals, reverse=True)
		avgMin = sum(zVals[:nMinZ]) / nMinZ
		avgMax = sum(revSorted[:nMaxZ]) / nMaxZ
		outVals.append( [avgMin, avgMax] )

	return outVals


def getMinMaxZPositionsAtomGroups(traj, groupIndices, minZ=None, maxZ=None):
	""" Returns the minimum and maximum z-positions for a group of atoms in every currStep in traj. Example use-case is seeing evolution of a surface 
	
	Args:
		traj: (TrajectoryBase object)
		groupIndices: (int list) Indices for atoms in this group
		minZ: (float) The minimum z-value allowed. Useful for looking at regions above/below a surface
		maxZ: (float) The maximum z-value allowed

	Returns
		outVals: iter of len-2 iters, Each element has [minZ,maxZ]
 
	"""
	outVals = list()

	for currStep in traj:
		if (minZ is None) and (maxZ is None):
			currZ = [x[2] for idx,x in enumerate(currStep.unitCell.cartCoords) if idx in groupIndices]
		elif (minZ is not None) and maxZ is None:
			currZ = [x[2] for idx,x in enumerate(currStep.unitCell.cartCoords) if idx in groupIndices and x[2]>minZ]
		elif (minZ is None) and (maxZ is not None):
			currZ = [x[2] for idx,x in enumerate(currStep.unitCell.cartCoords) if idx in groupIndices and x[2]<maxZ]
		elif (minZ is not None) and (maxZ is not None):
			currZ = [x[2] for idx,x in enumerate(currStep.unitCell.cartCoords) if idx in groupIndices and x[2]<maxZ and x[2]>minZ]
		else:
			raise ValueError("This should never happen")


		currVals = [ min(currZ), max(currZ) ]
		outVals.append(currVals)
	return outVals



def getAverageZPositionForAtomGroup(traj, groupIndices):
	""" Returns the average z-position for a group of atoms for every step in traj
	
	Args:
		traj: (TrajectoryBase object)
		groupIndices: (int list) Indices for atoms in this group
			 
	Returns
		avPos: (iter of floats) Each element is the average z-position of the group for one timestep
 
	"""
	outVals = list()

	for currStep in traj:
		zVals = [ x[-2] for idx,x in enumerate(currStep.unitCell.cartCoords) if idx in groupIndices ]
		avZVal = sum(zVals)/len(zVals)
		outVals.append( avZVal )
	return outVals
