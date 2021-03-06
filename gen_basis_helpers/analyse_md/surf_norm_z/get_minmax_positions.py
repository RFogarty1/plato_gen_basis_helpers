


def getMinMaxZPositionsAtomGroups(traj, groupIndices, minZ=None, maxZ=None):
	""" Returns the minimum and maximum z-positions for a group of atoms in every currStep in traj. Example use-case is seeing evolution of a surface for example
	
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




