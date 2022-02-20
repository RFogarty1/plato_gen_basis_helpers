
import itertools as it

from .. import get_indices_from_geom_impl as getIndicesImplHelp


def getMgIndicesByLayer(inpGeom, eleKey="Mg", distTol=1, nLayers=6, expNumbPerLayer=36):
	""" Gets Mg atoms grouped by surface layers, with the top layer in position [0]
	
	Args:
		inpGeom: (plato_pylib UnitCell object) Contains the geometry
		eleKey: (str) The surface element (Default="Mg")
		distTol: (float) Tolerance in z for determining which layer an atom is in (default=1)
		nLayers: (int) Expected number of layers (default=6)
		expNumbPerLayer: (int) Expected number of Mg per surface layer

	Returns
		outIndices: (iter of int-iters) Each element is one layer. Each element of that is one atom index for that layer. e.g. outIndices[1][2] means the 3rd atom in the 2nd layer
	
	Raises:
		AssertionError: If number of atoms in a layer != 36
	"""
	currArgs = [inpGeom, [eleKey]]
	currKwargs = {"distTol":distTol, "nLayers": nLayers}
	outIndices = getIndicesImplHelp.getIndicesGroupedBySurfLayerForFirstNLayers(*currArgs, **currKwargs)

	assert all([len(x)==expNumbPerLayer for x in outIndices])
	return outIndices


def getHydroxylIndicesByElimination(inpGeom, minAngle=90, maxAngle=120, maxOH=1.2, expNumbWater=144, expNumbHydroxyl=36):
	""" Gets hydroxyl indices by process of elimination; if oxygen/hy ISNT part of a water molecule then we assume there part of a hydroxyl. Dangerous to use outside this specific dataset
	
	Args:
		inpGeom: (plato_pylib UnitCell object)
		minAngle: (float) Minimum H-O-H angle for a water molecule
		maxAngle: (float) Maximum H-O-H angle for a water molecule
		maxOH: (float) Maximum O-H distance for a water molecule
		expNumbWater: (int) The number of water molecules you expect to find in your simulation
		expNumbHydroxyl: (int) The number of hydroxyl molecules you expect to find in your simulation

	Returns
		oxyIndices: (int-iter) Each element contains an oxygen index which is part of a hydroxyl
		hyIndices: (int-iter) Each element contains an oxygen index which is part of a hydroxyl
	
	Raises:
		AssertionError: If we find a number of water molecules OR hydroxyl doesnt match that expected

	"""
	waterOxyIndices, waterHyIndices = getWaterOxyAndHyIndicesFromGeom(inpGeom, minAngle=minAngle, maxAngle=maxAngle, maxOH=maxOH, expNumbWater=expNumbWater)
	chainedOxyIndices, chainedHyIndices = [x for x in it.chain(*waterOxyIndices)], [x for x in it.chain(*waterHyIndices)]

	cartCoords = inpGeom.cartCoords
	hydroxylOxyIndices = [idx for idx,coord in enumerate(cartCoords) if coord[-1].upper()=="O" and idx not in chainedOxyIndices]
	hydroxylHyIndices = [idx for idx,coord in enumerate(cartCoords) if coord[-1].upper()=="H" and idx not in chainedHyIndices]

	assert len(hydroxylOxyIndices)==expNumbHydroxyl, "Expected {} hydroxyl; found {}".format( expNumbHydroxyl, len(hydroxylOxyIndices) )
	assert len(hydroxylHyIndices)==expNumbHydroxyl, "Expected {} hydroxyl; found {}".format( expNumbHydroxyl, len(hydroxylHyIndices) )

	return hydroxylOxyIndices, hydroxylHyIndices

def getWaterOxyAndHyIndicesFromGeom(inpGeom, minAngle=90, maxAngle=120, maxOH=1.2, expNumbWater=144):
	""" Get indices for water separated into oxy/hy indices from input geometry
	
	Args:
		inpGeom: (plato_pylib UnitCell object) Geometry
		minAngle: (float) Minimum H-O-H angle for a water molecule
		maxAngle: (float) Maximum H-O-H angle for a water molecule
		maxOH: (float) Maximum O-H distance for a water molecule
		expNumbWater: (int) The number of water molecules you expect to find in your simulation

	Returns
		oxyIndices: (iter of len-1 iters) Each element contains the oxygen index for one water
		hyIndices: (iter of len-2 iters) Each element contains the hydrogen indices for one water
	
	Raises:
		AssertionError: If we find a number of water molecules that doesnt match that expected
	"""
	#1) Get water indices
	waterIdxGetter = getIndicesImplHelp.GetWaterMoleculeIndicesFromGeomStandard(minAngle=minAngle, maxAngle=maxAngle, maxOH=maxOH)
	waterIndices = waterIdxGetter.getIndicesFromInpGeom(inpGeom)
	assert len(waterIndices)==expNumbWater, "Expected {} water molecules; found {}".format(expNumbWater,len(waterIndices))

	#2) Convert into separate oxy/hy indices
	cartCoords = inpGeom.cartCoords
	oxyIndices, hyIndices = list(), list()

	for currIndices in waterIndices:
		currOxyList, currHList = list(), list()
		for atomIdx in currIndices:
			if cartCoords[atomIdx][-1].upper()=="O":
				currOxyList.append(atomIdx)
			if cartCoords[atomIdx][-1].upper()=="H":
				currHList.append(atomIdx)
		oxyIndices.append(currOxyList)
		hyIndices.append(currHList)

	return oxyIndices, hyIndices





