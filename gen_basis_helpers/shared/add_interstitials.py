
import copy
import itertools as it
import math

import plato_pylib.utils.supercell as supCellHelp

def addSingleInterToHcpBulkGeom(bulkCell, site, ele=None, strat=None):
	""" Add an interstitial to a hcp bulk cell. WARNING: Only tested for perfect elemental structures
	
	Args:
		bulkCell: plato_pylib UnitCell object, containing a hcp geometry
		site: str, denotes position interstitial should be added. e.g. "tetrahedral", "octahedral", "basal_octahedral"
		ele: str (optional), Element for the interstitial; if left blank the function will use the first element it encounters in the structure
		strat: str (optional), Str denoting the strategy to use when deciding where to put the interstitial. Currently only option is "utest"; which puts it in a place which makes it easy to check for the test code
 
	Returns
		None; works in place
 
	Raises:
		ValueError: If bulkCell angles inconsistent with a hcp cell
	"""
	_checkInpCellAnglesConsistentWithHcp(bulkCell)
	if ele is None:
		ele = bulkCell.fractCoords[0][-1]

	if site.lower() == "tetrahedral":
		_addSingleTetrahedralInterstitialToHcpBulkGeom(bulkCell, ele, strat)
	elif site.lower() == "basal_tetrahedral":
		_addSingleBasalTetrahedralInterstitialToHcpBulkGeom(bulkCell, ele, strat)
	elif site.lower() == "octahedral":
		_addSingleOctahedralInterstitialToHcpBulkGeom(bulkCell, ele, strat)
	elif site.lower() == "basal_octahedral":
		_addSingleBasalOctahedralInterstitialToHcpBulkGeom(bulkCell, ele, strat)
	elif site.lower() == "basal_crowdion":
		_addSingleBasalCrowdionInterstitialToHcpBulkGeom(bulkCell, ele, strat)
	elif site.lower() == "basal_split":
		_addSingleBasalSplitInterstitialToHcpBulkGeom(bulkCell, ele, strat)
	else:
		raise ValueError("{} is an invalid value for site variable".format(site.lower()))

def _checkInpCellAnglesConsistentWithHcp(inpCell):
	expAngles = [90,90,120]
	actAngles = inpCell.getLattAnglesList()
	diffs = [abs(exp-act) for exp,act in it.zip_longest(expAngles,actAngles)]
	if not all([x<1e-3 for x in diffs]):
		raise ValueError("Angles {} expected but found {}".format(expAngles,actAngles))

def _addSingleTetrahedralInterstitialToHcpBulkGeom(bulkCell, ele, strat=None):
	#Figure out the atom we centre around
	topAtomIdx, topAtomCoord = _getAtomIdxAndCoordsOfAtomToCentreAround(bulkCell, strat)

	
	nearestOutOfPlaneCoords = _getNearestNebCoordsOutOfZPlane(bulkCell,topAtomIdx) #0 is atom idx
	vectorToNearestOOPNeighbour = [x-y for x,y in it.zip_longest( nearestOutOfPlaneCoords, topAtomCoord )]

	nearestNebDistance =  _getDistTwoVectors([0,0,0], vectorToNearestOOPNeighbour)
	cVector = [0,0,1]
	vectAngle = _getAngleTwoVectors(cVector, vectorToNearestOOPNeighbour)
	zDisp = (0.5*nearestNebDistance) / math.cos(math.radians(vectAngle))

	#Get co-ords for the new atom
	newCoord = copy.deepcopy(topAtomCoord)
	newCoord[-1] += zDisp
	newCartCoord = newCoord + [ele]

	_addAtomCartCoordsToInpCell(bulkCell, newCartCoord)

def _addSingleBasalTetrahedralInterstitialToHcpBulkGeom(bulkCell, ele, strat=None):
	topAtomIdx, topAtomCoord = _getAtomIdxAndCoordsOfAtomToCentreAround(bulkCell, strat)

	nearestOutOfPlaneCoords = _getNearestNebCoordsOutOfZPlane(bulkCell,topAtomIdx) 
	vectorToNearestOOPNeighbour = [x-y for x,y in it.zip_longest( nearestOutOfPlaneCoords, topAtomCoord )]
	nearestNebDistance =  _getDistTwoVectors([0,0,0], vectorToNearestOOPNeighbour)
	cVector = [0,0,1]
	vectAngle = _getAngleTwoVectors(cVector, vectorToNearestOOPNeighbour)
	zDisp = nearestNebDistance*math.cos(math.radians(vectAngle))

	newCoord = copy.deepcopy(topAtomCoord)
	newCoord[-1] += zDisp
	newCartCoord = newCoord + [ele]

	_addAtomCartCoordsToInpCell(bulkCell, newCartCoord)

def _addSingleOctahedralInterstitialToHcpBulkGeom(bulkCell, ele, strat=None):
	topAtomIdx, topAtomCoord = _getAtomIdxAndCoordsOfAtomToCentreAround(bulkCell, strat)
	outCoord = _getCoordsForSingleOctahedralInterstitialToHcpBulkGeom(bulkCell, topAtomIdx, topAtomCoord)
	_addAtomCartCoordsToInpCell(bulkCell, outCoord + [ele])

def _addSingleBasalOctahedralInterstitialToHcpBulkGeom(bulkCell, ele, strat=None):
	topAtomIdx, topAtomCoord = _getAtomIdxAndCoordsOfAtomToCentreAround(bulkCell, strat)
	octaCoord = _getCoordsForSingleOctahedralInterstitialToHcpBulkGeom(bulkCell, topAtomIdx, topAtomCoord)
	zDisp = octaCoord[2] - topAtomCoord[2]
	octaCoord[2] += zDisp #The octahedral displacement was half was we need for the basal octahedral
	_addAtomCartCoordsToInpCell(bulkCell, octaCoord + [ele])

def _addSingleBasalCrowdionInterstitialToHcpBulkGeom(bulkCell, ele, strat=None):
	topAtomIdx, topAtomCoord = _getAtomIdxAndCoordsOfAtomToCentreAround(bulkCell, strat)
	nearestInPlaneCoord = _getNearestNebCoordsInPlane(bulkCell, topAtomIdx)
	vectToNearestInPlane = [b-a for a,b in it.zip_longest(topAtomCoord,nearestInPlaneCoord)]
	interPos = [x+(0.5*d) for x,d in it.zip_longest(topAtomCoord,vectToNearestInPlane)]
	_addAtomCartCoordsToInpCell(bulkCell, interPos + [ele])

def _addSingleBasalSplitInterstitialToHcpBulkGeom(bulkCell, ele, strat=None):
	topAtomIdx, topAtomCoord = _getAtomIdxAndCoordsOfAtomToCentreAround(bulkCell, strat)
	nearestInPlaneCoord = _getNearestNebCoordsInPlane(bulkCell, topAtomIdx, inclImages=False)
	bondVector = [b-a for a,b in it.zip_longest(topAtomCoord, nearestInPlaneCoord)]
	nnDist = _getLenOneVector(bondVector)
	newBondDist = (2*nnDist)/3
	displacement = nnDist - newBondDist 
	displaceVector = [displacement*(x/nnDist) for x in bondVector]
	newAtomCoord = [x+d for x,d in it.zip_longest(topAtomCoord,displaceVector)]
	
	#New displace one atom in place and add the other
	cartCoords = bulkCell.cartCoords
	cartCoords[topAtomIdx] = [x-d for x,d in it.zip_longest(topAtomCoord,displaceVector)] + [cartCoords[topAtomIdx][-1]]
	cartCoords.append(newAtomCoord + [ele])
	bulkCell.cartCoords = cartCoords

def _getCoordsForSingleOctahedralInterstitialToHcpBulkGeom(bulkCell, topAtomIdx, topAtomCoord):

	#Need to get 3 atoms in plane which form a triangle [i assume these atoms are exactly in plane with the first atom]
	nearestInPlaneCoords = _getNearestNebCoordsInPlane(bulkCell, topAtomIdx)
	vectToNearestInPlane = [x-y for x,y in it.zip_longest(nearestInPlaneCoords,topAtomCoord)]
	lenNearestNeb = _getLenOneVector(vectToNearestInPlane)

	#Get atomC coordinates by rotating the x/y directions of AB vector by 60 degrees (using a standard rotation matrix operator)
	aToCVector = [0,0,0]
	aToCVector[0] = vectToNearestInPlane[0]*math.cos(math.radians(60)) - vectToNearestInPlane[1]*math.sin(math.radians(60))
	aToCVector[1] = vectToNearestInPlane[0]*math.sin(math.radians(60)) + vectToNearestInPlane[1]*math.cos(math.radians(60))
	atomCCoords = [a+b for a,b in it.zip_longest(topAtomCoord,aToCVector)]
	atomDCoords = [a+b for a,b in it.zip_longest(atomCCoords,vectToNearestInPlane)]

	#Our parralelogram of atoms can be thought of as two triangles; the centroid of one is directly above
	#an atom in the next plane, while the other is directly above an octahedral interstitial site (i.e. the one we want)
	centroidA = [(a+b+c)/3 for a,b,c in it.zip_longest(topAtomCoord, nearestInPlaneCoords, atomCCoords)]
	centroidB = [(a+b+c)/3 for a,b,c in it.zip_longest(nearestInPlaneCoords, atomCCoords, atomDCoords)]
	nearestCentrA = _getNearestOutOfPlaneCoordsToPointForInpCell(bulkCell, centroidA)
	nearestCentrB = _getNearestOutOfPlaneCoordsToPointForInpCell(bulkCell, centroidB)
	nearestOOPDistCentroidA = _getDistTwoVectors(centroidA, nearestCentrA)
	nearestOOPDistCentroidB = _getDistTwoVectors(centroidB, nearestCentrB)
	topCentroid = centroidA if nearestOOPDistCentroidA > nearestOOPDistCentroidB else centroidB

	#Now we need to find the centroid of these thre atoms and displace downwards halfway to the next plane
	nearestOutOfPlaneCoords = _getNearestNebCoordsOutOfZPlane(bulkCell,topAtomIdx)
	vectorToNearestOOPNeighbour = [x-y for x,y in it.zip_longest( nearestOutOfPlaneCoords, topAtomCoord )]
	distNearestOutOfPlane = _getLenOneVector( vectorToNearestOOPNeighbour )
	cVector = [0,0,1]
	angleAC = _getAngleTwoVectors( vectorToNearestOOPNeighbour, cVector )
	planeSpacing = math.cos(math.radians(angleAC)) * distNearestOutOfPlane
	zDisp = 0.5*planeSpacing

	#Put all this together to get the next co-ordinate
	outCoord = [x for x in topCentroid]
	outCoord[-1] += zDisp
	return outCoord

def _getAtomIdxAndCoordsOfAtomToCentreAround(bulkCell,strat=None):
	cartCoords = bulkCell.cartCoords
	if strat=="utest":
		topAtomIdx = 0
	else:
		topAtomIdx = _getCentralAtomIdxInInpCell(bulkCell)
	topAtomCoord = cartCoords[topAtomIdx][:3]
	return topAtomIdx, topAtomCoord

def _addAtomCartCoordsToInpCell(inpCell, atomCoord):
	cartCoords = inpCell.cartCoords	
	cartCoords.append(atomCoord)
	inpCell.cartCoords = cartCoords


def _getNearestOutOfPlaneCoordsToPointForInpCell(inpCell, pointCoords, zCartTol=1e-2):
	#Step 1 = create the relevant supercell
	superCell = supCellHelp.superCellFromUCell(inpCell,[3,3,3])
	otherPoints = [x[:3] for x in superCell.cartCoords]
	filteredCoords = _getAllCoordsInDiffZPlaneAsInpCartCoord(pointCoords[:3], otherPoints)
	nearestPointIdx = _getIdxOfNearestPointFromListOfCoords(pointCoords[:3], filteredCoords)
	return filteredCoords[nearestPointIdx]


def _getNearestNebDistanceOutOfZPlane(inpCell, atomIdx, zCartTol=1e-2):
	""" Returns the nearest neighbour distance for a given atom in a unit cell (including periodic images by default) EXCLUDING neighbours which are in the same z-plane (within a tolerance)
	
	Args:
		inpCell: plato_pylib UnitCell object
		atomIdx: int, index of the atom you want to look for nearest neighboru distances for. 0 means the 1st atom in the list from inpCell.fractCoords (or inpCell.cartCoords)
		zFractTol, float (Optional), Maximum Difference in z cartesian co-ordinates which allows two atoms to be thought of as in the same plane (Default = 1e-2) 
	"""
	inpAtomCoords = inpCell.cartCoords[atomIdx][:3]
	coords = _getNearestNebCoordsOutOfZPlane(inpCell, atomIdx, zCartTol)
	dist = _getDistTwoVectors(coords, inpAtomCoords)
	return dist

def _getNearestNebCoordsInPlane(inpCell, atomIdx, zCartTol=1e-1, inclImages=True):
	#Step 1 = create the relevant supercell
	inpCartCoord = inpCell.cartCoords[atomIdx][0:3]
	if inclImages:
		superCell = supCellHelp.superCellFromUCell(inpCell,[3,3,3])
	else:
		superCell = supCellHelp.superCellFromUCell(inpCell,[1,1,1])
	newCartCoords = superCell.cartCoords[atomIdx][0:3]
	coordDiffs = [abs(x-y) for x,y in it.zip_longest(inpCartCoord,newCartCoords)]
	assert all([x<1e-5 for x in coordDiffs])


	#Step 2 = figure out the index (within the NEW supercell)
	filteredCoords = _getAllNeighbourCoordsInSameZPlane(superCell.cartCoords, atomIdx, zCartTol)
	idxInFilteredList = _getIdxOfNearestPointFromListOfCoords(inpCartCoord, filteredCoords)

	#Step 3 = return the co-ordinates for that index
	return filteredCoords[idxInFilteredList]


def _getNearestNebCoordsOutOfZPlane(inpCell, atomIdx, zCartTol=1e-2):

	#Step 1 = create the relevant supercell
	inpCartCoord = inpCell.cartCoords[atomIdx][0:3]
	superCell = supCellHelp.superCellFromUCell(inpCell,[3,3,4])
	newCartCoords = superCell.cartCoords[atomIdx][0:3]
	coordDiffs = [abs(x-y) for x,y in it.zip_longest(inpCartCoord,newCartCoords)]
	assert all([x<1e-5 for x in coordDiffs])

	#Step 2 = figure out the index (within the NEW supercell)
	filteredIndices = _getAllIndicesForNeighboursInDiffZPlane(inpCell.cartCoords, atomIdx, zCartTol) #TODO: Not actually needed
	filteredCoords = _getAllNeighbourCoordsInDiffZPlane(inpCell.cartCoords, atomIdx, zCartTol)
	idxInFilteredList = _getIdxOfNearestPointFromListOfCoords(inpCartCoord, filteredCoords)

	#Step 3 = return the co-ordinates for that index
	return filteredCoords[idxInFilteredList]

def _getAllIndicesForNeighboursInDiffZPlane(inpCoords, atomIdx, zCartTol=1e-2):
	inpCartCoord = inpCoords[atomIdx][:3]

	filteredIndices = list()
	for idx,coord in enumerate(inpCoords):
		xyz = coord[0:3]
		if idx!=atomIdx:
			zDisp = abs(xyz[-1] - inpCartCoord[-1])
			if (zDisp > zCartTol):
				filteredIndices.append(idx)
	return filteredIndices


def _getAllNeighbourCoordsInSameZPlane(inpCoords, atomIdx,zCartTol=1e-1):
	inpCartCoord = inpCoords[atomIdx][:3]

	filteredCoords = list()
	for idx,coord in enumerate(inpCoords):
		xyz = coord[0:3]
		if idx!=atomIdx:
			zDisp = abs(xyz[-1] - inpCartCoord[-1])
			if (zDisp < zCartTol):
				filteredCoords.append(xyz)
	return filteredCoords

def _getAllCoordsInDiffZPlaneAsInpCartCoord(inpCartCoord, allCoords, zCartTol=1e-2):
	inpCoords = allCoords
	filteredCoords = list()
	for idx,coord in enumerate(inpCoords):
		xyz = coord[0:3]
		zDisp = abs(xyz[-1] - inpCartCoord[-1])
		if (zDisp > zCartTol):
			filteredCoords.append(xyz)
	return filteredCoords

def _getAllNeighbourCoordsInDiffZPlane(inpCoords, atomIdx,zCartTol=1e-2):
	inpCartCoord = inpCoords[atomIdx][:3]

	filteredCoords = list()
	for idx,coord in enumerate(inpCoords):
		xyz = coord[0:3]
		if idx!=atomIdx:
			zDisp = abs(xyz[-1] - inpCartCoord[-1])
			if (zDisp > zCartTol):
				filteredCoords.append(xyz)
	return filteredCoords

#Assumes periodicity; guess i could put in a switch for that
#Currently private; since it will likely be moved somewhere else soon
def _getNearestNebDistanceForAtomInUcell(inpCell, atomIdx):
	""" Returns the nearest neighbour distance for a given atom in a unit cell (including periodic images by default)
	
	Args:
		inpCell: plato_pylib UnitCell object
		atomIdx: int, index of the atom you want to look for nearest neighboru distances for. 0 means the 1st atom in the list from inpCell.fractCoords (or inpCell.cartCoords)
 
	"""
	#Create a 2x2x2 supercell and check the atom indices havent been rearanged
	cartCoords = inpCell.cartCoords[atomIdx][0:3]
	superCell = supCellHelp.superCellFromUCell(inpCell,[2,2,2])
	newCartCoords = superCell.cartCoords[atomIdx][0:3]
	coordDiffs = [abs(x-y) for x,y in it.zip_longest(cartCoords,newCartCoords)]
	assert all([x<1e-5 for x in coordDiffs])

	#Figure out all distances
	otherAtomCoords = [superCell.cartCoords[x][0:3] for x in range( len(superCell.cartCoords) ) if x!=atomIdx]
	return _getNearestDistanceToPointFromListOfCoords(cartCoords, otherAtomCoords)

def _getCentralAtomIdxInInpCell(bulkCell):
	lattVects= bulkCell.lattVects
	centralPoint = [0,0,0]
	for lVect in lattVects:
		centralPoint = [x+l for x,l in it.zip_longest(centralPoint,lVect)]
	
	atomCoords = [x[:3] for x in bulkCell.cartCoords]
	return _getIdxOfNearestPointFromListOfCoords(centralPoint, atomCoords)

def _getNearestDistanceToPointFromListOfCoords(inpPoint, otherPoints):
	nearestIdx = _getIdxOfNearestPointFromListOfCoords(inpPoint,otherPoints)
	return _getDistTwoVectors(inpPoint,otherPoints[nearestIdx])
	
def _getIdxOfNearestPointFromListOfCoords(inpPoint, otherPoints):
	minDistance = _getDistTwoVectors(inpPoint, otherPoints[0])
	outIdx = 0
	for idx,x in enumerate(otherPoints):
		currDist = _getDistTwoVectors(inpPoint,x)
		if (currDist < minDistance):
			minDistance = currDist
			outIdx = idx
	return outIdx

def _getLenOneVector(vectA):
	return math.sqrt( sum([x**2 for x in vectA]) )

def _getDistTwoVectors(vectA,vectB):
	sqrDiff = [ (a-b)**2 for a,b in it.zip_longest(vectA,vectB) ]
	return math.sqrt( sum(sqrDiff) )

def _getAngleTwoVectors(vectA,vectB):
	normFactorA = math.sqrt( sum( [x**2 for x in vectA] ) )
	normFactorB = math.sqrt( sum( [x**2 for x in vectB] ) )

	normA = [x/normFactorA for x in vectA]
	normB = [x/normFactorB for x in vectB]

	dotProd = sum( [a*b for a,b in it.zip_longest(normA,normB)] )
	return math.degrees( math.acos(dotProd) )

