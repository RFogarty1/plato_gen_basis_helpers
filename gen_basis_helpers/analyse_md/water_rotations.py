
import itertools as it

import numpy as np

from . import binned_res as binResHelp
from . import calc_dists as calcDistsHelp
from . import calc_distrib_core as calcDistribCoreHelp

from ..shared import simple_vector_maths as vectHelp


class CalcStandardWaterOrientationDistribOptions(calcDistribCoreHelp.CalcDistribOptionsBase):

	def __init__(self, binResObj, waterIndices, angleType="roll", checkEdges=True):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			waterIndices: (iter of len-3 iters) Each element contains the indices of atoms in one water molecule of interest
			angleType: (str - roll, pitch or azimuth) The angle type we want to bin. roll/pitch/azimuth correspond to rotations around standard x,y,z axes respectively
			checkEdges: (Bool) If True raise error if the edges of bins go beyond the domains of the angleType

		Raises:
			ValueError: If checkEdges=True and the bin edges go beyond angle domains
		"""
		self.distribKey = "adf"
		self.domainTol = 1 #Allow bins to be up to 1 degrees outside the domain
		self.binResObj = binResObj
		self.checkEdges = checkEdges
		self.waterIndices = waterIndices

		#These are probably SLIGHTLY off; e.g. could be -90<x<=90 instead of -90<=x<=90
		self._angleTypeToDomainMap = {"roll":[-90,90], "pitch":[-90,90], "azimuth":[-180,180]}

		self.angleType = angleType #IMPORTANT to do this after setting checkEdges

	@property
	def angleType(self):
		return self._angleType

	@angleType.setter
	def angleType(self, val):
		newDomain = self._angleTypeToDomainMap[val]
		self._checkBinEdgesWithinDomain(newDomain)
		self._angleType = val

	@property
	def domain(self):
		return self._angleTypeToDomainMap[self.angleType]

	def _checkBinEdgesWithinDomain(self, domain):
		binResHelp._checkBinEdgesWithinDomain(self.binResObj, domain, self.domainTol)


def populateWaterOrientationDistribsFromOptionsObjs(inpTraj, optsObjs):
	""" Populates bins in opts objs using inpTraj
	
	Args:
		inpTraj: (TrajectoryInMemory object)
		optsObjs: (iter of CalcStandardWaterOrientationDistribOptions objects)
			 
	Returns
		Nothing; works in place
 
	"""
	#Get relevant iters out
	binObjs = [x.binResObj for x in optsObjs]
	domains = [x.domain for x in optsObjs]
	waterIndices = [x.waterIndices for x in optsObjs]
	angleTypes = [x.angleType for x in optsObjs]

	#Create single binners
	angleTypeToIdx = {"roll":0, "pitch":1, "azimuth":2}
	singleBinners = list()
	for binObj, indices, angleType in it.zip_longest(binObjs, waterIndices, angleTypes):
		currKwargs = {"resBins":binObj, "indices":indices, "angleIdx":angleTypeToIdx[angleType]}
		currBinner = _WaterRotationAngleBinnerFixedIndices(**currKwargs)
		singleBinners.append( currBinner )

	#Create multi binners
	multiBinnerObj = _WaterRotationAngleMultiBinnerFixedIndices(singleBinners)

	#bin
	nSteps = len(inpTraj.trajSteps)
	for currStep in inpTraj:
		multiBinnerObj.updateCountsFromTrajStep(currStep)

	#Get the pdf and adf values
	for binObj,indices,domain in it.zip_longest(binObjs,waterIndices,domains):
		nAngles = len(indices)*nSteps
		calcDistribCoreHelp._addPdfAndAdfToBinObj(binObj, nAngles, domain)



class _WaterRotationAngleMultiBinnerFixedIndices():

	def __init__(self, singleBinners):
		self.singleBinners = singleBinners

	def updateCountsFromTrajStep(self, trajStep):
		#Get all the UNIQUE sets of indices, then figure out how the non-unique list maps to this
		allIndices = [idxList for idxList in it.chain(*[x.indices for x in self.singleBinners])]
		uniqueIndices = [ list(a) for a in set([tuple(x) for x in allIndices]) ]
		indicesToUniqueIndicesMap = _getMapFromIterToUniqueVals(allIndices, uniqueIndices)


		#OLD (...and faster) way
		#Calculate all angles required (and no more)
		allAngles = getWaterStandardRotationAnglesForInpCell(trajStep.unitCell, uniqueIndices)

#		#divide + conquer-like way (weirdly slower....)
#		allAngles = list()
#		for uWaterIdxs in uniqueIndices:
#			currAngles = getWaterStandardRotationAnglesForInpCell(trajStep.unitCell, [uWaterIdxs])
#			assert len(currAngles)==1
#			allAngles.append(currAngles[0])

		#Pass calculated angles to the singleBinners of interest
		totalIdx = 0
		for binner in self.singleBinners:
			nIndices = len(binner.indices)
			currAngles = [allAngles[indicesToUniqueIndicesMap[idx]] for idx in range(totalIdx,totalIdx+nIndices)]
			binner.updateCountsFromCalcdAngles(currAngles)
			totalIdx += nIndices

def _getMapFromIterToUniqueVals(duplIter, uniqueIndices):
	""" Assuming the same set of values are in both iterA and iterB; this will loop duplIter and for each element get the position of the relevant one in uniqueIndices. 
	
	Args:
		duplIter: (iter) This contains a group of values, where some may be duplicated
		uniqueIndices: (iter) This contains a group of values with no duplications; all values in duplIter should be here
 
	Returns
		 outMap: (len-n iter of ints) n is len(duplIter) For each element in duplIter, this tells us the relevant index in uniqueIndices
 
	NOTES:
		We use equality to map between iters; so this will only work for iters with elements where equality works (i.e. not floats; since they have numerical imprecision which has be accounted for when using equality)

	"""
	outMap = list()
	for val in duplIter:
		outIdx = uniqueIndices.index(val)
		outMap.append(outIdx)
	return outMap

class _WaterRotationAngleBinnerFixedIndices():

	def __init__(self, resBins=None, indices=None, angleIdx=None):
		""" Initializer
		
		Args:
			resBins: (BinnedResultsStandard)
			indices: (iter of len-3 iters)
			angleIdx: (int) Either 0, 1 or 2. These refer to thetaX/thetaY/thetaZ in our standard angle system respectively
				 
		"""
		self.resBins = resBins
		self.indices = indices
		self.angleIdx = angleIdx

	def updateCountsFromCalcdAngles(self, calcdAngles):
		relAngles = [x[self.angleIdx] for x in calcdAngles]
		binResHelp.binCountsFromOneDimDataSimple(relAngles, self.resBins)


def getWaterStandardRotationAnglesForInpCell(inpCell, waterIndices):
	""" Gets the standard rotational angles for water molecules in inpCell. These are the Tait-Bryan angles [theta_x,theta_y,theta_z] / [roll, pitch, azimuth] where rotations are around the axes [1,0,0], [0,-1,0], [0,0,1]
	
	Args:
		inpCell: (plato_pylib unitCell)
		waterIndices: (iter of len-3 iters) Each element contains indices of 3 atoms making up a water molecule
			 
	Returns
		outVals: (iter of len-3 iters) Each contains [theta_x, theta_y, theta_z] values
 
	"""
	rotationMatrices = _getWaterStandardRotationMatricesFromInpCell(inpCell, waterIndices)
	outAngles = _getWaterStandardRotationCoordsFromMatrices(rotationMatrices)
	return outAngles


def _getWaterStandardRotationMatricesFromInpCell(inpCell, waterIndices):
	#1) Reorder so we get [O,H,H] in each case
	orderedIndices = list()
	fractCoords = inpCell.fractCoords
	for waterIdxList in waterIndices:
		oIndices = [idx for idx in waterIdxList if fractCoords[idx][-1].upper()=="O"]
		hIndices = [idx for idx in waterIdxList if fractCoords[idx][-1].upper()=="H"]
		assert len(oIndices) == 1
		assert len(hIndices) == 2
		orderedIndices.append( oIndices + hIndices )

	#Commented out code is the original way; was about 2x slower for a real system
#	#X) Get the minimal matrix linking oxygen positions to hydrogen positions; Note this involves computing A LOT of things i dont need
#	#Get oxy/hydrogen indices separately; and get maps to get their idx to position in this list
#	sortedMinimalOxyIndices = sorted(list(set([x[0] for x in orderedIndices])))
#	sortedMinimalHyIndices = sorted(list(set( it.chain(*[[x[1],x[2]] for x in orderedIndices]) )))
#
#	oxyIdxToMatrixRow = {geomIdx:idx for idx,geomIdx in enumerate(sortedMinimalOxyIndices)}
#	hyIdxToMatrixRow  = {geomIdx:idx for idx,geomIdx in enumerate(sortedMinimalHyIndices)}
#
#	ohVectMatrix = calcDistsHelp.getNearestImageVectorMatrixBasic(inpCell, indicesA=sortedMinimalOxyIndices, indicesB=sortedMinimalHyIndices)
#
#	#2) Convert water indices list to [ [OH],[OH] ] UNIT vectors
#	ohVectors = list()
#	for waterIdxList in orderedIndices:
#		oxyGeomIdx, hyGeomIdxA, hyGeomIdxB = waterIdxList
#		oxyRow, hyColA, hyColB = oxyIdxToMatrixRow[oxyGeomIdx], hyIdxToMatrixRow[hyGeomIdxA], hyIdxToMatrixRow[hyGeomIdxB]
#		
#		vectA = vectHelp.getUnitVectorFromInpVector( ohVectMatrix[oxyRow][hyColA] )
#		vectB = vectHelp.getUnitVectorFromInpVector( ohVectMatrix[oxyRow][hyColB] )
#		ohVectors.append( [vectA,vectB] )


	#NOTE: SOME VERY inefficient function reuse in getNearestImageNeb thing
#	#x.2) We just do nearest image co-ords one at a time to get the vectors
	ohVectors = list()
	cartCoords = inpCell.cartCoords
	for waterIdxList in orderedIndices:
		oxyIdx, hyIdxA, hyIdxB = waterIdxList
		oxyCoord = cartCoords[oxyIdx]
		nearestCoordA = calcDistsHelp.getNearestImageNebCoordsBasic(inpCell, oxyCoord, cartCoords[hyIdxA])
		nearestCoordB = calcDistsHelp.getNearestImageNebCoordsBasic(inpCell, oxyCoord, cartCoords[hyIdxB])
		vectA = vectHelp.getUnitVectorFromInpVector([ y-x for x,y in zip(oxyCoord, nearestCoordA[:3])])
		vectB = vectHelp.getUnitVectorFromInpVector([ y-x for x,y in zip(oxyCoord, nearestCoordB[:3])])
		ohVectors.append( [vectA,vectB] )

	#3) Get the rotation matrices from the oh vectors
	outMatrices = list()
	for vectPair in ohVectors:
		currMatrix = _getStandardRotationMatrixFromTwoOHVectors(vectPair[0],vectPair[1])
		outMatrices.append(currMatrix)

	return outMatrices


def _getWaterStandardRotationCoordsFromMatrices(rotMatrices):

	outAngles = list()
	for rotMatrix in rotMatrices:
		currAngles = vectHelp.getStandardRotationAnglesFromRotationMatrix(rotMatrix)
		#Deal with roll; need it to stay in the +-90 degrees domain
		if currAngles[0] > 90:
			currAngles[0] -= 180
		elif currAngles[0] < -90:
			currAngles[0] += 180

		outAngles.append(currAngles)

	return outAngles

def _getStandardRotationMatrixFromTwoOHVectors(vectA, vectB):
	#Get the standard orientation; where the bisector points along the x axis
	vectAngle = vectHelp.getAngleTwoVectors(vectA, vectB)
	stdVectA = [1,0,0]
	stdVectB = [1,0,0]
	stdVectC = [0,0,1]
	vectHelp.applyRotationAroundAxisToInpCoords([stdVectA], [0,0,1], vectAngle*0.5)
	vectHelp.applyRotationAroundAxisToInpCoords([stdVectB], [0,0,1], -1*vectAngle*0.5)

	#The rotation matrix rotates the standard orientation to the actual vectors; thus can simply use matrix inversion to find it
	useVectA, useVectB = (vectHelp.getUnitVectorFromInpVector(x) for x in [vectA,vectB])
	useVectC = np.cross( np.array(useVectA), np.array(useVectB) ) #up/down is arbitrary; "wrong" direction will simply cause roll to be +180

	stdArray = np.array([stdVectA,stdVectB, stdVectC]).transpose() #We want a different vector in each column; hence transpose is correct
	finalArray = np.array([useVectA,useVectB, useVectC]).transpose()
	rotMatrix = finalArray @ np.linalg.inv(stdArray)

	return rotMatrix




