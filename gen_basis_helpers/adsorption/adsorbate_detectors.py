
import copy
import types

import numpy as np

from . import parse_from_geoms as parseFromGeomBase

from ..shared import cart_coord_utils as cartHelp

class DetectH2OAdsorbatesFromInpGeomStandard(parseFromGeomBase.AdsorbatesFromInpGeom):

	def __init__(self, minBondLength=0.1, maxBondLength=2.0):
		""" Initializer
		
		Args:
			minBondLength: (float) Minimum bondlength between O-H.
			maxBondLength: (float) Maximum separation for O-H to be considered as bonded
				 
		"""
		self.minBondLength = minBondLength
		self.maxBondLength = maxBondLength

	def getAdsorbateObjsFromInpGeom(self, inpGeom):
		#1) Filter out any non O/H atoms, I may have to reverse later if i want to implement minDistOtherNebs like in the H2 detector
		# This step SHOULD speed things up A LOT though for many systems
		useGeom = copy.deepcopy(inpGeom)
		cartCoords = [x for x in inpGeom.cartCoords if x[-1].upper()=="H" or x[-1].upper()=="O"]
		useGeom.cartCoords = cartCoords

		#Get distance matrices
		centralCoords, imageCoords = cartHelp._getCentralAndImageCoordsFromInpCell(useGeom)
		dMatrices = _getDistanceMatricesFromCentralAndImageCoords(centralCoords, imageCoords)
		imageCentralDistMatrix = np.transpose(dMatrices.centralImage)

		#Get indices of each adsorbate based SOLELY on the distance criteria
		centralO = [idx for idx,coords in enumerate(centralCoords) if coords[-1].upper()=="O"]
		imageO = [idx for idx,coords in enumerate(imageCoords) if coords[-1].upper()=="O"]
		centralH = [idx for idx,coords in enumerate(centralCoords) if coords[-1].upper()=="H"]
		imageH   = [idx for idx,coords in enumerate(imageCoords)   if coords[-1].upper()=="H"]
	

		#Step 1) Get data for the CENTRAL oxygen atoms
		outIndices = list() # each element must be [centralIndices,imageIndices]	
		for oIdx in centralO:
			centralHNebs = self._getHAtomIndicesWithinDistTol(oIdx, centralH, dMatrices.centralCentral)
			imageHNebs = self._getHAtomIndicesWithinDistTol(oIdx, imageH, dMatrices.centralImage)
			currIndices = [ [oIdx]+centralHNebs, imageHNebs ]
			if len(currIndices[0]) + len(currIndices[1])==3:
				outIndices.append(currIndices)
	
		#Step 2) Get data for the IMAGE oxygen atoms
		#image neigbhours needed; but len(centralNebs) >=1 required
		for oIdx in imageO:
			centralHNebs = self._getHAtomIndicesWithinDistTol(oIdx, centralH, imageCentralDistMatrix)
			imageHNebs = self._getHAtomIndicesWithinDistTol(oIdx, imageH, dMatrices.imageImage)
			currIndices = [ centralHNebs, [oIdx]+imageHNebs ]
			if len(currIndices[0])>=1:
				if len(currIndices[0]) + len(currIndices[1]) == 3:
					outIndices.append(currIndices)


		#Step 3) Turn the index lists into adsorbate objects
		outObjs = list()
		for outIdxs in outIndices:
			centralIndices,imageIndices = outIdxs
			currCentCoords = [ copy.deepcopy(centralCoords[x]) for x in centralIndices]
			currImageCoords = [ copy.deepcopy(imageCoords[x]) for x in imageIndices ]
			currObj = types.SimpleNamespace( geom=currCentCoords+currImageCoords )
			outObjs.append(currObj)

		return outObjs


	def _getHAtomIndicesWithinDistTol(self, atomIdx, hIndices, distMatrix):
		outIndices = list()
		for colIdx,dist in enumerate(distMatrix[atomIdx]):
			if (distMatrix[atomIdx][colIdx] <=self.maxBondLength) and (distMatrix[atomIdx][colIdx] >self.minBondLength):
				if (colIdx in hIndices):
					outIndices.append(colIdx)
		return outIndices
					


class DetectH2AdsorbatesFromInpGeomStandard(parseFromGeomBase.AdsorbatesFromInpGeom):

	def __init__(self, minBondLength=0.1, maxBondLength=1.6,  minDistOtherNebs=1.42):
		""" Initializer
		
		Args:
			minBondLength: (float) Minimum bondlength between H-H
			maxBondLength: (float) Maximum separation for two H to be considered as bonded
			minDistOtherNebs: (float) H2 will only be considered H2 if each H atom has no other neighbours closer than minDistOtherNebs

		"""
		self.minBondLength = minBondLength
		self.maxBondLength = maxBondLength
		self.minDistOtherNebs = minDistOtherNebs 

	def getAdsorbateObjsFromInpGeom(self, inpGeom):
		#Get distance matrices for central-central and central-image cells
		centralCoords, imageCoords = cartHelp._getCentralAndImageCoordsFromInpCell(inpGeom)
		dMatrices = _getDistanceMatricesFromCentralAndImageCoords(centralCoords, imageCoords)

		#Get indices of h-atoms in central cell
		centralH = [idx for idx,coords in enumerate(centralCoords) if coords[-1].upper()=="H"]
		imageH   = [idx for idx,coords in enumerate(imageCoords)   if coords[-1].upper()=="H"]

		#Step 1 = Get all the pairs of H2 which satisfy distance criteria
		h2CentralImageIdxPairs = list()
		h2CentralCentralIdxPairs = list()
		for rowIdx in centralH:
			centralHNebs, centralHNebs = list(), list()
			for colIdx in centralH:
				if colIdx > rowIdx:
					if (dMatrices.centralCentral[rowIdx][colIdx]<=self.maxBondLength):
						if (dMatrices.centralCentral[rowIdx][colIdx]>=self.minBondLength):
							h2CentralCentralIdxPairs.append( [rowIdx,colIdx] )
			for colIdx in imageH:
				if colIdx > rowIdx:
					if (dMatrices.centralImage[rowIdx][colIdx]<=self.maxBondLength):
						if (dMatrices.centralImage[rowIdx][colIdx]>=self.minBondLength):
							h2CentralImageIdxPairs.append( [rowIdx,colIdx] )


		#Step 2 = filter this list based on self.minDistOtherNebs
		outCentralCentralPairs, outCentralImagePairs = list(), list()
		for centPair in h2CentralCentralIdxPairs:
			idxA,idxB = centPair
			otherNebsA  = [idx for idx,dist in enumerate(dMatrices.centralCentral[idxA]) if (idx!=idxA) and (idx!=idxB) and (dist<self.minDistOtherNebs)]
			otherNebsA += [idx for idx,dist in enumerate(dMatrices.centralImage[idxA]) if (dist<self.minDistOtherNebs)]
			otherNebsB  = [idx for idx,dist in enumerate(dMatrices.centralCentral[idxB]) if (idx!=idxA) and (idx!=idxB) and (dist<self.minDistOtherNebs)]
			otherNebsB += [idx for idx,dist in enumerate(dMatrices.centralImage[idxB]) if (dist<self.minDistOtherNebs)]
			if len(otherNebsA)+len(otherNebsB)==0:
				outCentralCentralPairs.append(centPair)

		for centPair in h2CentralImageIdxPairs:
			idxA, idxB = centPair
			otherNebsA  = [idx for idx,dist in enumerate(dMatrices.centralCentral[idxA]) if (idx!=idxA)  and (dist<self.minDistOtherNebs)]
			otherNebsA += [idx for idx,dist in enumerate(dMatrices.centralImage[idxA]) if (dist<self.minDistOtherNebs) and idx!=idxB]
			otherNebsB  = [idx for idx,dist in enumerate(dMatrices.imageImage[idxB]) if (idx!=idxB) and (dist<self.minDistOtherNebs)]
			centImageIndices = [x[idxB] for x in dMatrices.centralImage]
			otherNebsB += [idx for idx,dist in enumerate(centImageIndices) if (dist<self.minDistOtherNebs) and idx!=idxA]
			if len(otherNebsA)+len(otherNebsB)==0:
				outCentralImagePairs.append(centPair)

		#Step 3 = Convert central index pairs to adsorbate objects
		outObjs = list()
		for centPair in outCentralCentralPairs:
			idxA, idxB = centPair
			currObj  = types.SimpleNamespace( geom=[ copy.deepcopy(centralCoords[idxA]),
			                                         copy.deepcopy(centralCoords[idxB]) ] )
			outObjs.append(currObj)

		#Step 4 = Convert central-image index pairs to adsorbate objects
		for centImagePair in outCentralImagePairs:
			idxA, idxB = centImagePair
			currObj = types.SimpleNamespace( geom=[ copy.deepcopy(centralCoords[idxA]),
			                                        copy.deepcopy(imageCoords[idxB]) ] )
			outObjs.append(currObj)

		return outObjs

def _getDistanceMatricesFromCentralAndImageCoords(centralCoords, imageCoords):
	distMatrixCentral = cartHelp._getDistMatrixForSetOfCoords(centralCoords)
	distMatrixImageImage = cartHelp._getDistMatrixForSetOfCoords(imageCoords)
	distMatrixCentralAndImages = cartHelp._getDistMatrixBetweenTwoSetsOfSeparateCoords(centralCoords, imageCoords)
	return types.SimpleNamespace( centralCentral=distMatrixCentral, imageImage=distMatrixImageImage, centralImage=distMatrixCentralAndImages )



