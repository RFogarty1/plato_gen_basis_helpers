
import copy
import types
from . import parse_from_geoms as parseFromGeomBase

from ..shared import cart_coord_utils as cartHelp

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
		distMatrixCentral = cartHelp._getDistMatrixForSetOfCoords(centralCoords)
		distMatrixImageImage = cartHelp._getDistMatrixForSetOfCoords(imageCoords)
		distMatrixCentralAndImages = cartHelp._getDistMatrixBetweenTwoSetsOfSeparateCoords(centralCoords, imageCoords)
	

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
					if (distMatrixCentral[rowIdx][colIdx]<=self.maxBondLength):
						if (distMatrixCentral[rowIdx][colIdx]>=self.minBondLength):
							h2CentralCentralIdxPairs.append( [rowIdx,colIdx] )
			for colIdx in imageH:
				if colIdx > rowIdx:
					if (distMatrixCentralAndImages[rowIdx][colIdx]<=self.maxBondLength):
						if (distMatrixCentralAndImages[rowIdx][colIdx]>=self.minBondLength):
							h2CentralImageIdxPairs.append( [rowIdx,colIdx] )


		#Step 2 = filter this list based on self.minDistOtherNebs
		outCentralCentralPairs, outCentralImagePairs = list(), list()
		for centPair in h2CentralCentralIdxPairs:
			idxA,idxB = centPair
			otherNebsA  = [idx for idx,dist in enumerate(distMatrixCentral[idxA]) if (idx!=idxA) and (idx!=idxB) and (dist<self.minDistOtherNebs)]
			otherNebsA += [idx for idx,dist in enumerate(distMatrixCentralAndImages[idxA]) if (dist<self.minDistOtherNebs)]
			otherNebsB  = [idx for idx,dist in enumerate(distMatrixCentral[idxB]) if (idx!=idxA) and (idx!=idxB) and (dist<self.minDistOtherNebs)]
			otherNebsB += [idx for idx,dist in enumerate(distMatrixCentralAndImages[idxB]) if (dist<self.minDistOtherNebs)]
			if len(otherNebsA)+len(otherNebsB)==0:
				outCentralCentralPairs.append(centPair)

		for centPair in h2CentralImageIdxPairs:
			idxA, idxB = centPair
			otherNebsA  = [idx for idx,dist in enumerate(distMatrixCentral[idxA]) if (idx!=idxA)  and (dist<self.minDistOtherNebs)]
			otherNebsA += [idx for idx,dist in enumerate(distMatrixCentralAndImages[idxA]) if (dist<self.minDistOtherNebs) and idx!=idxB]
			otherNebsB  = [idx for idx,dist in enumerate(distMatrixImageImage[idxB]) if (idx!=idxB) and (dist<self.minDistOtherNebs)]
			centImageIndices = [x[idxB] for x in distMatrixCentralAndImages]
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

