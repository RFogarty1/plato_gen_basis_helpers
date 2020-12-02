
import copy
import itertools as it
import types

import numpy as np

from . import parse_from_geoms as parseFromGeomBase

from ..shared import plane_equations as planeEqnHelp
from ..shared import cart_coord_utils as cartHelp


class DetectSimpleAtomicAdsorbateFromInpGeom(parseFromGeomBase.AdsorbatesFromInpGeom):

	def __init__(self, eleKey, caseSensitive=True, postFilterFuncts=None):
		""" Initializer
		
		Args:
			eleKey: (str) The element symbol for the adsorbate
			caseSensitive: (Bool) Whether to use case sensitivity when looking for the eleKey symbol
			postFilterFuncts: (iter of f(inpGeom,outAds)->outAds)) These functions are applied IN ORDER after the main function to get adsorbates. The original use was to allow only adsorbates either "above" or "below" the surface to be returned
		"""
		self._eleKey = eleKey #May want to do this for multiple eleKeys later; hence hiding the eleKey variable for now
		self.caseSensitive = caseSensitive
		self.postFilterFuncts = list() if postFilterFuncts is None else list(postFilterFuncts)

	def _getFilteredOutputAdsorbates(self, inpGeom, outAds):
		for funct in self.postFilterFuncts:
			outAds = funct(inpGeom, outAds)
		return outAds

	def getAdsorbateObjsFromInpGeom(self, inpGeom):
		#1) Get the cartesian co-ords
		outCartCoords = list()
		for cartCoord in inpGeom.cartCoords:
			if self.caseSensitive:
				if cartCoord[-1] == self._eleKey:
					outCartCoords.append(cartCoord)
			else:
				if cartCoord[-1].upper() == self._eleKey.upper():
					outCartCoords.append(cartCoord)

		#2) Convert into adsorbate objects
		outAdsObjs = [types.SimpleNamespace(geom=copy.deepcopy([x])) for x in outCartCoords]

		return self._getFilteredOutputAdsorbates(inpGeom, outAdsObjs)


class DetectH2OAdsorbatesFromInpGeomStandard(parseFromGeomBase.AdsorbatesFromInpGeom):

	def __init__(self, minBondLength=0.1, maxBondLength=2.0, postFilterFuncts=None):
		""" Initializer
		
		Args:
			minBondLength: (float) Minimum bondlength between O-H.
			maxBondLength: (float) Maximum separation for O-H to be considered as bonded
			postFilterFuncts: (iter of f(inpGeom,outAds)->outAds)) These functions are applied IN ORDER after the main function to get adsorbates. The original use was to allow only adsorbates either "above" or "below" the surface to be returned
				 
		"""
		self._eqTol = 1e-5
		self.minBondLength = minBondLength
		self.maxBondLength = maxBondLength
		self.postFilterFuncts = list() if postFilterFuncts is None else list(postFilterFuncts)

	def _getFilteredOutputAdsorbates(self, inpGeom, outAds):
		for funct in self.postFilterFuncts:
			outAds = funct(inpGeom, outAds)
		return outAds

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

		return self._getFilteredOutputAdsorbates(inpGeom, outObjs)

	def getCentralCentralAdsorbateIndices(self, useGeom):
		centralCoords, imageCoords = cartHelp._getCentralAndImageCoordsFromInpCell(useGeom)
		distMatrixCentral = cartHelp._getDistMatrixForSetOfCoords(centralCoords)
		centralO = [idx for idx,coords in enumerate(centralCoords) if coords[-1].upper()=="O"]
		centralH = [idx for idx,coords in enumerate(centralCoords) if coords[-1].upper()=="H"]
		outIndices = list()
		for oIdx in centralO:
			centralHNebs = self._getHAtomIndicesWithinDistTol(oIdx, centralH, distMatrixCentral)
			if len(centralHNebs)==2:
				currIndices = [oIdx] + centralHNebs
				outIndices.append(currIndices)
		return outIndices

	def _getHAtomIndicesWithinDistTol(self, atomIdx, hIndices, distMatrix):
		outIndices = list()
		for colIdx,dist in enumerate(distMatrix[atomIdx]):
			if (distMatrix[atomIdx][colIdx] <=self.maxBondLength) and (distMatrix[atomIdx][colIdx] >self.minBondLength):
				if (colIdx in hIndices):
					outIndices.append(colIdx)
		return outIndices


	def __eq__(self, other):
		eqTol = min(self._eqTol,other._eqTol)
		numbAttrs = ["minBondLength","maxBondLength"]

		#Check the numerical attributes are similar
		for attr in numbAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if abs(valA-valB) > eqTol:
				return False

		#Check post-filter functs are the same. NOTE: I dont neccesarily expect these to have
		#proper __eq__ implementations, so if any post-filter functions are present objects will probably
		#just always compare unequal
		if len(self.postFilterFuncts) != len(other.postFilterFuncts):
			return False
		else:
			for fA,fB in zip(self.postFilterFuncts,other.postFilterFuncts):
				if fA!=fB:
					return False

		return True
					

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



class FilterToAtomsAboveOrBelowSurf():

	def __init__(self, surfDetector, top=True):
		""" Initializer. Note: I'm not totally sure this will work if your surface crosses PBCs (which it wont in electrode-centric cells)
		
		Args:
			surfDetector: (SurfaceAtomsFromInpGeom object) Used to figure out which atoms make up the surface for an input geometry
			top: (Bool) True means filter to atoms ABOVE the surface (defined as those with >c than the middle of the surface). False means filter to the opposite (those below surface/not above it)
				 
		"""
		self.surfDetector = surfDetector
		self.top = top

	def __call__(self, inpGeom, outAds):

		#1) Get plane equation for the centre of the surface
		useCell = copy.deepcopy(inpGeom)
		surfCoords = self.surfDetector(useCell)
		useCell.cartCoords = surfCoords
		topSurfPlaneEqn = cartHelp.getPlaneEqnForOuterSurfaceAtoms(useCell)
		bottomSurfPlaneEqn = cartHelp.getPlaneEqnForOuterSurfaceAtoms(useCell, top=False)

		#2) Filter to co-ords either above/below the surface
		filteredObjs = list()
		for adsObj in outAds:
			distFromTopPlane = _getDistOfClosestAdsorbateAtomToPlane(topSurfPlaneEqn, adsObj, inpGeom)
			distFromBottomPlane = _getDistOfClosestAdsorbateAtomToPlane(bottomSurfPlaneEqn, adsObj, inpGeom)
			if distFromTopPlane<distFromBottomPlane:
				if self.top:
					filteredObjs.append(adsObj)
			else:
				if not self.top:
					filteredObjs.append(adsObj)
		return filteredObjs


def _getCentralPlaneEqnFromTopAndBottomSurfPlaneEqns(topPlaneEqn,bottomPlaneEqn):
	paramsTop, paramsBot = topPlaneEqn.coeffs, [x*-1 for x in bottomPlaneEqn.coeffs]

	tol = 1e-3
	for pTop, pBot in zip(paramsTop[:3],paramsBot[:3]):
		if abs(pTop-pBot) > tol:
			raise ValueError("a,b,c (top) = {}: a,b,c(bot) = {}. These SHOULD be equal".format(paramsTop,paramsBot))

	outDValue = (paramsBot[-1] + paramsTop[-1]) / 2
	outParams = [x for x in paramsTop[:3]] + [outDValue]
	return planeEqnHelp.ThreeDimPlaneEquation(*outParams)
	

def _getDistOfClosestAdsorbateAtomToPlane(planeEqn, adsObj, inpCell):
	secondPlaneCoeffs = [-1*x for x in planeEqn.coeffs]
	secondPlaneEqn = planeEqnHelp.ThreeDimPlaneEquation(*secondPlaneCoeffs)
	adsGeoms = [x[:3] for x in adsObj.geom]
	minAdsDist = -1

	for atomGeom in adsGeoms:
		distA = _getClosestDistOfAtomToPlane(planeEqn, atomGeom, inpCell)
		distB = _getClosestDistOfAtomToPlane(secondPlaneEqn, atomGeom, inpCell)
		currAdsDist = distA if distA<distB else distB
		if (currAdsDist<minAdsDist) or (minAdsDist<0):
			minAdsDist = currAdsDist

	return minAdsDist

def _getClosestDistOfAtomToPlane(planeEqn, atomGeom, inpCell):
	surfVector = inpCell.lattVects[-1]
	coordA = atomGeom[:3]
	coordUp = [x+t for x,t in it.zip_longest(coordA,surfVector)]
	coordDn = [x-t for x,t in it.zip_longest(coordA,surfVector)]

	minAdsDist = -1
	for coords in [coordA,coordUp,coordDn]:
		currAdsDist = planeEqn.getDistanceOfPointFromPlane(coords)
		if (currAdsDist<minAdsDist) or (minAdsDist<0):
			minAdsDist = currAdsDist
	return minAdsDist
