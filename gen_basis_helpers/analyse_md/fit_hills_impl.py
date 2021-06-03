
import copy
import itertools as it

import scipy.optimize

from . import analyse_metadyn_hills as aMetaHillsHelp


def fitHillsToOneDimDataSimple(posVsPot, nGaus=1, scaleAll=1):
	""" Simple function for getting metadynamics hills to fit to a potential
	
	Args:
		posVsPot: (iter of len-2 iters) Positions vs potential, which we want to fit to
		nGaus: (int, Optional) Number of Gaussians to fit
		scaleAll: (float, Optional) Width parameter for the Gaussian hills
 
	Returns
		fitObj: Data structure describing aspects of the fit
		outHills: (MetadynHillsInfo) Contains
 
	"""
	#Get start params 
	inpX = [ x[0] for x in posVsPot]
	targY = [ x[1] for x in posVsPot ]
	scalesAll = [scaleAll for unused in range(nGaus)]
	startParams = getStartParamsToFitSimple(inpX, targY, nGaus=nGaus, scales=scalesAll)

	#Get objective function
	objFunctToFit = getFitHillsObjFunctSimple(inpX, targY)

	#Fit objective function
	simpleFitA = scipy.optimize.minimize(objFunctToFit,startParams)

	#Output relevant data
	outParams = simpleFitA.x
	outHillsObj = _getHillsInfoObjFromOneDimParams(outParams)

	return simpleFitA, outHillsObj


def getStartParamsToFitSimple(inpX, targY, nGaus=1, scales=None):
	""" Gets some reasonable-ish start parameters for fitting Gaussian functions
	
	Args:
		inpX: (iter of floats) Input x-values
		targY: (iter of floats) Target y-values
		nGaus: (int) Number of Gaussians to use
		scales: (Optional, iter of floats) The scale parameter to use to start with. Default is to use 1 for each hill
 
	Returns
		outParams: (iter of floats) Every 3 elements is one Gaussian hill [posA,scaleA,heightA,posB,scaleB,heightB,....]

	Notes:
		a) Heights are determined as max(Y)/nGaus
		b) Positions are centred at the x-value corresponding to max(Y) 

	"""
	assert len(inpX)==len(targY)
	scales = [1 for x in range(nGaus)] if scales is None else scales
	assert len(scales) == nGaus

	#Get the maximum height and its position
	taggedYData = [ [idx,val] for idx,val in enumerate(targY) ]
	maxHeightIdx = max(  taggedYData, key=lambda x:x[1] )[0]
	heightsAll = targY[maxHeightIdx] / nGaus
	posAll = inpX[maxHeightIdx]

	#Combine all parameters
	outParams = list()
	for gauIdx in range(nGaus):
		currParams = [posAll, scales[gauIdx], heightsAll]
		outParams.extend(currParams)
	return outParams


def getFitHillsObjFunctSimple(inpX, targY, objFunct=None):
	""" Creates an objective function for fitting Gaussian Hills to a target potential
	
	Args:
		inpX: (iter of floats) These are the values to evaluate the Gaussians at
		targY: (iter of floats) The target values for the potential at inpX values
		objFunct: f(actY, targY)->objVal Function to transform target Y and actual Y into an objective function value. Default is to use a mean-squared deviation function ( average( sum((yTarg-yAct**2)) ) )
 
	Returns
		outObjFunct: f(params)->objFunctVal. Maps params to hillsInfoObject then evaluates potential at inpX and returns objFunct(actY,targY) 
 
	"""
	objFunct = _msdObjFunct if objFunct is None else objFunct
	def outObjFunct(params):
		hillsObj = _getHillsInfoObjFromOneDimParams(params)
		groupedGaussians = hillsObj.createGroupedHills()
		actYVals = groupedGaussians.evalFunctAtVals([ [x] for x in inpX ])
		outVal = objFunct(targY,actYVals)
		return outVal
	return outObjFunct


def _msdObjFunct(targVals, actVals):
	sdVals = list()
	for tVal, aVal in it.zip_longest(targVals, actVals):
		currSD = (tVal-aVal)**2
		sdVals.append( currSD )
	return sum(sdVals)/len(sdVals)


def _getParamsFromOneDimHillsInfoObj(hillsInfoObj):
	""" Gets parameters for fitting from a hillsInfoObj; assumes that all Gaussians on the hills info object 
	
	Args:
		hillsInfoObj: (MetadynHillsInfo object)
			 
	Returns:
		outParams: (iter of floats) Order is [posA,scaleA,heightA,posB,scaleB,heightB,....]
 
	"""
	positions, scales, heights = [getattr(hillsInfoObj,key) for key in ["positions","scales","heights"]]
	outVals = list()
	for pos,scale,height in it.zip_longest(positions, scales, heights):
		assert len(pos)==1
		assert len(scale)==1
		assert len(height)==1
		currVals = [pos[0], scale[0], height[0]]
		outVals.extend(currVals)
	return outVals

def _getHillsInfoObjFromOneDimParams(inpParams):
	""" Gets MetadynHillsInfo object from a set of parameters in a defined order
	
	Args:
		inpParams: (iter of floats) Order is [posA,scaleA,heightA,posB,scaleB,heightB,....] 
			 
	Returns
		outHillsObj: (MetadynHillsInfo object) All times are set to zero
 
	"""
	idx = 0
	positions, scales, heights = list(), list(), list()
	while idx < len(inpParams):
		pos, scale, height = inpParams[idx], inpParams[idx+1], inpParams[idx+2]
		positions.append([pos])
		scales.append([scale])
		heights.append([height])
		idx += 3
	times = [0 for x in range(len(positions))]
	outKwargs = {"times":times, "positions":positions, "scales":scales, "heights":heights}
	outObj = aMetaHillsHelp.MetadynHillsInfo(**outKwargs)
	return outObj

def getStandardPreparedPesDataForFitToFlatten(inpData, dataRange):
	""" When given an estimate for a potential-energy surface, this gets the data transformed such that transformed + inpData will add to give a flat potential energy surface. 
	
	Args:
		inpData: (iter of len-2 iters)
		dataRange: (len-2 iters or None) [minX,maxX] This determines the range of x-values to include. If set to None then all data is included. If limited, then values outside will be set to a constant value (the constant value is taken from the closest data points outside the range in inpData )
			 
	Returns
		transformedData: (iter of len-2 iters) The y-values in this + those in inpData should lead to a flat PES (though NOT at energy=0). ALSO attempts to mirror the PES; really just best to look at some output to see the transformations this does
 
	"""
	outData = copy.deepcopy(inpData)
	_setDataToConstantValueOutsideOfRange(outData, dataRange)
	outData = _getPesDataMirroredAroundFirstIdx(outData)
	outData = _getPesDataFlippedWithMinvalZero(outData)
	return outData


def _setDataToConstantValueOutsideOfRange(inpData, dataRange):
	""" Modifies y-values of inpData such that values outside of dataRange are set to a constant value. The constant value is taken as the value from the nearest data point
	
	Args:
		inpData: (iter of len-2 iters) 
		dataRange: (len-2 iter) [minX,maxX] Set either value to None to not 
			 
	Returns
		Nothing; works in place
 
	"""
	minX, maxX = dataRange

	taggedData = [ [idx,data] for idx,data in enumerate(inpData) ]
	if minX is not None:
		closestIdx = sorted(taggedData, key=lambda x: abs(x[1][0]-minX))[0][0]
		constY = inpData[closestIdx][1]
		for idx,data in enumerate(inpData):
			if data[0]<minX:
				inpData[idx][1] = constY

	if maxX is not None:
		closestIdx = sorted(taggedData, key=lambda x: abs(x[1][0]-maxX))[0][0]
		constY = inpData[closestIdx][1]
		for idx, data in enumerate(inpData):
			if data[0]>maxX:
				inpData[idx][1] = constY



def _getPesDataMirroredAroundFirstIdx(inpData):
	""" Gets the PES (Potential Energy Surface) mirrored around the first index; Generally expect that to represent a minimum or maximum in x (or close enough). Since we usually only probe one direction in metadynamics this is a sensible guess for the potential to apply in the other direction
	
	Args:
		inpData: (iter of len-2 iters)
			 
	Returns
		outData: (iter of len-2 iters) Includes data mirrored in x. An example shows this best: inpData=[ [0,1], [1,2], [2,4] ] leads to outData= [ [0,1], [1,2], [2,4], [-1,2], [-2,4] ] 
 
	"""
	outData = copy.deepcopy(inpData)
	mirrorX = inpData[0][0]
	for data in inpData[1:]:
		startDeltaX = data[0] - mirrorX
		newDeltaX = -1*startDeltaX
		currY = data[1]
		currX = mirrorX + newDeltaX
		outData.append( [currX,currY] )

	return outData


def _getPesDataFlippedWithMinvalZero(inpData):
	""" Gets the PES (Potential Energy Surface) data whereby the trough is turned into a peak; this is done by multiplying y by -1 and shifting upwards(or technically could be downwards) such that the the minimum value is zero. This should get the potential that needs to be applied to flatten the PES (WITHOUT adding a negative potential).
	
	Args:
		inpData: (iter of len-2 iters)
			 
	Returns
		outData: (iter of len-2 iters). See description for whats in it

	NOTES:
		This flips around y=0. Which should always be fine
 
	"""
	outData = [ [x,y*-1] for x,y in inpData ]
	shiftY = -1*min([data[1] for data in outData])

	for idx,data in enumerate(outData):
		outData[idx][1] += shiftY

	return outData
		


