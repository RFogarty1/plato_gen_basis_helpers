
import math



def calcStandardErrorOfMeanForUncorrelatedData(stdDev, nSamples):
	""" Gets the standard error of the mean (essentially the standard deviation of the mean value; used to estimate its error margins) assuming sampling points are uncorrelated

	Standard.Error = stdDev / sqrt(nSamples)
	
	Args:
		stdDev: (float) The standard deviation of the sample points
		nSamples: (int) The number of sample points used
			 
	Returns
		standardError: (float)
 
	"""
	return stdDev/math.sqrt(nSamples)


def calcVarianceOfData(inpData, besselCorr=False):
	""" Calculate the variance for a set of data
	
	Args:
		inpData: (iter of floats)
		besselCorr: (Bool) If true we multiply by n/n-1 instead of n; this leads to an unbiased estimate
			 
	Returns
		variance: (float) The variance of the input data
 
	Raises:
		 Errors
	"""
	if besselCorr:
		raise NotImplementedError("")
	mean = sum(inpData)/len(inpData)
	outVal = (1/len(inpData)) * sum([(x-mean)**2 for x in inpData])
	return outVal


def getStatsFromBlockingDataUpToMaxOrder(inpData, maxOrder):
	""" Divides data into blocks of blockSize and calculates mean (should be constant between blocksizes assuming not much trimmed) and standard error. This should help in estimating the standard error for the mean for a correlated data set (see e.g. https://doi.org/10.1063/1.457480 or 10.1002/jcc.20746)
	
	Args:
		inpData: (iter of floats)
		maxOrder: (int) We divide each block into two for each order. So maxOrder=0 just gets stats for input data unmodified; maxOrder=1 includes data where we average each pair of data points; maxOrder=4 we average every 4 data points etc.
			 
	Returns
		outDicts: (iter of dicts) Contains various info for each order. Some shown belowish

	"""
	outDicts = list()
	currOrder = 0
	currBlock = inpData

	while currOrder <= maxOrder:
		#Sort out this order
		currVar = calcVarianceOfData(currBlock) / (len(currBlock)-1)
		currStdDev = math.sqrt(currVar)
		currStdDevError = currStdDev * (1 / math.sqrt( 2*(len(currBlock)-1) ))
		currMean = sum(currBlock)/len(currBlock)
		currDict = {"order":currOrder, "mean":currMean, "mean_std_dev":currStdDev, "mean_std_dev_std_dev":currStdDevError}
		outDicts.append(currDict)

		#Block the data
		currBlock = _getDataDividedIntoTwoBlocks(currBlock)

		#Break if we cant block data further
		if len(currBlock) < 2:
			break
		currOrder += 1
	return outDicts


def _getDataDividedIntoTwoBlocks(inpData):
	outData = list()
	idx = 0
	while idx<len(inpData)-1:
		outData.append( (inpData[idx]+inpData[idx+1])/2 )
		idx += 2
	return outData


#Moving averages: At time of writing these are effectively tested in analyse_thermo
def getSimpleMovingAverage(inpData):
	""" Gets the moving average of inpData at each point
	
	Args:
		inpData: (iter of floats) The input data
			 
	Returns
		movingAvg: (iter of floats). Each data point is sum(prevPoints)/len(prevPoints) [i.e. the moving average]
 
	Raises:
		 Errors
	"""
	currSum = 0
	outVals = list()
	for idx,val in enumerate(inpData, start=1):
		currSum += val
		currVal = currSum / idx
		outVals.append(currVal)
	return outVals


def getCentralWindowAverageFromIter(inpIter, widthEachSide, fillValWhenNotCalc=None):
	""" Gets a central-window moving average for the input data (each mean value is calculated from n values either side) 
	
	Args:
		inpIter: (iter of Numbers) We take the moving averages of these numbers
		widthEachSide: (int) Number of data points to take each side when calculating the mean
		fillValWhenNotCalc: The value we output when we cant calculate a central moving average (e.g. we cant get a moving average for the first or last data points)

	Returns
		outVals: (iter of float) Central moving average. Each data point is an average of 2*widthEachSide + 1 data points
 
	"""
	stack = list()
	outVals = list()
	lenIter = len(inpIter)
	reqStackSize = 2*widthEachSide + 1

	#Initialise our stack
	stack = [x for x in inpIter[:widthEachSide]]

	#Calculate moving averages
	for idx,val in enumerate(inpIter):

		#Deal with the stack
		if len(stack) == reqStackSize:
			stack.pop(0)
		
		if idx < lenIter-widthEachSide:
			stack.append(inpIter[idx+widthEachSide])

		#Calculate moving average
		if len(stack) < reqStackSize:
			outVals.append(fillValWhenNotCalc)
		else:
			outVals.append( sum(stack)/len(stack) )
		
	return outVals




