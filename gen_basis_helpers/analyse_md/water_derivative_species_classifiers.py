
import itertools as it

import numpy as np

from . import classifier_objs as clasObjsBase




class _WaterDerivativeDistanceOnlyClassifierGeneric(clasObjsBase.ClassifierBase):

	def __init__(self, oxyIndices, hyIndices, maxOHDist=1.3, nNebs=2, execCount=0):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) All oxygen atoms you want to check for being part of the water derivative
			hyIndices: (iter of ints) All hydrogen atoms you want to check for being part of the water derivative
			maxDistOH: (float) Maximum distance between O-H for the hydrogen to be considered bonded to the OH 
			nNebs: (int) The number of hydrogen that need to be in range (e.g. 3 for hydronium, 2 for water, 1 for hydroxyl)
			execCount: (int) Used to track how many times .classify is called; used as a safety check when using "byReference" classifiers

		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.maxOHDist = maxOHDist
		self.nNebs = nNebs 
		self.execCount = execCount

	def classify(self, sparseMatrixCalculator):
		oxyIndices,hyIndices = self._assignAllRelevantHydrogenToOxygens(sparseMatrixCalculator)

		outIndices = [ list(), list() ]
		for oxyIdxs, hyIdxs in it.zip_longest(oxyIndices, hyIndices):
			if len(hyIdxs)==self.nNebs:
				outIndices[0].append(oxyIdxs)
				outIndices[1].append(hyIdxs)

		self.storedClassifyResult = outIndices
		self.execCount += 1

		return outIndices


	def _assignAllRelevantHydrogenToOxygens(self, sparseMatrixCalculator):
		#Use sorted indices; makes plenty easier
		sortedOxyIndices, sortedHyIndices = sorted(self.oxyIndices), sorted(self.hyIndices)

		outOxyList = [ [x] for x in sortedOxyIndices ]
		outHyList = [ list() for x in sortedOxyIndices ]

		#Assign[Note: This might be very slow....]
		distMatrix = sparseMatrixCalculator.outDict["distMatrix"]
		for hyIdx in self.hyIndices:
			currDists = distMatrix[hyIdx][np.array(self.oxyIndices)] #Worried about this but...
			if np.min(currDists) < self.maxOHDist:
				currIdx = np.argmin(currDists)
				outHyList[currIdx].append(hyIdx)

		return outOxyList, outHyList

