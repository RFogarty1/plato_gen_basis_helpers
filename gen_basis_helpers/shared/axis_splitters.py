
import copy
import itertools as it

import numpy as np

import matplotlib.pyplot as plt


#TODO: Move this to a different file
class PostPlotFunctionBase():
	""" Base class (PostPlotFunctionBase) representing a function which can override the "_modifyAxisPlot" method of the DataPlotterBase class. Only requirement is that it is a callable object which takes ONLY a DataPlotterBase object as the input arg

	"""

	def __call__(self, plotterInstance):
		raise NotImplementedError("")

class SingleAxisSplitterBase(PostPlotFunctionBase):
	""" Base class (SingleAxisSplitterBase) for objects with the function of splitting a single axis (x or y) into multiple axes

	"""


	@property
	def axType(self):
		""" Type of axis we're splitting; either "x" or "y"
		"""
		return self._axType

	@axType.setter
	def axType(self,val):
		self._axType = val


#Attrs:
#1) positions
#2) splitAxes (to stop infinite recursion)
class SingleAxisSplitterTemplate(SingleAxisSplitterBase):

	def __init__(self, positions, limits, yLabelPos=None, xLabelPos=None, extraKwargDicts=None):
		""" Initializer
		
		Args:
			positions (iter of len-2 iters): Each denotes the fractional position of one of the split axes
			limits (iter of len-2 iters): Each denotes the axis limits of one of the split axes

		"""
		self.positions = positions
		self.limits = limits
		self.yLabelPos = yLabelPos
		self.xLabelPos = xLabelPos
		self.extraKwargDicts = extraKwargDicts
		self.splitAxes = False


	def __call__(self, plotterInstance):
		 #Short-circuit if we've already split the axes; lets us recall plotter w/o infinite recursion
		#As this is called by plotterInstance
		if self.splitAxes is True:
			outVal = None
		else:
			outVal = self._splitAxes(plotterInstance)
		return outVal


	#TODO: This is where the main template code goes
	def _splitAxes(self, plotterInstance):
		origAxis = self._getOrigAxisHandle(plotterInstance)
		newAxes = self._createNewAxesFromOrigHandle(origAxis)
		self._setOrigAxisProps(origAxis, plotterInstance)
		self._addPlotToAxes(newAxes, plotterInstance)
		self._modNewAxesHook(newAxes)
		return [origAxis] + newAxes

	def _getOrigAxisHandle(self, plotterInstance):
		if plotterInstance.axHandle is not None:
			origAx = plotterInstance.axHandle
		else:
			origAx = plt.gca()
		return origAx

	def _getFigureHandleFromAxHandle(self, axHandle):
		return axHandle.figure

	def _createNewAxesFromOrigHandle(self, axHandle):
		figHandle = self._getFigureHandleFromAxHandle(axHandle)
		newAxes = list()
		fullPositions = self._getFullNewAxesPositionsFromOrigAx(axHandle)
		for currPos in fullPositions:
			currAx = figHandle.add_axes(currPos)
			newAxes.append(currAx)
		self.splitAxes = True
		return newAxes

	def _getFullNewAxesPositionsFromOrigAx(self, axHandle):
		newAxesReducedPositions = self._getNewAxesPositions(axHandle)
		origPosition = axHandle.get_position().bounds
		outPositions = list()
		for pos in newAxesReducedPositions:
			currFullPos = self._getFullAxisPosFromReducedPosAndOrigFullPosition(pos, origPosition)
			outPositions.append(currFullPos)
		return outPositions

	def _getNewAxesPositions(self, axHandle):
		startPosition = self._getOrigAxisReducedPositionsFromInpHandle(axHandle)
		outPositions = list()
		for fractPos in self.positions:
			currOutPos = self._getNewPosFromOldAndFractPos(startPosition, fractPos)
			outPositions.append( currOutPos )
		return outPositions

	def _addPlotToAxes(self, axes, plotterInstance):
		kwargDicts = self._getKwargDictsForNewAxes()
		newPlotter = copy.deepcopy(plotterInstance)
		for currAx,kwargDict in it.zip_longest(axes,kwargDicts):
#			print("kwargDict[legend] = {}".format(kwargDict["legend"]))
#			newPlotter.legend = kwargDict["legend"]
			newPlotter.createPlot(axHandle=currAx, **kwargDict)

	def _getNewPosFromOldAndFractPos(self, oldPos, fractPos):
		origStartPos, origLength = oldPos
		newStartPos = origStartPos + (fractPos[0]*origLength)
		newLength = (fractPos[1]-fractPos[0])*origLength
		return [newStartPos, newLength]

	#Fix label
	def _setOrigAxisProps(self, origAxis, plotterInstance):
		#Getting rid of all the borders in the original case
		origAxis.spines['top'].set_visible(False)
		origAxis.spines['right'].set_visible(False)
		origAxis.spines['bottom'].set_visible(False)
		origAxis.spines['left'].set_visible(False)	

		#Make sure the plot data doesnt appear anymore and that we no longer use any of these tick markers
		maxXVal = max( [max(np.array(x[:,0])) for x in plotterInstance.data] )
		origAxis.set_xlim([maxXVal+1, maxXVal+2])
		origAxis.set_yticks([])
		origAxis.set_yticklabels([])
		origAxis.set_xticks([])
		origAxis.set_xticklabels([])
		origAxis.legend().set_visible(False)

    #Want to know which is the furthest down (for y splitting) or furthest to the left (for x splitting) axis
	def _getFirstAxisHandleIdx(self):
		outIdx, outFractVal = 0, self.positions[0][0]
		for currIdx,currPos in enumerate(self.positions):
			if currPos[0] < outFractVal:
				outIdx = currIdx
		return outIdx

	def _setLabelPositions(self, origAxis):
		if self.yLabelPos is not None:
			outAx = origAxis.get_yaxis()
			origY = outAx.get_label()._y
			outAx.set_label_coords(self.yLabelPos, origY)
		if self.xLabelPos is not None:
			outAx = origAxis.get_xaxis()
			origX = outAx.get_label()._x
			outAx.set_label_coords(origX, self.xLabelPos)

		
	#LIKELY TO BE OVERRIDEN
	def _getKwargDictsForNewAxes(self):
		sharedDict = self._getSharedKwargDictForNewAxes()
		return [sharedDict for x in self.positions]

	def _getSharedKwargDictForNewAxes(self):
		outDict = dict()
		outDict["ylabel"] = None
		outDict["xlabel"] = None
		outDict["showTitle"] = False
		outDict["legend"] = False
		return outDict

	def _modNewAxesHook(self, newAxes):
		pass

	#Will differ for x/y axes. reducedPos is [startPos, length] for jsut the axis (x or y) of interest
	# while origFullPos is the full definition of position for the original axis
	def _getFullAxisPosFromReducedPosAndOrigFullPosition(self, reducedPos, origFullPos):
		raise NotImplementedError("")

	def _getOrigAxisReducedPositionsFromInpHandle(self, axHandle):
		raise NotImplementedError("")


class SingleXAxisSplitter(SingleAxisSplitterTemplate):

	@property
	def axType(self):
		return "x"

	def _getFullAxisPosFromReducedPosAndOrigFullPosition(self, reducedPos, origFullPos):
		outPos = list(origFullPos)
		outPos[0], outPos[2] = reducedPos
		return outPos

	def _getOrigAxisReducedPositionsFromInpHandle(self, axHandle):
		fullPos = axHandle.get_position().bounds
		return [fullPos[0], fullPos[2]]


	def _getKwargDictsForNewAxes(self):
		outList = list()
		for currLim in self.limits:
			currDict = self._getSharedKwargDictForNewAxes()
			currDict["xLim"] = currLim
			outList.append( currDict )

		#Update with any user-kwarg dicts
		if self.extraKwargDicts is not None:
			for currDict, currExtra in it.zip_longest(outList,self.extraKwargDicts):
				currDict.update(currExtra)

		return outList

	def _setOrigAxisProps(self, origAxis, plotterInstance):
		super()._setOrigAxisProps(origAxis, plotterInstance)
		origAxis.set_xticks([])
		self._setLabelPositions(origAxis)

	def _modNewAxesHook(self, newAxes):
		firstAxHandleIdx = self._getFirstAxisHandleIdx()
		for idx,x in enumerate(newAxes):
			if idx != firstAxHandleIdx:
				x.set_yticks([])
				x.set_yticklabels([])


class SingleYAxisSplitter(SingleAxisSplitterTemplate):

	@property
	def axType(self):
		return "y"

	def _getFullAxisPosFromReducedPosAndOrigFullPosition(self, reducedPos, origFullPos):
		outPos = list(origFullPos)
		outPos[1], outPos[3] = reducedPos
		return outPos

	def _getOrigAxisReducedPositionsFromInpHandle(self, axHandle):
		fullPos = axHandle.get_position().bounds
		return [fullPos[1], fullPos[3]]

	def _getKwargDictsForNewAxes(self):
		outList = list()
		for currLim in self.limits:
			currDict = self._getSharedKwargDictForNewAxes()
			currDict["yLim"] = currLim
			outList.append( currDict )

		#Update with any user-kwarg dicts
		if self.extraKwargDicts is not None:
			for currDict, currExtra in it.zip_longest(outList,self.extraKwargDicts):
				currDict.update(currExtra)

		return outList

	def _setOrigAxisProps(self, origAxis, plotterInstance):
		super()._setOrigAxisProps(origAxis, plotterInstance)
		origAxis.set_yticks([])
		self._setLabelPositions(origAxis)

	#May not be the best place for this; liekly want usr options to overwrite all else...
	def _modNewAxesHook(self, newAxes):
		firstAxHandleIdx = self._getFirstAxisHandleIdx()
		for idx,x in enumerate(newAxes):
			if idx != firstAxHandleIdx:
				x.set_xticks([])
				x.set_xticklabels([])

			

