
import contextlib
import itertools as it

import numpy as np
import matplotlib.pyplot as plt
from ..shared import misc_utils as misc


def addSetOfKwargDescriptorsToClass(inpCls):
	for attrName in inpCls.registeredKwargs:
		setattr(inpCls, attrName, DescriptorBasic(attrName))
	return inpCls


def _addKwargToSetAndMakeBasicDescriptor(attrName, inpCls, kwargSet):
	lenOrig = len(kwargSet)
	
	kwargSet.add(attrName)
	setAttr(inpCls, attrName, DescriptorBasic(attrName))
	
	if lenOrig!= lenEnd:
		raise ValueError("{} seems to be a duplicated attribute".format(attrName))


class DescriptorBasic():

	def __init__(self, attrName):
		self.name = attrName 
		self.val = None

	def __get__(self, instance, owner):
		return self.val
  
	def __set__(self, instance, value):
		self.val = value


@contextlib.contextmanager
def temporarilySetDataPlotterRegisteredAttrs(instance, modKwargDict):
	registeredKwargs = instance.registeredKwargs
	origKwargDict = {k:getattr(instance,k) for k in instance.registeredKwargs}
	for key in modKwargDict:
		if key in registeredKwargs:
			setattr(instance,key,modKwargDict[key])
		else:
			raise KeyError("{} is not a registeredKwargs. Registered Kwargs are {}".format(key,registeredKwargs))

	try:
		yield
	finally:
		for key in origKwargDict:
			setattr(instance, key, origKwargDict[key])




#@addSetOfKwargDescriptorsToClass
class DataPlotterBase():

	registeredKwargs = set()
	registeredKwargs.add("xLim")
	registeredKwargs.add("yLim")
	registeredKwargs.add("xlabel")
	registeredKwargs.add("ylabel")
	registeredKwargs.add("plotFunct")
	registeredKwargs.add("showTitle")
	registeredKwargs.add("titleStr")
	registeredKwargs.add("show")
	registeredKwargs.add("legend")
	registeredKwargs.add("dataLabels")
	registeredKwargs.add("axHandle") #If none we create a new figure; else we just plot on the provided axis handle
	registeredKwargs.add("data")

	def __init__(self, **kwargs):
		#First initialise all arguments
		for key in self.registeredKwargs:
			setattr(self,key,None)

		#Hook to set any default values for attributes
		self._setDefaultInitAttrs()

		for key in kwargs:
			if key in self.registeredKwargs:
				setattr(self,key,kwargs[key])
			else:
				raise KeyError("{} is an invalid keyword.\n Available kwargs are {}".format(key , self.registeredKwargs))


	#Taken from CreatorWithResetableKwargsTemplate
	def _setDefaultInitAttrs(self):
		""" Hook - overiding is optional. This is used to explicitly set any attributes to default values upon initiatiation. Any values set here will be overwritten if the relevant kwarg is passed to the constructor (including if you explicitly pass None with the kwarg)
		"""
		pass



	def _updateAttrsFromKwargs(self, **kwargs):
		for key in self.registeredKwargs:
			if key in kwargs.keys():
				setattr(self,key,kwargs[key])

	def _getDictForAllRegisteredAttrKwargs(self):
		outDict = {k:getattr(self,k) for k in self.registeredKwargs}
		return outDict



	def createPlot(self, plotData=None, **kwargs):
		""" Takes data in plotData argument and creates a plot 
		
		Args:
			plotData: list of input data. Each list entry is one data series (made from two columns). e.g. [ (np.array(xDataA,yDataA)), (np.array(xDataB,yDataB)) ]
			kwargs: keyword-arguments in same format as when assigning to the object attributres (keys in self.registeredKwargs). These are assigned to the object for solely this function call
		Returns
			Handle to the overall figure
		"""

		with temporarilySetDataPlotterRegisteredAttrs(self,kwargs):


			toPlot = self._getToPlotData(plotData)
			outFig = self._getOutfigHandleAndSetCurrentAxis()
	
			for idx,pData in enumerate(toPlot):
				try:
					currLabel = self.dataLabels[idx]
				except TypeError:
					currLabel = "data set {}".format(idx)

				self._plotSingleDataSeries(pData,currLabel, idx)
				
			if self.xLim is not None:
				plt.xlim(self.xLim)
	
			if self.yLim is not None:
				plt.ylim(self.yLim)
	
			if self.xlabel is not None:
				plt.xlabel(self.xlabel)
	
			if self.ylabel is not None:
				plt.ylabel(self.ylabel)
	
			if (self.showTitle is not None):
				if (self.showTitle is not False):
					plt.title(self.titleStr)

			if self.legend is True:
				plt.legend()

			#Allow toPlot to be accesible when modifying the plot
			with misc.temporarilySetInstanceAttrs(self, {"data":toPlot}):
				self._modifyAxisPlot() #Hook for doing further modifications

		return outFig


	def _getOutfigHandleAndSetCurrentAxis(self):
		if self.axHandle is None:
			outFig = plt.figure()
			outFig.add_subplot(111)
		else:
			outFig = None
			plt.sca(self.axHandle)	
		return outFig

	def _getToPlotData(self, plotData=None):
		if plotData is not None:	
			toPlot = [np.array(x) for x in plotData]
		else:
			if self.data is not None:
				toPlot = [np.array(x) for x in self.data]
			else:
				toPlot = list()
		return toPlot


	def _plotSingleDataSeries(self, inpData, label, idx):
		plotFunct = self._getPlotFunction()
		plotFunct( inpData[:,0], inpData[:,1], label=label )

	def _getPlotFunction(self):
		if self.plotFunct is None:
			return plt.plot
		elif self.plotFunct=="scatter":
			return plt.scatter
		elif self.plotFunct == "plot":
			return plt.plot
		else:
			raise ValueError("{} is an invalid value for plotFunct".format(self.plotFunct))


	def _modifyAxisPlot(self):
		""" Hook - provides a function to overwrite to modify the plot after basic things have been done. The current axis will be set correctly, so you can just use plt.something without specifying the axis to work on """
		pass


#Recommended class to use. Less flexible than Base but has a couple of painful to code features inbuilt 
class DataPlotterStandard(DataPlotterBase):
	""" Class containing a createPlot(self, plotData, **kwargs) used to plot a specific type of graph from data of a specific format

	Attributes:
		registeredKwargs: This contains a set of all keyword arguments associated with the class. These can be passed to the constructor
		                  or createPlot. If passed to createPlot they will only be set temporarily (the class state will be the same before/after the call).
		                  NOTE: Some of these are defined in the class, but some are only added upon initialisation. Hence you should query an INSTANCE rather 
		                  than the class
	"""

	def __init__(self, **kwargs):
		self.registeredKwargs.add("lineStyles")
		self.registeredKwargs.add("lineColors")
		self.registeredKwargs.add("lineMarkers")
		self.registeredKwargs.add("lineMarkerSizes")
		self.registeredKwargs.add("removeXTicks") #Currently takes priority over xTickLabels
		self.registeredKwargs.add("removeYTicks")
		self.registeredKwargs.add("mapPlotDataFunct")
		self.registeredKwargs.add("xTickLabels")
		self.registeredKwargs.add("xTickValues")
		self.registeredKwargs.add("xTickLabelRotation")

		super().__init__(**kwargs)
	
		if self.mapPlotDataFunct is None:
			self.mapPlotDataFunct = lambda x: x
	

	def createPlot(self, plotData=None, **kwargs):
		if plotData is None:
			plotData = self.data
		plotData = self.mapPlotDataFunct(plotData)
		outFig = super().createPlot(plotData, **kwargs)


		self._changeLineStylesIfNeeded(outFig, **kwargs)
		self._changeColorsIfNeeded(outFig, **kwargs)
		self._changeLineMarkersIfNeeded(outFig, **kwargs)
		self._changeLineMarkerSizesIfNeeded(outFig, **kwargs)

		self._sortXTickLabels(plt.gca())

		if self.removeYTicks:
			plt.yticks([])

		if self.removeXTicks:
			plt.xticks([])

		#Changes to lineStyles/lineColors means the legend needs remaking (if it was present)
		if kwargs.get("legend",False):
			plt.legend()
		if "legend" not in kwargs:
			if self.legend:
				plt.legend() 

		return outFig 

	#TODO: Refactor to remove duplication with the change line colors function
	def _changeLineStylesIfNeeded(self, outFig, **kwargs):
		def lineStyleModFunct(inpLine, value):
			inpLine.set_linestyle(value)

		self.changeLineProp(outFig, "lineStyles", lineStyleModFunct, **kwargs)


	def _changeColorsIfNeeded(self, outFig, **kwargs):
		def colorLineModFunct(inpLine, value):
			inpLine.set_color(value)

		self.changeLineProp(outFig, "lineColors", colorLineModFunct, **kwargs)

	def _changeLineMarkersIfNeeded(self, outFig, **kwargs):
		def lineMarkerModFunct(inpLine, value):
			inpLine.set_marker(value)
		
		self.changeLineProp(outFig, "lineMarkers", lineMarkerModFunct, **kwargs)


	def _changeLineMarkerSizesIfNeeded(self, outFig, **kwargs):
		def lineMarkerSizeModFunct(inpLine,value):
			inpLine.set_markersize(value)

		self.changeLineProp(outFig, "lineMarkerSizes", lineMarkerSizeModFunct, **kwargs)

	def changeLineProp(self, outFig, propName, propSetFunct, inclLegend=True, **kwargs):
		""" Changes a property for each line in the plot based on the setter function
		
		Args:
			outFig: Plot of interest (assumed to have 1 axis). This is ONLY used as a method of obtaining the axis of interest, using outFig.get_axes()[0]. If self.axHandle is set then it is used in preference to outFig[Using self.axHandle is RECOMMENDED]. 
			propName: Name of property as stored on dataPlotter (e.g. lineStyles). Used to extract values we want to set lines to
			propSetFunct: Function with interface (inpLine, value). Used to set individual lines to individual values (e.g. the style of 1 line)
			kwargs: Values to temporarily set on object (must be registered Kwargs) 
				 
		"""
		#NOTE: The way axHandle is got, through kwargs rather thab self causes problems for inherited classes. TODO: Sort this
		with misc.fragile(temporarilySetDataPlotterRegisteredAttrs(self,kwargs)):
			axHandle = kwargs.get("axHandle",None) if kwargs.get("axHandle",None) is not None else outFig.get_axes()[0]
			usedProp = getattr(self,propName)
			if usedProp is None:
				raise misc.fragile.Break
			dataLines = axHandle.get_lines()
			linePropVals = it.cycle( usedProp )
			for idx, (handle, propVal) in enumerate( zip(dataLines,linePropVals) ):
				propSetFunct(handle,propVal)

			#Need to modify the property for the line in the legend(assuming one is present)
			if self.legend and inclLegend:
				legLineHandles = axHandle.get_legend().get_lines()
				propSetFunct( legLineHandles[idx], propVal )


	def _sortXTickLabels(self, axis):
		if self.xTickLabels is None:
			return None
		if self.xTickValues is None:
			raise ValueError("xTickValues needs to set if xTickValues is set")

		putXTickLabelsOnPlot(axis, self.xTickValues, self.xTickLabels, rotation=self.xTickLabelRotation)




def putXTickLabelsOnPlot(inpAxis, tickVals, plotLabels, rotation=0):
	if rotation is None:
		rotation = 0

	inpAxis.set_xticks(tickVals)
	inpAxis.set_xticklabels(plotLabels)
	for x in inpAxis.get_xticklabels():
		x.set_rotation(rotation)


