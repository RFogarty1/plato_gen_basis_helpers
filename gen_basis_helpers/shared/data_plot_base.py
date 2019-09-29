
import contextlib

import numpy as np
import matplotlib.pyplot as plt


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


	def __init__(self, **kwargs):
		for key in self.registeredKwargs:
			setattr(self,key,None)

		for key in kwargs:
			if key in self.registeredKwargs:
				setattr(self,key,kwargs[key])
			else:
				raise KeyError("{} is an invalid keyword.\n Available kwargs are {}".format(key , self.registeredKwargs))


	def _updateAttrsFromKwargs(self, **kwargs):
		for key in self.registeredKwargs:
			if key in kwargs.keys():
				setattr(self,key,kwargs[key])

	def _getDictForAllRegisteredAttrKwargs(self):
		outDict = {k:getattr(self,k) for k in self.registeredKwargs}
		return outDict


	def createPlot(self, plotData, **kwargs):
		""" Takes data in plotData argument and creates a plot 
		
		Args:
			plotData: list of input data. Each list entry is one data series (made from two columns). e.g. [ (np.array(xDataA,yDataA)), (np.array(xDataB,yDataB)) ]
			kwargs: keyword-arguments in same format as when assigning to the object attributres (keys in self.registeredKwargs). These are assigned to the object for solely this function call
		Returns
			Handle to the overall figure
		"""
		#Backup the state from before the function call. TODO: Need to use a context manager to do this more safely
		with temporarilySetDataPlotterRegisteredAttrs(self,kwargs):
	
			toPlot = [np.array(x) for x in plotData]
	
			outFig = plt.figure()
			outFig.add_subplot(111)
	
			plotFunct = self._getPlotFunction()


	
			for idx,pData in enumerate(toPlot):
				try:
					currLabel = self.dataLabels[idx]
				except TypeError:
					currLabel = "data set {}".format(idx)


				plotFunct( pData[:,0], pData[:,1], label=currLabel )
	
			if self.xLim is not None:
				plt.xlim(self.xLim)
	
			if self.yLim is not None:
				plt.ylim(self.yLim)
	
			if self.xlabel is not None:
				plt.xlabel(self.xlabel)
	
			if self.ylabel is not None:
				plt.ylabel(self.ylabel)
	
			if self.showTitle is not None:
				plt.title(self.titleStr)

			if self.legend is True:
				plt.legend()

		return outFig


	def _getPlotFunction(self):
		if self.plotFunct is None:
			return plt.plot
		elif self.plotFunct=="scatter":
			return plt.scatter
		elif self.plotFunct == "plot":
			return plt.plot
		else:
			raise ValueError("{} is an invalid value for plotFunct".format(self.plotFunct))










