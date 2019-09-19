
import numpy as np
import matplotlib.pyplot as plt
from ..shared import data_plot_base as basePlotter



#TODO: Probably put the kwarg registration in the base class file in future, so the dataPlotter interface is pretty fixed
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


#WARNING: THIS IS PROBABLY IMPLEMENTED VERY VERY WRONGLY AND STUPIDLY
class DescriptorBasic():

	def __init__(self, attrName):
		self.name = attrName 
		self.val = None

	def __get__(self, instance, owner):
		return self.val
  
	def __set__(self, instance, value):
		self.val = value



#Class needed to be defined last due to the decorator
@addSetOfKwargDescriptorsToClass
class DataPlotterConvergers(basePlotter.DataPlotterBase):

	registeredKwargs = set()
	registeredKwargs.add("xLim")
	registeredKwargs.add("yLim")
	registeredKwargs.add("xlabel")
	registeredKwargs.add("ylabel")
	registeredKwargs.add("plotFunct")
	registeredKwargs.add("showTitle")
	registeredKwargs.add("titleStr")
	registeredKwargs.add("show")
	def __init__(self, **kwargs):
		for key in self.registeredKwargs:
			if key in kwargs:
				setattr(self,key,kwargs[key])
			else:
				setattr(self,key,None)


	@classmethod
	def fromDefaultPlusKwargs(cls, **kwargs):
		inpKwargs = dict()
		inpKwargs["xlabel"] = "Convergence parameter"
		inpKwargs["ylabel"] = "Property Value / units"
		inpKwargs["showTitle"] = True
		inpKwargs["titleStr"] = "Convergence Plot"
		inpKwargs["plotFunct"] = "scatter"
		inpKwargs.update(kwargs) #user-defined kwargs overwrite defualt ones where relevant
		return cls(**inpKwargs)


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
		startKwargs = self._getDictForAllRegisteredAttrKwargs()
		self._updateAttrsFromKwargs(**kwargs)

		toPlot = [np.array(x) for x in plotData]

		outFig = plt.figure()
		outFig.add_subplot(111)

		plotFunct = self._getPlotFunction()

		for pData in toPlot:
			plotFunct( pData[:,0], pData[:,1] )

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

		#Reset state to before the function call
		self._updateAttrsFromKwargs(**startKwargs)

		return outFig


	def _getPlotFunction(self):
		if self.plotFunct is None:
			return plt.plot
		elif self.plotFunct=="scatter":
			return plt.scatter
		else:
			raise ValueError("{} is an invalid value for plotFunct".format(self.plotFunct))

