
import numpy as np

class ThermoDataInterface():
	""" Goal  of this object is to store thermodynamic data (e.g. pressure/temperature) for varying steps in the simulation

	"""

	@property
	def props(self):
		""" The properties available. Order not gauranteed to be consistent between calls """
		raise NotImplementedError("")

	def getPropsArray(self, props):
		""" Get an array with desired properties
		
		Args:
			props: (iter of str) Each str should correspond to an entry in props
				 
		Returns
			outArray: nxlen(props) np array. Each column corresponds to data of one attribute in props
	 
		"""
		raise NotImplementedError("")


#This implementation assumes we can just read all the data in one go (i.e. that it fits into memory). However, the interface should be
#general enough for when this isnt the case
class ThermoDataStandard(ThermoDataInterface):

	def __init__(self, propDataDict):
		""" Initializer
		
		Args:
			propDataDict: (dict) Keys are strs like "temp", "pressure" etc. Values are iters containing NUMERIC values. All should be the same length; though this wont be checked
	 
		"""
		self._eqTol = 1e-5
		self.dataDict = dict(propDataDict)


	@classmethod
	def fromStdKwargs(cls, step=None, time=None, temp=None, eTotal=None, eKinetic=None,
	                  ePot=None, pressure=None):
		outDict = {"step":step, "time":time, "temp":temp, "eTotal":eTotal, "eKinetic":eKinetic,
		           "ePot":ePot, "pressure":pressure}

		dictKeys = list(outDict.keys())
		for key in dictKeys:
			if outDict[key] is None:
				outDict.pop(key)

		return cls(outDict)

	@property
	def props(self):
		return self.dataDict.keys()


	def getPropsArray(self, props):
		inpList = [self.dataDict[prop] for prop in props]
		outVals = np.array( (inpList) ).transpose()
		return outVals

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)
		#compare properties, make sure order doesnt matter
		ourProps = self.props
		otherProps = other.props
		if len(ourProps)!=len(otherProps):
			return False

		for prop in ourProps:
			if prop not in otherProps:
				return False

		#Get the numerical values are equivalent
		arrayA = self.getPropsArray(ourProps)
		arrayB = other.getPropsArray(ourProps)
		if np.allclose(arrayA,arrayB, atol=eqTol) is False:
			return False


		return True




