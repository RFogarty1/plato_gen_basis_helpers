

class GrimmeDispersionCorrOptsCP2K():

	def __init__(self, corrType=None, excludeKindsD3=None, paramFile=None, refFunctional=None, printDFTD=True):
		""" Initializer
		
		Args:
			corrType: (str) Maps to CP2K Type keyword; DFTD2, DFTD3, DFTD3(BJ)
			excludeKindsD3: (str-list) List of atom kinds to exclude from the D3 correction
			paramFile: (str, path) The parameter file to use (default="dftd3.dat")
			refFunctional: (str) Corresponds to the reference_functional keyword, should be equal to the xc-functional you use. Makes CP2K use the optimal scaling factors for given functional
			printDFTD: (Bool, default=True) Whether to turn on the flag to print to a *.dftd file

		"""
		self.listedAttrs = ["corrType", "excludeKindsD3", "paramFile", "printDFTD", "refFunctional"]
		self.corrType = None if corrType is None else corrType
		self.excludeKindsD3 = None if excludeKindsD3 is None else excludeKindsD3
		self.paramFile = "dftd3.dat" if paramFile is None else paramFile
		self.refFunctional = None if refFunctional is None else refFunctional
		self.printDFTD = None if printDFTD is None else printDFTD

	@property
	def modPyCP2KDict(self):
		outDict = dict()
		for currAttr in self.listedAttrs:
			outDict["grimme_disp_" + currAttr] = getattr(self, currAttr)
		return outDict


class NonLocalDispersionsCorrOptsCP2K():

	def __init__(self, corrType=None, cutoff=None, kernelFileName=None, verboseOutput=False):
		""" Initializer
		
		Args:
			corrType: (str) "DRSLL", "LMKLL" or "RVV10" dependning on non-local correlation functional
			cutoff: Optionally specify the cutoff value to use (eV)
			kernelFileName: (str) Name of the kernel data file Default="vdW_kernel_table.dat"
			verboseOutput: (Bool,Default=False) Whether to turn on verbose printing
		"""
		self.listedAttrs = ["corrType", "cutoff", "kernelFileName", "verboseOutput"]
		self.corrType = None if corrType is None else corrType
		self._cutoff = None if cutoff is None else cutoff
		self.kernelFileName = "vdW_kernel_table.dat" if kernelFileName is None else kernelFileName
		self.verboseOutput = None if verboseOutput is None else verboseOutput

	@property
	def modPyCP2KDict(self):
		outDict = dict()
		for currAttr in self.listedAttrs:
			outDict["disp_nl_" + currAttr] = getattr(self, currAttr)
		return outDict

	@property
	def cutoff(self):
		if self._cutoff is None:
			return None
		else:
			return "[eV] {}".format(self._cutoff)

	@cutoff.setter
	def cutoff(self, val):
		self._cutoff = val
		

