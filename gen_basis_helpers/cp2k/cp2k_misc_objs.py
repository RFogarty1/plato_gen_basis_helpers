

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
		self.paramFile = None if paramFile is None else paramFile
		self.refFunctional = None if refFunctional is None else refFunctional
		self.printDFTD = None if printDFTD is None else printDFTD

	@property
	def modPyCP2KDict(self):
		outDict = dict()
		for currAttr in self.listedAttrs:
			outDict["grimme_disp_" + currAttr] = getattr(self, currAttr)
		return outDict

