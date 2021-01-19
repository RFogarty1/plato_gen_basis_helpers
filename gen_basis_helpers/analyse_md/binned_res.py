


class BinnedResultsStandard():
	""" Simple class for storing results of binning

	Attributes:
		binCentres: (iter of floats) Each value represents the centre of a bin
		binEdges: (iter of floats) Each value represents the edge of a bin (left->right). Same ordering as binCentres
		binVals: (dict) Keys describe a binned property (e.g. "rdf") while values are iters of floats, represnting the value of each bin

	"""

	def __init__(self, binCentres, binEdges, binVals):
		""" Initializer
		
		Args:
			binCentres: (iter of floats) Each value represents the centre of a bin
			binEdges: (iter of floats) Each value represents the edge of a bin (left->right). Same ordering as binCentres
			binVals: (dict) Keys describe a binned property (e.g. "rdf") while values are iters of floats, represnting the value of each bin
				 
		"""
		self._eqTol = 1e-5
		self.binCentres = binCentres
		self.binEdges = binEdges
		self.binVals = binVals

	@classmethod
	def fromConstantBinWidth(cls, binCentres, binWidth, binVals):
		binEdges = [x-(0.5*binWidth) for x in binCentres]
		binEdges.append( binCentres[-1] + (0.5*binWidth) )
		return cls(binCentres, binEdges, binVals)

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)

		#Numerical attributes (but not the dictionary vals)
		numbAttrs = ["binCentres", "binEdges"]
		for attr in numbAttrs:
			valsA, valsB = getattr(self, attr), getattr(other,attr)
			if len( valsA ) != len( valsB ):
				return False
			diffVals = [abs(x-y) for x,y in zip(valsA,valsB)]
			if any([x>eqTol for x in diffVals]):
				return False

		#Check the binVals are equal
		binValsA, binValsB = getattr(self, "binVals"), getattr(other, "binVals")
		if binValsA.keys() != binValsB.keys():
			return False

		for key in binValsA.keys():
			valsA, valsB = binValsA[key], binValsB[key]
			if len(valsA) != len(valsB):
				return False
			diffVals = [abs(x-y) for x,y in zip(valsA,valsB)]
			if any([x>eqTol for x in diffVals]):
				return False

		return True





