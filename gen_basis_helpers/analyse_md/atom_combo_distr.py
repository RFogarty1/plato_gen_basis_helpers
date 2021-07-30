
import itertools as it

import numpy as np

from . import binned_res as binResHelp
from . import calc_distrib_core as calcDistrCoreHelp
from . import calc_dists as calcDistsHelp
from . import calc_radial_distrib_impl as calcRadImpl
from . import water_combo_distrs as waterComboDistrHelp

from ..shared import plane_equations as planeEqnHelp


def getAtomicComboDistrBinsFromOptsObjs(inpTraj, optsObjsGroups):
	""" Gets populated NDimensionalBinnedResults objects from a trajectory and options objects
	
	Args:
		inpTraj: (TrajectoryInMemory object)
		optsObjsGroups: (iter of iter of calcOptions objects). Currently these can have "CalcRdfOptions" (with minDistAToB=True) or "" or "CalcPlanarRdfOptions" option objs (a mixture is fine). 
			 
	Returns
		outRes: (iter of NDimensionalBinnedResults) One of these per element in optsObjs
 
	Raises:
		ValueError: Should be raised if first set of indices are not all the same within one element of optsObjsGroups (e.g. all objects in optsObjsGroups[1] must be the same; but they can differ between optsObjsGroups[0] and optsObjsGroups[1]
	"""

	#1) Check options groups are consistent; each should have the same set of indicesA/indices
	for optObjGroup in optsObjsGroups:
		_checkIndicesConsistentInOptsObjGroup(optObjGroup)

	#2) Setup object to calculate basic matrices (e.g. dist matrix) and one to get bin values from this
	sparseMatrixCalculator = _getMatrixCalculatorForOptsObjsGroups(optsObjsGroups)
	binValGetters = list()
	for group in optsObjsGroups:
		currBinValGetter = _GetMultiDimValsToBinFromSparseMatrices.fromOptsObjs(group)
		binValGetters.append( currBinValGetter )

	#3) Setup the bin objects
	outBinObjs = [_getBinObjForOptsObjGroup(group) for group in optsObjsGroups]

	#4) Loop over trajectory + bin values
	for trajStep in inpTraj:
		currGeom = trajStep.unitCell
		sparseMatrixCalculator.calcMatricesForGeom(currGeom)
		for binValGetter,binObj in it.zip_longest(binValGetters,outBinObjs):
			currBinVals = binValGetter.getValsToBin(sparseMatrixCalculator)
			binObj.addBinValuesToCounts(currBinVals)

	#Get the normalised counts
	nSteps = len(inpTraj.trajSteps)
	for binObj in outBinObjs:
		waterComboDistrHelp._attachCountsNormalisedByNStepsToBin(binObj, nSteps)


	return outBinObjs


def _checkIndicesConsistentInOptsObjGroup(optsObjGroup):
	primaryIndices = [_getPrimaryIndicesFromOptObj(optObj) for optObj in optsObjGroup]
	if any([x!=primaryIndices[0] for x in primaryIndices]):
		raise ValueError("Primary indices not consistent within optsObjGroup")


def _getPrimaryIndicesFromOptObj(optObj):
	try:
		outIndices = optObj.indicesA
	except AttributeError:
		outIndices = optObj.indices

	return outIndices

def _getBinObjForOptsObjGroup(optsObjGroup):
	binEdges = [x.binResObj.binEdges for x in optsObjGroup]
	outBinObj = binResHelp.NDimensionalBinnedResults(binEdges)
	outBinObj.initialiseCountsMatrix()
	return outBinObj


def _getMatrixCalculatorForOptsObjsGroups(optsObjsGroups):
	""" Gets a _SparseMatrixCalculator object when given optsObjsGroups. This handles calculating all the basic matrices needed to calculate any properties we are binning (in as efficient a way as possible, hopefully)
	
	Args:
		optsObjsGroups: (iter of iter of calcOptions objects). Currently these can have "CalcRdfOptions" (with minDistAToB=True) or "" or "CalcPlanarRdfOptions" option objs (a mixture is fine)
			 
	Returns
		matrixCalculator: (_SparseMatrixCalculator)
 
	"""
	#Order doesnt actually matter for distIndicesTo and distIndicesFor
	distIndicesTo, distIndicesFrom = list(), list()

	planeDistEquations, planeDistIndices = list(), list()

	#This loop is sorta annoying...I might be able to chain the opts objects between groups without any real bad effect
	for group in optsObjsGroups:
		for optObj in group:
			if isinstance(optObj, calcDistrCoreHelp.CalcRdfOptions):
				distIndicesFrom.append(optObj.indicesA)
				distIndicesTo.append(optObj.indicesB)
				if optObj.minDistAToB is False:
					raise ValueError("")

			if isinstance(optObj, calcRadImpl.CalcPlanarRdfOptions):
				planeDistIndices.append( optObj.indices )
				currPlaneEqn = _getDefaultPlaneEqn() if optObj.planeEqn is None else optObj.planeEqn
				planeDistEquations.append( currPlaneEqn )


	#Create the matrix calculator
	currKwargs = {"distIndicesFrom":distIndicesFrom, "distIndicesTo":distIndicesTo,
	             "planeDistEquations":planeDistEquations, "planeDistIndices":planeDistIndices}

	#Set any empty lists to None; makes it easier to test
	currKwargs = {k:v for k,v in currKwargs.items() if v!=list()}

	return _SparseMatrixCalculator(**currKwargs)


#TODO: Can add hook for conditional calculators such as H-O-H angles when O-O dists < x
class _SparseMatrixCalculator():
	""" Purpose of this class is to calculate sparsely populated versions of any matrices needed """

	def __init__(self, distIndicesFrom=None, distIndicesTo=None, planeDistEquations=None, planeDistIndices=None):
		""" Initializer. NOT MEANT TO BE CALLED DIRECTLY
		
		Args:
			distIndicesFrom: (iter of iter of ints) Each is a list of indices to calcultate distances FROM. Each matches with iter in distIndicesTo
			distIndicesTo: (iter of iter of ints) Each is a list of indices to calcualte distances TO. Each matches with iter in distIndicesFrom
			planeDistEquations: (iter of ThreeDimPlaneEquation objects) Each should match with planeDistIndices. These are the plane equations to calculate distance from. Code will probably assume ALL are different
			planeDistIndices: (iter of iter of ints) Each is a list of indices for atoms to calculate distance of plane to. Each iter of ints should match with planeDistEquations

		"""
		self.distIndicesFrom = distIndicesFrom
		self.distIndicesTo = distIndicesTo
		self.planeDistEquations = planeDistEquations
		self.planeDistIndices = planeDistIndices

		#For the equality tests
		self._directCmpAttrs = ["distIndicesFrom", "distIndicesTo", "planeDistIndices", "planeDistEquations"]

	def calcMatricesForGeom(self, inpGeom):
		""" Calculates the matrices that we require to get values to bin
		
		Args:
			inpGeom: (plato_pylib UnitCell object)
				 
		Returns
			Nothing; modifies this object in place	
	 
		Notes:
			a) This also calculates values needed (for the most part), but returns sparsely-populated (memory-inefficient) arrays effectively of maximum size (e.g. NAtoms*NAtoms for dist-matrices). These matrices may contain a large number of np.nan values

		Populates:
			distsMatrix: (nxn matrix) N is the number of atoms in the cell. Only actually contains values that are needed; others are np.nan 
			planarDistsMatrix: (nxm matrix) First index is the planeEquation index, second is the atom index

		"""
		self._populatePlaneDistMatrices(inpGeom)
		self._populateDistMatrix(inpGeom)

	#TODO: If it matters (speed wise) factor the part out where i figure out which indices are needed
	def _populatePlaneDistMatrices(self, inpGeom):
		if self.planeDistEquations is None:
			return None

		#Step 1) Figure out the minimal number of things for us to calculate
		uniquePlaneEquations = self.uniquePlaneEquations
		outIndices = [set() for x in uniquePlaneEquations] #
		for currPlaneEqn, currIndices in it.zip_longest(self.planeDistEquations, self.planeDistIndices):
			#Find the relevant plane equation to use
			matched = [x==currPlaneEqn for x in uniquePlaneEquations]
			currIdx = [idx for idx,boolVal in enumerate(matched) if boolVal is True]
			assert len(currIdx)==1
			currIdx = currIdx[0]
			outIndices[currIdx].update( set(currIndices) )

		#Step 2) Actually calculate them
		outMatrices = list()
		for planeEqn, currIndices in it.zip_longest(uniquePlaneEquations, outIndices):
			useIndices = sorted(list(currIndices)) #Shouldnt be neccesary but seems safer in general to pass an ordered list
			currMatrix = calcDistsHelp.calcDistancesFromSurfPlaneForCell(inpGeom, indices=currIndices, planeEqn=planeEqn, sparseMatrix=True)
			outMatrices.append(currMatrix)

		#Step 3) Put them in an accesible place
		self.planarDists = np.array(outMatrices)

	def _populateDistMatrix(self, inpGeom):
		if self.distIndicesFrom is None:
			return None

		#1) Figure out minimal index lists for from and to
		#This will still lead to some pointlessly calculated values but hopefully wont be TOO slow
		fromIndices, toIndices = list(), list()
		for fromIdxs, toIdxs in it.zip_longest(self.distIndicesFrom, self.distIndicesTo):
			fromIndices.extend(fromIdxs)
			toIndices.extend(toIdxs)

		fromIndices = sorted(list(set(fromIndices)))
		toIndices = sorted(list(set(toIndices)))

		#2) Actually calculate the matrix
		currKwargs = {"indicesA":fromIndices, "indicesB":toIndices, "sparseMatrix":True}
		distMatrix = calcDistsHelp.calcDistanceMatrixForCell_minImageConv(inpGeom, **currKwargs)
		
		#3) Put it somewhere accesible
		self.distMatrix = distMatrix


	def __eq__(self,other):

		for attr in self._directCmpAttrs:
			valsA, valsB = getattr(self, attr), getattr(other, attr)
			if valsA != valsB:
				return False

		return True

	@property
	def uniquePlaneEquations(self):
		if self.planeDistEquations is None:
			return None
		outVals = list()
		for planeEqn in self.planeDistEquations:
			if all( [planeEqn!=x for x in outVals] ):
				outVals.append(planeEqn)
		return outVals


class _GetMultiDimValsToBinFromSparseMatrices():
	""" Purpose is in the function getValsToBin() and gets a bunch of values from matrices. One of these should be initialised per iter of calcOptions objects """

	def __init__(self, oneDimObjs):
		self.oneDimObjs = oneDimObjs


	def getValsToBin(self, sparseMatrixCalc):
		""" Takes a bunch of matrices (e.g. distance matrix, horizontal matrix etc.) and returns an iter of values ready to bin
		
		Args:
			sparseMatrixCalc: (_SparseMatrixCalculator object) This stores all matrices that we calculate (e.g. distance matrices)
				 
		Returns
			valsToBin: (iter of len-n iters) 
	 
		"""
		oneDimBinVals = [ x.getValsToBin(sparseMatrixCalc) for x in self.oneDimObjs ]
		multiBinVals = list()

		for binVals in it.zip_longest(*oneDimBinVals):
			multiBinVals.append( binVals )

		return multiBinVals 


	@classmethod
	def fromOptsObjs(cls, optsObjs):
		""" Alternative initializer
		
		Args:
			optsObjs: (iter of options objects, e.g. CalcPlanarRdfOptions)
				 
		"""
		outObjs = list()
		for obj in optsObjs:
			if isinstance(obj, calcDistrCoreHelp.CalcRdfOptions):
				currObj = _MinDistsGetOneDimValsToBin.fromOptsObjs(obj)
			elif isinstance(obj, calcRadImpl.CalcPlanarRdfOptions):
				currObj = _PlanarDistsGetOneDimValsToBin.fromOptsObjs(obj)
			else:
				raise ValueError("")

			outObjs.append(currObj)

		return cls(outObjs)


#This needs subclassing for EACH type of options obj
#So may as well have an overwritable .fromOptionsObj; coupled tightly but its backend code anyway so.....
class _GetOneDimValsToBinFromSparseMatricesBase():
	""" Purpose of this is to get a one-dimensional set of values to bin from the _SparseMatrixCalculator object. One of these should correspond to a single calcOptions object """

	def __init__(self):
		raise NotImplementedError("")

	def getValsToBin(self, sparseMatrixCalculator):
		""" Takes a set of basic matrices (e.g. dist matrices) and returns a 1-dimensional set of values to bin
		
		Args:
			sparseMatrixCalc: (_SparseMatrixCalculator object) This stores all matrices that we calculate (e.g. distance matrices)
				 
		Returns
			valsToBin: (iter of floats) The values to bin
	 
		"""
		raise NotImplementedError("")

	@classmethod
	def fromOptsObj(cls, optsObj):
		raise NotImplementedError("")



class _PlanarDistsGetOneDimValsToBin(_GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, planeEqn, planeDistIndices):
		""" Initializer
		
		Args:
			planeEqn: (ThreeDimPlaneEquation object) The plane equation we want to calculate distances from. Need this to match up to one of the uniquePlaneEquations on the sparseMatrixCalculator we pass to getValsToBin
			planeDistIndices: (iter of ints) The ints of atoms to calculate planar-atom distances
				 
		"""
		self.planeEqn = planeEqn
		self.planeDistIndices = planeDistIndices

	def getValsToBin(self, sparseMatrixCalculator):
		#Get the index of the relevant plane-equation 
		planeEqns = sparseMatrixCalculator.uniquePlaneEquations
		boolVals = [self.planeEqn==x for x in planeEqns]
		assert len([x for x in boolVals if x is True])==1 , "Couldnt find exactly one match for this plane equation"
		planeEqnIdx = boolVals.index(True)

		#
		relPlanarDists = sparseMatrixCalculator.planarDists[planeEqnIdx]
		outVals = [ relPlanarDists[idx] for idx in self.planeDistIndices ]
		return outVals


	@classmethod
	def fromOptsObjs(cls, optsObj):
		""" Alternative initializer
		
		Args:
			optsObj: (CalcPlanarRdfOptions object) Contains the relevant options
				 
		"""
		planeEqn = _getDefaultPlaneEqn() if optsObj.planeEqn is None else optsObj.planeEqn
		planeDistIndices = optsObj.indices
		return cls(planeEqn, planeDistIndices)


class _MinDistsGetOneDimValsToBin(_GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, fromIndices, toIndices):
		""" Initializer
		
		Args:
			fromIndices: (iter of ints) The indices of atoms we calculate distances from. Our output bin values will be len(fromIndices)
			toIndices: (iter of ints) The indices of atoms we calculate distances to
 
		"""
		self.fromIndices = fromIndices
		self.toIndices = toIndices

	def getValsToBin(self, sparseMatrixCalculator):
		relevantMatrix = sparseMatrixCalculator.distMatrix
		outVals = list()
		for idx in self.fromIndices:
			currDists = relevantMatrix[idx][:]
			outVals.append( np.nanmin(currDists[self.toIndices]) )

		return outVals

	@classmethod
	def fromOptsObjs(cls, optsObj):
		""" Alternative initializer
		
		Args:
			optsObj: (CalcRdfOptions object ) Contains the relevant options.
				 
		"""
		if optsObj.minDistAToB is False:
			raise ValueError("")
		fromIndices = optsObj.indicesA
		toIndices = optsObj.indicesB
		return cls(fromIndices, toIndices)


def _getDefaultPlaneEqn():
	return planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0)






