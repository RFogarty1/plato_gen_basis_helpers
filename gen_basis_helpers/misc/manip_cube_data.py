
import itertools as it

import numpy as np

import gen_basis_helpers.shared.simple_vector_maths as vectHelp

def getCubeDataObjFromDictSimple(cubeDict):
	""" Returns a CubeDataSimple instance (useful for various analysis) from cubeDict (standard parsed format)
	
	Args:
		cubeDict: (dict) Dictionary obtained from parsing a cube file with plato_pylib parseCubeFile function
			 
	Returns
		cubeDataSimple: (CubeDataSimple instance) This contains positions and values for the input cube data. Useful for doing various analyses
 
	"""
	vectors = [ cubeDict[key] for key in ["step_x", "step_y", "step_z"] ]
	values = cubeDict["data_grid"]
	origin = cubeDict["origin"]
	atomCoords = [ coord+[atNumber] for coord,atNumber in it.zip_longest(cubeDict["atomic_coords"], cubeDict["atomic_numbers"]) ]  

	outArgs = [vectors, values]
	outKwargs = {"origin":origin,"atomCoords":atomCoords}
	return CubeDataSimple(*outArgs, **outKwargs)


class CubeDataSimple():
	""" Class for storing data from cube files """

	def __init__(self, vectors, values, origin=None, atomCoords=None):
		""" Initializer
		
		Args:
			vectors: [vectA, vectB, vectC]. These are "step_x"/"step_y"/"step_z" and they define the individual cubes
			values: (iter of len-3 float iters)
			origin: (len-3 iter) The origin of the co-ordinate system. Default is [0,0,0]
			atomCoords: (iter of len-4 iters) [x,y,z, atomicNumber] Last entry can probably also be a str (e.g. "Mg")
 
		"""
		self.vectors = vectors
		self.values = values
		self.origin = [0,0,0] if origin is None else origin
		self.atomCoords = atomCoords

		#Deal with equality attrs
		self._arrayAttrs = ["vectors", "values", "origin"]

	@property
	def lenEachDim(self):
		""" Gets us the number of cells in each dimension; returns a list """
		nDims = len(self.vectors)
		outVals = list()
		outVals.append( len(self.values) )

		if nDims>=2:
			outVals.append( len(self.values[0]) )

		if nDims>=3:
			outVals.append( len(self.values[0][0]) )

		return outVals

	@property
	def centres(self):
		""" Gets array of centres for the geometric "cubes" up to three-dimensions
				 
		Returns
			outCentres: (iter of len-n iters) Each element contains the central position of the parallelotope represented by the relevant cell in values. e.g. outCentres[2][1][0] gives the centre for self.values[2][1][0]
	 
		"""
		edges = np.array(self.edges)
		outCentres = np.zeros( tuple( list(edges.shape[:-2]) + [3] ) )

		for comboIdx in np.ndindex( tuple(outCentres.shape[:-1]) ):
			startVertex, endVertex = edges[tuple(comboIdx)][0], edges[tuple(comboIdx)][-1]
			outCentre = [ start + (0.5*(end-start)) for start,end in it.zip_longest(startVertex,endVertex)]
			outCentres[tuple(comboIdx)] = outCentre


		return outCentres


	@property
	def edges(self):
		""" Gets array of edges for the geometric "cubes" (really paralellotopes in the general case) up to three-dimensions
				 
		Returns
			outEdges: (iter of len-n iters) Each element contains edges for the cube represented by the relevant cell in values. e.g. outEdges[2][1][0] gives the edges for self.values[2][1][0]

		Edges:
			a) Each element in outEdges returns an iter; length depends on len(self.vectors) and is 2^{n} (1-d has 2 edges, 2-d has 4 edges, 3-d has 8 edges)
	 		b) Each co-ordinate is a len-3 iter [x,y,z]
			c) The order of elements is left->right (vectA), lower->upper(vectB), bottom->above(vectC).
			For example, for the 3-d case the first four elements are for the bottom of the paralelotope, while the latter four are the top of the paralellotope (they are the first four plus self.vectors[-1]). Of the first four elements, the first two are the origin and origin+vectA. The second two are just vectB+ the value in the first two. 
		"""

		#Initialise the output
		nDims = len(self.vectors)
		edgesPer = 2**(len(self.vectors))

		nIndicesAll = list()

		nIndicesAll.append( len(self.values) ) #1-d; always present 
		if nDims >=2:
			nIndicesAll.append( len(self.values[0]) )
		if nDims>=3:
			nIndicesAll.append( len(self.values[0][0]) )
		
		outEdges = np.zeros( (*nIndicesAll, edgesPer, 3) )

		#combo index will be nDims long, specifying the cell position in terms of the self.vectors
		for comboIdx in np.ndindex( outEdges.shape[:-2] ):
			#Get the current start position
			startPos = [x for x in self.origin]
			for idx, nSteps in enumerate(comboIdx):
				currVector = self.vectors[idx]
				startPos = [x+(t*nSteps) for x,t in it.zip_longest(startPos, currVector)]

			#Setup curr edges
			currEdges = list()

			#Get edges along vectA (ALWAYS runs)
			posLeftLowerBelow = [x for x in startPos]
			posRightLowerBelow = [x+t for x,t in it.zip_longest(startPos, self.vectors[0])]
			currEdges.append( posLeftLowerBelow)
			currEdges.append( posRightLowerBelow)

			#If >=2 dim we get previous translated along vectB
			if nDims>=2:
				posLeftUpperBelow = [x+t for x,t in it.zip_longest(startPos, self.vectors[1])]
				posRightUpperBelow = [x+t for x,t in it.zip_longest(posRightLowerBelow, self.vectors[1])]
				currEdges.append( posLeftUpperBelow )
				currEdges.append( posRightUpperBelow )

			#If 3 dim we get previous translated along vectC
			if nDims>=3:
				posLeftLowerAbove = [x+t for x,t in it.zip_longest(posLeftLowerBelow, self.vectors[2])]
				posLeftUpperAbove = [x+t for x,t in it.zip_longest(posLeftUpperBelow, self.vectors[2])]
				posRightLowerAbove = [x+t for x,t in it.zip_longest(posRightLowerBelow, self.vectors[2])]
				posRightUpperAbove = [x+t for x,t in it.zip_longest(posRightUpperBelow, self.vectors[2])]
			
				currEdges.append( posLeftLowerAbove )
				currEdges.append( posRightLowerAbove )
				currEdges.append( posLeftUpperAbove )
				currEdges.append( posRightUpperAbove )

			outEdges[comboIdx] = currEdges

		return outEdges

#		return outEdges.tolist()


	def __eq__(self, other):
		#Treat most attrs just as np arrays
		for attr in self._arrayAttrs:
			valA, valB = np.array(getattr(self,attr)), np.array(getattr(other,attr))
			if valA.shape != valB.shape:
				return False

			if not np.allclose( np.array(valA), np.array(valB) ):
				return False

		#Treat atom coords a bit differently; the 4th idx in each row may be str or int
		#Hence the pure np-array method isnt really good for it
		atCoordsA, atCoordsB = getattr(self,"atomCoords"), getattr(other,"atomCoords")
		if (atCoordsA is None) and (atCoordsB is None):
			pass
		elif (atCoordsA is None) or (atCoordsB is None):
			return False
		else:
			coordArrA, coordArrB = np.array([x[:3] for x in atCoordsA]), np.array([x[:3] for x in atCoordsB])
			atomValsA, atomValsB = [ x[-1] for x in atCoordsA ], [ x[-1] for x in atCoordsB ] 

			if not np.allclose( coordArrA, coordArrB ):
				return False

			if atomValsA != atomValsB:
				return False

		return True


#Implementation basically stolen from the relevant binned_res.py code
def getLowerDimCubeData_keepSingleIdxForRemovedDims(inpCubeData, keepDims, useIdxOtherDims):
	""" Takes a cubeData instance and returns another with lower dimensionality. 
	
	Args:
		cubeData: (CubeDataSimple)
		keepDims: (iter of ints) The dimensions to keep (e.g. [0,2] may the input for turning 3-D data into 2-D data)
		useIdxOtherDims: (iter of ints) The index of cubes to use for the dimension we throw away (ordered by dimension which we chuck)
 
	Returns
		outData: (CubeDataSimple) Same as input but with lower dimensionality (less vectors). Dimensions are removed by taking a single value for each; e.g. if removing the z-dimension we will take the value at a constant chosen z-value (though note the vectors we remove wont always match up with x/y/z)
 
	NOTE:
		We may or may not be taking view of inpCubeData here but I make zero promises it will stay this way; so its a bad idea to modify any data after doing this (i.e. inpCubeData should be thought as read-only really).

	"""
	#1)Get list of indices we keep
	nDims = len(inpCubeData.vectors)
	keepBools = list()
	for idx in range(nDims):
		currVal = True if idx in keepDims else False
		keepBools.append(currVal)

	#2) Figure out the slice we need
	outSlice = list()
	lenEachDim = inpCubeData.lenEachDim
	idxInUseIdx = 0

	for idx,keep in enumerate(keepBools):
		if keep:
			outSlice.append( slice( lenEachDim[idx] ) )
		else:
			currIdx = useIdxOtherDims[idxInUseIdx]
			outSlice.append( currIdx )
			idxInUseIdx += 1

	#Figure out the output grid vals and vectors
	outGridVals = ( np.array( inpCubeData.values) [tuple(outSlice)] ).tolist()
	outVectors = [vect for vect,boolVal in it.zip_longest(inpCubeData.vectors,keepBools) if boolVal]

	#Create the output object
	outOrigin, outAtomic = inpCubeData.origin, inpCubeData.atomCoords
	outArgs = [outVectors, outGridVals]
	outKwargs = {"origin":outOrigin, "atomCoords":outAtomic}

	return CubeDataSimple(*outArgs, **outKwargs)


def getIdxClosestToInpPointForCubeDataObj(inpPoint, inpDataObj):
	""" Gets the index closest to inpPoint when given a CubeDataSimple instance. This uses the centre of the paralellotopes stored in inpDataObj
	
	Args:
		inpPoint: (len-3 float iter) [x,y,z] for input point
		inpDataObj: (CubeDataSimple)
 
	Returns
		outIdx: (len-n int iter) where n is the number of dimensions
 
	NOTES:
		a) I havent properly tested the equi-distant edge case. It should prioritise high indices by default

	"""
	centres = np.array(inpDataObj.centres)

	#1) Find the distances from centres to each 
	allIndices = [tuple(comboIdx) for comboIdx in np.ndindex(centres.shape[:-1])]
	distMatrix = np.zeros( centres.shape[:-1] ) #centres final index is a len-3 array; whereas we only want a float for the dist

	for comboIdx in allIndices:
		currCentre = centres[comboIdx]
		distMatrix[comboIdx] = vectHelp.getDistTwoVectors(currCentre, inpPoint)

	#Loop over the matrix to find the index of the minimum
	outIdx = tuple([0 for x in range(len(distMatrix.shape))])
	outDist = distMatrix[outIdx]

	for comboIdx in allIndices:
		if distMatrix[comboIdx]<outDist:
			outIdx = comboIdx
			outDist = distMatrix[comboIdx]

	return outIdx


