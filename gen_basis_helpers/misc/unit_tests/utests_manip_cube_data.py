
import copy
import unittest
import unittest.mock as mock

import numpy as np

import gen_basis_helpers.misc.manip_cube_data as tCode



class TestGetCubeDataObjFromCubeDict(unittest.TestCase):

	def setUp(self):
		self._createTestObjs()

	def _createTestObjs(self):
		outDict = dict()
		
		#The xy co-ords
#		outDict["n_x"], outDict["n_y"], outDict["n_z"] = 2, 3, 4
		outDict["origin"] = [0,0,2]
		outDict["step_x"] = [0.8 , 0.0, 0.0]
		outDict["step_y"] = [-0.4, 0.7, 0.0]
		outDict["step_z"] = [0.0 , 0.0, 5.0]

		#Atom coords
		outDict["atomic_coords"] = [ [5,2,1] ]
		outDict["atomic_numbers"] = [ 10 ]

		#vals (1,2,1 shape)
		gridValsA = [ [ [2] ], [ [3] ] ]
		outDict["data_grid"] = gridValsA

		self.cubeDict = outDict		

	def _loadExpectedObjA(self):
		vectors = [ self.cubeDict["step_x"], self.cubeDict["step_y"], self.cubeDict["step_z"] ]  
		values = self.cubeDict["data_grid"]
		origin = self.cubeDict["origin"]
		atomCoords = [  [5,2,1,10] ] 

		args = [vectors, values]
		kwargs = {"origin":origin, "atomCoords":atomCoords}

		return tCode.CubeDataSimple(*args,**kwargs)

	def _runTestFunct(self):
		return tCode.getCubeDataObjFromDictSimple(self.cubeDict)

	def testExpectedCaseA(self):
		expObj = self._loadExpectedObjA()
		actObj = self._runTestFunct()
		self.assertEqual(expObj, actObj)


class TestCubeDataProps(unittest.TestCase):

	def setUp(self):
		self.vectA = [0.8 , 0.0, 0.0]
		self.vectB = [-0.4, 0.7, 0.0]
		self.vectC = [0.0 , 0.0, 5.0]
		self.vectors = [self.vectA, self.vectB, self.vectC]
		self.origin = [0,0,0]
		self.atomCoords = [ [2,2,1,8] ]

		self.vals = np.zeros( (1,2,1) ).tolist()

		self._createTestObjs()

	def _createTestObjs(self):
		currArgs = [self.vectors, self.vals]
		currKwargs = {"origin":self.origin, "atomCoords":self.atomCoords}
		self.testObjA = tCode.CubeDataSimple(*currArgs,**currKwargs)

	def testExpectedEdges_fullMatrix_originAtZero_threeDims(self):

		#left->right (vectA dir), then lower->upper(vectB dir), below->above(vectC dir)
		expEdgesA = [ [0,0,0], [0.8,0,0], [-0.4,0.7,0], [0.4,0.7,0], 
		              [0,0,5], [0.8,0,5], [-0.4,0.7,5], [0.4,0.7,5] ] 

		expEdgesB = [ [-0.4,0.7,0], [0.4,0.7,0], [-0.8,1.4,0], [0,1.4,0],
		              [-0.4,0.7,5], [0.4,0.7,5], [-0.8,1.4,5], [0,1.4,5] ]

		#expEdgesA should be at [0][0][0], expEdgesB should be at [0][1][0]
		expEdgesTot = [[  [expEdgesA] ,  [expEdgesB] ]]
		actEdgesTot = self.testObjA.edges

		self.assertTrue( np.allclose( np.array(expEdgesTot), np.array(actEdgesTot) ) )


	def testExpectedEdges_twoDim_originZero(self):
		self.vectors = [self.vectA, self.vectB]
		self._createTestObjs()

		expEdgesA = [ [0,0,0], [0.8,0,0], [-0.4,0.7,0], [0.4,0.7,0] ]
		expEdgesB = [ [-0.4,0.7,0], [0.4,0.7,0], [-0.8,1.4,0], [0,1.4,0] ]
		
		expEdgesTot = [ [ expEdgesA, expEdgesB] ]
		actEdgesTot = self.testObjA.edges

		self.assertTrue( np.allclose( np.array(expEdgesTot), np.array(actEdgesTot) ) )


	def testExpectedEdges_oneDim_originZero(self):
		self.vectors = [self.vectB]
		self.vals = np.zeros( (2) ).tolist()

		self._createTestObjs()

		expEdgesA = [ [0,0,0], [-0.4,0.7,0] ]
		expEdgesB = [ [-0.4,0.7,0], [-0.8,1.4,0] ]

		expEdgesTot = [ expEdgesA, expEdgesB ]
		actEdgesTot = self.testObjA.edges

		self.assertTrue( np.allclose( np.array(expEdgesTot), np.array(actEdgesTot) ) )


	def testExpectedEdges_oneDim_originNonZero(self):
		self.vectors = [self.vectB]
		self.vals = np.zeros( (2) ).tolist()
		self.origin = [2,3,4]
		self._createTestObjs()

		expEdgesA = [ [2,3,4], [1.6, 3.7, 4] ]
		expEdgesB = [ [1.6,3.7,4], [1.2,4.4,4] ] 

		expEdgesTot = [ expEdgesA, expEdgesB ]
		actEdgesTot = self.testObjA.edges

		self.assertTrue( np.allclose( np.array(expEdgesTot), np.array(actEdgesTot) ) )


	def testExpectedCentres_threeDimCaseA(self):

		#Centre if the midpoint of the diagonal vector
		expCentreA = [ ((b-a)/2)+a for a,b in zip( [ 0  ,0  ,0], [0.4,0.7,5] )]
		expCentreB = [ ((b-a)/2)+a for a,b in zip( [-0.4,0.7,0], [0,1.4,5] ) ]

		#expCentresA should be at [0][0][0], expEdgesB should be at [0][1][0]
		expCentresTot = [[  [expCentreA] ,  [expCentreB] ]]
		actCentresTot = self.testObjA.centres

		self.assertTrue( np.allclose( np.array(expCentresTot), np.array(actCentresTot) ) )

	def testExpectedCentres_twoDimCaseA(self):

		self.vectors = [self.vectA, self.vectB]
		self.vals = np.zeros( (1,2) ).tolist()
		self._createTestObjs()

		#exp Centres
		expCentreA = [ (0.5*(b-a))+a for a,b in zip( [0,0,0], [0.4,0.7,0] ) ] 
		expCentreB = [ (0.5*(b-a))+a for a,b in zip( [-0.4,0.7,0], [0,1.4,0] ) ]

		expCentresTot = [ [expCentreA, expCentreB] ]
		actCentresTot = self.testObjA.centres
		
		self.assertTrue( np.allclose( np.array(expCentresTot), np.array(actCentresTot) ) )


class TestGetCubeIdxClosestToPoint(unittest.TestCase):

	#Note centres (for default) are at:
#[[[0.2, 0.35, 0.0],   (0,0)
# [-0.2, 1.0499999999999998, 0.0]], (0,1)
# [[1.0, 0.35, 0.0], (1,0)
# [0.6000000000000001, 1.0499999999999998, 0.0]], (1,1)
# [[1.8000000000000003, 0.35, 0.0], (2,0)
# [1.4000000000000001, 1.0499999999999998, 0.0]]] (2,1)
	def setUp(self):
		self.vectA = [0.8 , 0.0, 0.0]
		self.vectB = [-0.4, 0.7, 0.0]
		self.vectors = [self.vectA, self.vectB]
		self.inpPoint = [0,0,0]


		self.vals = np.zeros( (3,2,1) ).tolist()

		self._createTestObjs()

	def _createTestObjs(self):
		currArgs = [self.vectors, self.vals]
		self.testObjA = tCode.CubeDataSimple(*currArgs)

	def _runTestFunct(self):
		args = [self.inpPoint, self.testObjA]
		return tCode.getIdxClosestToInpPointForCubeDataObj(*args)

	def testExpectedCaseA_pointAtOrigin(self):
		expIdx = (0,0)
		actIdx = self._runTestFunct()
		self.assertEqual(expIdx, actIdx)

	def testExpectedCaseB(self):
		self.inpPoint = [0.65, 1.0, 0.0]
		expIdx = (1,1)
		actIdx = self._runTestFunct()
		self.assertEqual(expIdx, actIdx)

	def testExpectedCaseC(self):
		self.inpPoint = [1.8, 0.3, 0.0]
		expIdx = (2,0)
		actIdx = self._runTestFunct()
		self.assertEqual( expIdx, actIdx )


class TestGetLowerDimensionCubeData_keepSingleBinIdx(unittest.TestCase):

	def setUp(self):
		self.keepDims = [0,1]
		self.useIdxOtherDims = [2]

		#Can use cubic for simplicity, since shape should have NO EFFECT on how this function works
		self.vectA = [1,0,0] 
		self.vectB = [0,1,0]
		self.vectC = [0,0,1]
		self.vectors = [self.vectA, self.vectB, self.vectC]
		self.vals = self._load2x2x3GridValsA()

		self._createTestObjs()

	def _createTestObjs(self):
		currArgs = [self.vectors, self.vals]
		self.testObj = tCode.CubeDataSimple(*currArgs)

	def _runTestFunct(self):
		currArgs = [self.testObj, self.keepDims, self.useIdxOtherDims]
		return tCode.getLowerDimCubeData_keepSingleIdxForRemovedDims(*currArgs)

	def _load2x2x3GridValsA(self):
		outArray = np.zeros( (2,2,3) )
		outArray[0][0] = [1, 2, 3]
		outArray[0][1] = [6, 3, 4] 
		outArray[1][0] = [5, 7, 9] 
		outArray[1][1] = [4, 2, 1]

		return outArray.tolist()

	def testThreeDimToTwoDimA(self):
		expOutArray = [ [3,4],
		                [9,1] ]
		vectors = [self.vectA, self.vectB]
		expOutObj = tCode.CubeDataSimple(vectors, expOutArray)
		actOutObj = self._runTestFunct()

		self.assertEqual(expOutObj, actOutObj)

	def testThreeDimToTwoDimB(self):
		self.keepDims = [0,2]
		self.useIdxOtherDims = [0]

		expOutArray = [ [1,2,3],
		                [5,7,9] ]
		vectors = [self.vectA, self.vectC]
		expOutObj = tCode.CubeDataSimple(vectors, expOutArray)
		actOutObj = self._runTestFunct()

		self.assertEqual(expOutObj, actOutObj) 

	def testThreeDimToOneDim(self):
		self.keepDims= [2]
		self.useIdxOtherDims = [0,1]

		expOutArray = [6,3,4]
		vectors = [self.vectC]
		expOutObj = tCode.CubeDataSimple(vectors, expOutArray)
		actOutObj = self._runTestFunct()

		self.assertEqual(expOutObj, actOutObj)


class TestCubeDataSimpleEquality(unittest.TestCase):

	def setUp(self):
		self.vectA = [0.8 , 0.0, 0.0]
		self.vectB = [-0.4, 0.7, 0.0]
		self.vectC = [0.0 , 0.0, 5.0]
		self.vectors = [self.vectA, self.vectB, self.vectC]
		self.origin = [0,0,0]
		self.atomCoords = [ [2,2,1,8] ]

		#We use 2x, 2y, 3z as our standard test 
		self.vals = self._load2x2x3GridValsA()
	
		self._createTestObjs()	

	def _load2x2x3GridValsA(self):
		outArray = np.zeros( (2,2,3) )
		outArray[0][0] = [1, 2, 3]
		outArray[0][1] = [6, 3, 4] 
		outArray[1][0] = [5, 7, 9] 
		outArray[1][1] = [4, 2, 1]

		return outArray.tolist()

	def _createTestObjs(self):
		currArgs = [self.vectors, self.vals]
		currKwargs = {"origin":self.origin, "atomCoords":self.atomCoords}
		self.testObjA = tCode.CubeDataSimple(*currArgs,**currKwargs)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self._createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA,objB)
		
	def testUnequalCompareUnequal_diffOrigin(self):
		objA = copy.deepcopy(self.testObjA)
		self.origin[-1] += 1
		self._createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)

	def testUnequalCompareUnequal_diffGridVals(self):
		objA = copy.deepcopy(self.testObjA)
		self.vals[0][1][1] += 2
		self._createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalCompareUnequal_diffNumbVectors(self):
		objA = copy.deepcopy(self.testObjA)
		self.vectors = [self.vectA, self.vectB]
		self._createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalCompareUnequal_diffVectB(self):
		objA = copy.deepcopy(self.testObjA)
		self.vectB[1] += 1
		self._createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalCompareUnequal_atomCoordsNoneAndPresent(self):
		objA = copy.deepcopy(self.testObjA)
		self.atomCoords = None
		self._createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalCompareUnequal_diffCoordVals(self):
		objA = copy.deepcopy(self.testObjA)
		self.atomCoords[0][1] += 1
		self._createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalCompareUnequal_diffCoordEleType(self):
		objA = copy.deepcopy(self.testObjA)
		self.atomCoords[0][-1] += 2
		self._createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalCompareUnequal_diffLenCoords(self):
		objA = copy.deepcopy(self.testObjA)
		self.atomCoords.append( [2,3,4,11] )
		self._createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)






