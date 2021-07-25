
import unittest
import unittest.mock as mock

import numpy as np

import gen_basis_helpers.shared.heatmap_plotter as tCode

class TestStandardMapFunction(unittest.TestCase):

	def setUp(self):
		edgesA = [1,2,3]
		edgesB = [3,4,5,6]
		self.createTestObjs()

	def createTestObjs(self):
		#Create the bin edge array
		self.binEdgeArray = np.zeros((2,3,2,2))

		self.binEdgeArray[0][0] = [ [1,2], [3,4] ]
		self.binEdgeArray[0][1] = [ [1,2], [4,5] ]
		self.binEdgeArray[0][2] = [ [1,2], [5,6] ]

		self.binEdgeArray[1][0] = [ [2,3], [3,4] ]
		self.binEdgeArray[1][1] = [ [2,3], [4,5] ]
		self.binEdgeArray[1][2] = [ [2,3], [5,6] ]

		#Create a very simple counts array
		self.countsArray = [ [1,2,3],
		                     [4,5,6] ]

	def _runTestFunct(self):
		return tCode.mapBinEdgesAndDataToXYZStandard(self.binEdgeArray, self.countsArray)

	def testExpectedOutputA(self):
		#figure out the X,Y we expect
		expX = np.array( [ [1,1,1,1],
		                   [2,2,2,2],
		                   [3,3,3,3] ] )

		expY = np.array( [ [3,4,5,6],
		                   [3,4,5,6],
		                   [3,4,5,6] ] )

		expZ = self.countsArray
		actX, actY, actZ = self._runTestFunct()

		self.assertTrue( np.allclose(expX,actX) )
		self.assertTrue( np.allclose(expY,actY) )
		self.assertTrue( np.allclose(expZ,actZ) )

	def testRaisesForThreeDimArray(self):
		self.binEdgeArray = np.zeros((2,3,4,2,2))
		with self.assertRaises(ValueError):
			self._runTestFunct()



class TestFourVertexMapFunction(unittest.TestCase):


	# [1,0,0] and [0,1,0] vectors would produce these i guess
	def setUp(self):
		self.verticesAA = [ [0,0,0], [1,0,0], [0,1,0], [1,1,0] ]
		self.verticesAB = [ [1,0,0], [2,0,0], [1,1,0], [2,1,0] ]
		self.verticesBA = [ [0,1,0], [1,1,0], [0,2,0], [1,2,0] ]
		self.verticesBB = [ [1,1,0], [2,1,0], [1,2,0], [2,2,0] ]

		self.verticesArray = [  [self.verticesAA, self.verticesBA],
		                        [self.verticesAB, self.verticesBB] ]

		self.ignoreDim = 2
		self.dataArray = [ [1,2], [3,4] ]

	def _runTestFunct(self):
		args = [self.verticesArray, self.dataArray]
		return tCode.mapFourVertexCoordsAndDataToXYZStandard(*args, ignoreDim=self.ignoreDim)

	def testExpectedCaseA(self):
		expX = [ [0,0,0], [1,1,1], [2,2,2] ]
		expY = [ [0,1,2], [0,1,2], [0,1,2] ]
		expZ = self.dataArray

		actX, actY, actZ = self._runTestFunct()

		self.assertTrue( np.allclose( np.array(expX), np.array(actX) ) )
		self.assertTrue( np.allclose( np.array(expY), np.array(actY) ) )
		self.assertTrue( np.allclose( np.array(expZ), np.array(actZ) ) )


