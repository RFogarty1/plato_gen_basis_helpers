
import copy
import math
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCell
import gen_basis_helpers.shared.add_interstitials as tCode

class TestFindNearestNebDistanceForUCell(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,9,8]
		self.lattAngles = [90,90,90]
		self.fractPositions = [ [0.5,0.5,0.5], [0.51,0.5,0.5] ]
		self.atomList = ["X" for x in self.fractPositions]
		self.testAtomIdx = 0
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"lattParams":self.lattParams, "lattAngles":self.lattAngles,
		             "fractCoords": self.fractPositions, "elementList":self.atomList}
		self.testCellA = uCell.UnitCell( **kwargDict )

	def _checkExpMatchesActual(self,expDist):
		actDist = tCode._getNearestNebDistanceForAtomInUcell(self.testCellA, self.testAtomIdx)
		self.assertAlmostEqual(expDist,actDist)

	def testCorrectForOneAtomCell(self):
		self.fractPositions = [ [0.0,0.0,0.0] ]
		self.atomList = ["X" for x in self.fractPositions]
		self.testAtomIdx = 0
		self.createTestObjs()
		expDist = min(self.lattParams)
		self._checkExpMatchesActual(expDist)

	def testCorrectForTwoAtomCellWithoutPeriodicityNEEDED(self):
		""" Test we get the expected value when two non-periodic image atoms are nearest neighbours """
		expDist = 0.01*10
		self._checkExpMatchesActual(expDist)

	def testCorrectForTwoAtomCellWithoutPeriodicityNEEDED_2ndIdx(self):
		expDist = 0.01*10
		self.testAtomIdx = 1
		self._checkExpMatchesActual(expDist)

	def testCorrectForTwoAtomCellWithImageAsNearestNeighbour(self):
		self.fractPositions = [ [0.99, 0.0, 0.0], [0.01, 0.0, 0.0] ]
		self.createTestObjs()
		expDist = 0.02*self.lattParams[0]
		self._checkExpMatchesActual(expDist)


class TestNearestNebDistanceOutOfZPlane(unittest.TestCase):

	def setUp(self):
		self.lattParams = [4,4,4]
		self.lattAngles = [90,90,90]
		self.fractPositions = [ [0.5,0.5,0.5], [0.51,0.5,0.5], [0.5,0.5,0.6] ]
		self.testIdx = 0
		self.zFractTol = 1e-2
		self.createTestObjs()

	def createTestObjs(self):
		atomList = ["X" for x in self.fractPositions]
		kwargDict = {"lattParams":self.lattParams, "lattAngles":self.lattAngles,
		             "fractCoords": self.fractPositions, "elementList":atomList}
		self.testCellA = uCell.UnitCell(**kwargDict)

	def _checkExpMatchesActual(self,expDist):
		actDist = tCode._getNearestNebDistanceOutOfZPlane(self.testCellA, self.testIdx)
		self.assertAlmostEqual(expDist,actDist)

	def testForSimpleCellWhereNearestNebIsInPlane(self):
		expDist = 0.1*self.lattParams[-1]
		self._checkExpMatchesActual(expDist)

class TestAddInterToHcpBulkGeom(unittest.TestCase):

	def setUp(self):
		self.lattParams = [2,2,3]
		self.lattAngles = [90,90,120]
		self.fractPositions = [ [0.0,0.0,0.0], [1/3,2/3,0.5] ]
		self.atomList = ["X" for x in self.fractPositions]
		self.site = "tetrahedral"
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"lattParams":self.lattParams, "lattAngles":self.lattAngles,
		             "fractCoords": self.fractPositions, "elementList":self.atomList}
		self.testCellA = uCell.UnitCell( **kwargDict )

	def testRaisesForIncorrectLatticeAngles(self):
		self.lattAngles[2] = 50
		self.createTestObjs()
		with self.assertRaises(ValueError):
			tCode.addSingleInterToHcpBulkGeom(self.testCellA, self.site)

	def testTetraInExpectedPositionA(self):
		nearestNebOutOfPlaneDist = tCode._getNearestNebDistanceOutOfZPlane(self.testCellA,0)
		nearestNebCartCoords = self.testCellA.cartCoords[1][:3]
		angle = math.degrees(math.acos(0.7924058156930613)) #Calculated "by hand"
		expZDisplace = (0.5*nearestNebOutOfPlaneDist) / math.cos(math.radians(angle)) 
		expOutCell = copy.deepcopy(self.testCellA)
		expCartCoords = expOutCell.cartCoords
		expNewCoord = [ expCartCoords[0][0], expCartCoords[0][1], expCartCoords[0][2] + expZDisplace ]
		expCartCoords.append( expNewCoord + ["X"] )
		expOutCell.cartCoords = expCartCoords
		tCode.addSingleInterToHcpBulkGeom(self.testCellA, self.site, ele="X", strat="utest")
		self.assertEqual(expOutCell,self.testCellA)

	def testTetraInExpectedPositionWithCOverAGreaterThanPerfect(self):
		self.lattParams[-1] = self.lattParams[0]*2
		self.createTestObjs()
		#Largely C+P from above
		nearestNebOutOfPlaneDist = tCode._getNearestNebDistanceOutOfZPlane(self.testCellA, 0)
		angle = 30
		expZDisplace = (0.5*nearestNebOutOfPlaneDist) / math.cos(math.radians(angle)) 
		expOutCell = copy.deepcopy(self.testCellA)
		expCartCoords = expOutCell.cartCoords
		expNewCoord = [ expCartCoords[0][0], expCartCoords[0][1], expCartCoords[0][2] + expZDisplace ]
		expCartCoords.append( expNewCoord + ["X"] )
		expOutCell.cartCoords = expCartCoords
		tCode.addSingleInterToHcpBulkGeom(self.testCellA, self.site, ele="X", strat="utest")
		self.assertEqual(expOutCell,self.testCellA)



