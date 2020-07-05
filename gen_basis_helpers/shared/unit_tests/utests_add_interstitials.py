
import copy
import itertools as it
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
		self.site = "tetrahedral"
		self.ele = "X"
		self.createTestObjs()

	def createTestObjs(self):
		self.atomList = ["X" for x in self.fractPositions]
		kwargDict = {"lattParams":self.lattParams, "lattAngles":self.lattAngles,
		             "fractCoords": self.fractPositions, "elementList":self.atomList}
		self.testCellA = uCell.UnitCell( **kwargDict )

	def testRaisesForIncorrectLatticeAngles(self):
		self.lattAngles[2] = 50
		self.createTestObjs()
		with self.assertRaises(ValueError):
			tCode.addSingleInterToHcpBulkGeom(self.testCellA, self.site)

	def testRaisesForIncorrectSite(self):
		with self.assertRaises(ValueError):
			tCode.addSingleInterToHcpBulkGeom(self.testCellA, "fake_site")

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

	def testTetraBasalInExpectedPositionCoverALessThanPerfect(self):
		self.site = "basal_tetrahedral"
		expOutCell = copy.deepcopy(self.testCellA)
		expCartCoords = expOutCell.cartCoords
		expNewCoord = expCartCoords[0][:2] + [expCartCoords[0][2]+0.5*self.lattParams[-1]] + ["X"]
		expCartCoords.append(expNewCoord)
		expOutCell.cartCoords = expCartCoords
		tCode.addSingleInterToHcpBulkGeom(self.testCellA, self.site, ele="X", strat="utest")
		self.assertEqual(expOutCell,self.testCellA)

	def testOctaInExpectedPositionCoverAlessThanPerfect(self):
		self.site = "octahedral"
		expOutCell = copy.deepcopy(self.testCellA)
		expCartCoords = expOutCell.cartCoords
		expNewCoord = self._getCoordForOctaInter()
		expCartCoords.append(expNewCoord)
		expOutCell.cartCoords = expCartCoords
		tCode.addSingleInterToHcpBulkGeom(self.testCellA, self.site, ele=self.ele, strat="utest")
		self.assertEqual(expOutCell,self.testCellA)

	def testBasalOctaInExpectedPositionCoverAlessThanPerfect(self):
		self.site = "basal_octahedral"
		expOutCell = copy.deepcopy(self.testCellA)
		expCartCoords = expOutCell.cartCoords
		expNewCoord = self._getCoordForOctaInter()
		expNewCoord[2] += 0.25*self.lattParams[2]
		expCartCoords.append(expNewCoord)
		expOutCell.cartCoords = expCartCoords
		tCode.addSingleInterToHcpBulkGeom(self.testCellA, self.site, ele=self.ele, strat="utest")
		self.assertEqual(expOutCell,self.testCellA)

	def testBasalCrowdionInExpectedPositionCoverAlessThanPerfect(self):
		self.site = "basal_crowdion"
		expOutCell = copy.deepcopy(self.testCellA)
		expCartCoords = expOutCell.cartCoords
		expNewCoord = [0.5*self.lattParams[0],0,0] + [self.ele]
		expCartCoords.append(expNewCoord)
		expOutCell.cartCoords = expCartCoords
		tCode.addSingleInterToHcpBulkGeom(self.testCellA, self.site, ele=self.ele, strat="utest")
		self.assertEqual(expOutCell, self.testCellA)

	def testBasalSplitInExpectedPositionCoverALessThanPerfect(self):
		self.site = "basal_split"
		#Need at least TWO atoms in-plane in the regular (non periodic) cell
		self.fractPositions = [[0.0, 0.0, 0.0],
		                       [0.16666666666666677, 0.6666666666666664, 0.4999999999999999],
		                       [0.5, 0.0, 0.0],
		                       [0.6666666666666666, 0.6666666666666664, 0.4999999999999999]]
		self.createTestObjs()

		expOutCell = copy.deepcopy(self.testCellA)
		expCartCoords = expOutCell.cartCoords
		displacement = (2/3)*0.5*0.5*self.lattParams[0]
		newCoordPos = copy.deepcopy(expCartCoords[0])
		newCoordPos[0] += displacement
		expCartCoords[0][0] -= displacement
		expCartCoords.append(newCoordPos)
		expOutCell.cartCoords = expCartCoords
		tCode.addSingleInterToHcpBulkGeom(self.testCellA, self.site, ele=self.ele, strat="utest")
		self.assertEqual(self.testCellA,expOutCell)

	def testCrowdionInExpectedPositionCoverALessThanPerfect(self):
		self.site = "crowdion"
		expOutCell = copy.deepcopy(self.testCellA)
		expCartCoords = expOutCell.cartCoords
		expNewCoord = [0.5*x for x in expCartCoords[1][:3]] + [self.ele]
		expCartCoords.append(expNewCoord)
		expOutCell.cartCoords = expCartCoords
		tCode.addSingleInterToHcpBulkGeom(self.testCellA, self.site, ele=self.ele, strat="utest")
		self.assertEqual(expOutCell, self.testCellA)

	def _getCoordForOctaInter(self):
		#Note the im assuming the start atom is at origin
		a=self.lattParams[0]
		atomAPos = [0,0,0]
		atomBPos = [a,0,0]
		atomCPos = [a*math.cos(math.radians(60)), a*math.sqrt(0.75), 0.0]
		centroidPos = [(v1+v2+v3)/3 for v1,v2,v3 in it.zip_longest(atomAPos,atomBPos,atomCPos)]
		zDisp = 0.25*self.lattParams[2]
		outPos = [x for x in centroidPos]
		outPos[-1] += zDisp
		outPos += [self.ele]
		return outPos



