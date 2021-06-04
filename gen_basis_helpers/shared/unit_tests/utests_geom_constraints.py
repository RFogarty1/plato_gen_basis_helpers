import copy
import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.geom_constraints as tCode

class TestCellConstraintsCls(unittest.TestCase):

	def setUp(self):
		self.angleConstraints = [False,False,False]
		self.lattParamConstraints = [False,False,False]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.CellConstraints(self.angleConstraints, self.lattParamConstraints)

	def testConstraintsPresentIsTrueWhenAnAngleIsFixed(self):
		self.angleConstraints[-1] = True
		self.createTestObjs()
		self.assertTrue( self.testObjA.constraintsPresent is True )

	def testConstraintsPresentIsTrueWhenLattParamIsFixed(self):
		self.lattParamConstraints[-1] = True
		self.createTestObjs()
		self.assertTrue( self.testObjA.constraintsPresent is True )

	def testConstraintsPresentIsFalseWhenNothingIsFixed(self):
		self.angleConstraints = [False for x in range(3)]
		self.lattParamConstraints = [False for x in range(3)]
		self.createTestObjs()
		self.assertTrue( self.testObjA.constraintsPresent is False )

	def testRaisesIfItersNotLenThree(self):
		self.assertTrue( len(self.angleConstraints)==3 )
		self.angleConstraints.append( False )
		with self.assertRaises(AttributeError):
			self.createTestObjs()

	def testFreeInitiationWorks(self):
		testObj = tCode.CellConstraints.initWithNoConstraints()
		self.assertTrue( testObj.constraintsPresent is False )

	def testTwoEqualObjsCompareEqual_noConstraints(self):
		objA = tCode.CellConstraints.initWithNoConstraints()
		objB = tCode.CellConstraints.initWithNoConstraints()
		self.assertEqual(objA,objB)

	def testTwoUnequalObjsCompareUnequal_anglesConstraints(self):
		objA = tCode.CellConstraints.initWithNoConstraints()
		objB = tCode.CellConstraints.initWithNoConstraints()
		objA.anglesToFix = [True,False,True]
		objB.anglesToFix = [True,False,False]
		self.assertNotEqual(objA,objB)

	def testToDictAndFromDictConsistent(self):
		self.angleConstraints = [True,True,True]
		self.createTestObjs()
		objA = copy.deepcopy(self.testObjA)
		dictA = objA.toDict()
		objB = tCode.CellConstraints.fromDict(dictA)
		self.assertEqual(objA, objB)

class TestGeomConstraintsClass(unittest.TestCase):

	def setUp(self):
		self.atomicPositionsConstraints = mock.Mock()
		self.cellConstraints = mock.Mock()
		self.atomicPositionsConstraints = mock.Mock()
		self.atomicPosConstraintsPresent = False
		self.cellConstraintsPresent = False
		self.createTestObjs()

	def createTestObjs(self):
		self.atomicPositionsConstraints.constraintsPresent = self.atomicPosConstraintsPresent
		self.cellConstraints.constraintsPresent = self.cellConstraintsPresent
		self.testObjA = tCode.GeomConstraints(self.atomicPositionsConstraints, self.cellConstraints)
		
	def testConstraintsPresentIsTrueWhenCellConstraintsOn(self):		
		self.cellConstraintsPresent = True
		self.createTestObjs()
		self.assertTrue( self.testObjA.constraintsPresent is True )

	def testInitWithNoConstraints(self):
		testObjA = tCode.GeomConstraints.initWithNoConstraints()
		self.assertTrue( testObjA.constraintsPresent is False )

	def testTwoEqualCompareEqual_noConstraints(self):
		objA = tCode.GeomConstraints.initWithNoConstraints()
		objB = tCode.GeomConstraints.initWithNoConstraints()
		self.assertEqual(objA,objB)

	def testTwoUnequalCompareUnequal_diffCellConstraints(self):
		cellConstrA, cellConstrB = mock.Mock(), mock.Mock()
		objA = tCode.GeomConstraints(self.atomicPositionsConstraints, cellConstrA)
		objB = tCode.GeomConstraints(self.atomicPositionsConstraints, cellConstrB)
		self.assertNotEqual(objA,objB)

	def testTwoUnequalCompareUnequal_diffAtomicPositionConstraints(self):
		atPosConstrA, atPosConstrB = mock.Mock(), mock.Mock()
		objA = tCode.GeomConstraints(atPosConstrA, self.cellConstraints)
		objB = tCode.GeomConstraints(atPosConstrB, self.cellConstraints)
		self.assertNotEqual(objA,objB)


class TestGeomConstraintsToAndFromDict(unittest.TestCase):

	def setUp(self):
		self.cellAnglesToFix = [False,True,False]
		self.atomIndicesToFixX = [3,5]
		self.createTestObjs()

	def createTestObjs(self):
		#cart constraints
		cartConstraints = list()
		for idx in self.atomIndicesToFixX:
			currConstraint = tCode.AtomicCartesianConstraint(idx, fixX=True)
			cartConstraints.append( currConstraint )

		#
		self.testObjA = tCode.GeomConstraints.initWithNoConstraints()
		self.testObjA.cellConstraints.anglesToFix = self.cellAnglesToFix
		self.testObjA.atomicPositionConstraints.atomicCartConstraints = cartConstraints

	def testConsistentCaseA(self):
		expObj = self.testObjA
		inpDict = expObj.toDict()
		actObj = tCode.GeomConstraints.fromDict(inpDict)
		self.assertEqual(expObj, actObj)


class TestAtomicPositionConstraintsClass(unittest.TestCase):

	def setUp(self):
		self.atomicCartConstraints = list()
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.AtomicPositionConstraints( self.atomicCartConstraints )

	def testConstraintsPresentFalseWhenNonePassed(self):
		self.assertFalse(self.testObjA.constraintsPresent)

	def testConstraintsPresentTrueIfOnePassed(self):
		testAtomIdx = 2
		testConstraint = tCode.AtomicCartesianConstraint(testAtomIdx, fixY=True)
		self.atomicCartConstraints.append(testConstraint)
		self.createTestObjs()
		self.assertTrue(self.testObjA.constraintsPresent)

	def testConstraintsPresentTrueIfBlankConstraintPassed(self):
		testAtomIdx = 3
		blankConstraint = tCode.AtomicCartesianConstraint(testAtomIdx)
		self.atomicCartConstraints.append(blankConstraint)
		self.createTestObjs()
		self.assertFalse(self.testObjA.constraintsPresent)

	def testTwoEqualObjsCompareEqual_noConstraints(self):
		objA = tCode.AtomicPositionConstraints.initWithNoConstraints()
		objB = tCode.AtomicPositionConstraints.initWithNoConstraints()
		self.assertEqual(objA,objB)

	def testTwoUnequalObjsCompareUnequal_diffNumberOfConstraints(self):
		atomicConstrA = tCode.AtomicCartesianConstraint(1,fixX=True)
		objA = tCode.AtomicPositionConstraints.initWithNoConstraints()
		objB = tCode.AtomicPositionConstraints(atomicCartConstraints=[atomicConstrA])
		self.assertNotEqual(objA, objB)

	def testTwoCompareEqual_oneWithBlankCartConstraints(self):
		atomicConstrA = tCode.AtomicCartesianConstraint(1) #not really a constraint
		objA = tCode.AtomicPositionConstraints.initWithNoConstraints()
		objB = tCode.AtomicPositionConstraints(atomicCartConstraints=[atomicConstrA])
		self.assertEqual(objA,objB)

	def testTwoCompareUnequal_differentAtomicConstraintObjs(self):
		atomicConstrA = tCode.AtomicCartesianConstraint(1,fixY=True)
		atomicConstrB = tCode.AtomicCartesianConstraint(1,fixZ=True)
		objA = tCode.AtomicPositionConstraints(atomicCartConstraints = [atomicConstrA])
		objB = tCode.AtomicPositionConstraints(atomicCartConstraints = [atomicConstrB])
		self.assertNotEqual(objA,objB)

	def testToAndFromDictConsistent(self):
		#Setup in effect
		testAtomIdx = 2
		testConstraintA = tCode.AtomicCartesianConstraint(testAtomIdx, fixY=True)
		testConstraintB = tCode.AtomicCartesianConstraint(testAtomIdx, fixZ=True)
		self.atomicCartConstraints.append(testConstraintA)
		self.atomicCartConstraints.append(testConstraintB)
		self.createTestObjs()

		#
		expObj = self.testObjA
		inpDict = expObj.toDict()
		actObj = tCode.AtomicPositionConstraints.fromDict(inpDict)
	
class TestAtomicCartConstraints(unittest.TestCase):

	def setUp(self):
		self.atomIdx = 2
		self.fixX = True
		self.fixY = True
		self.fixZ = True
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"fixX":self.fixX, "fixY":self.fixY, "fixZ":self.fixZ}
		self.testObjA = tCode.AtomicCartesianConstraint(self.atomIdx, **kwargDict)

	def testEqualityWorksForEqualObjs(self):
		testObjB = copy.deepcopy(self.testObjA)
		self.assertEqual(testObjB, self.testObjA)

	def testUnequalObjsCompareUnequal_diffConstraints(self):
		objA = copy.deepcopy(self.testObjA)
		self.fixY = not(self.fixY)
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)

	def testUnequalObjsCompareUnequal_diffAtomIdx(self):
		objA = copy.deepcopy(self.testObjA)
		self.atomIdx += 1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testConstraintsPresentReturnsTrueWhenAttrConstrained(self):
		self.fixX, self.fixY = False, False
		self.createTestObjs()
		self.assertTrue(self.testObjA.constraintsPresent)

	def testConstraintsPresentReturnsFalseWhenNoAttrConstrained(self):
		self.fixX, self.fixY, self.fixZ = False, False, False		
		self.createTestObjs()
		self.assertFalse(self.testObjA.constraintsPresent)

	def testToAndFromDictConsistent(self):
		expObj = self.testObjA
		inpDict = expObj.toDict()
		actObj = tCode.AtomicCartesianConstraint.fromDict(inpDict)
		self.assertEqual(expObj,actObj)



