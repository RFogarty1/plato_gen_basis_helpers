
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






