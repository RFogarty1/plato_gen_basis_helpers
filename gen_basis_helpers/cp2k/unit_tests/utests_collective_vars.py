
import copy
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.cp2k.method_register as methReg
import gen_basis_helpers.cp2k.collective_vars as tCode
import gen_basis_helpers.shared.plane_equations as planeEqnHelp

class TestMetaVarStandard(unittest.TestCase):

	def setUp(self):
		self.index = 4
		self.scale = 2
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.MetaVarStandard(index=self.index, scale=self.scale)
		self.pyCp2kObj = methReg.createCP2KObjFromMethodStr("spe_standard")

	def testExpectedCaseA(self):
		self.testObjA.addMetaVarToPyCp2kObj(self.pyCp2kObj)
		actMetaVar = self.pyCp2kObj.CP2K_INPUT.MOTION.FREE_ENERGY.METADYN.METAVAR_list[-1]
		actIdx = actMetaVar.Colvar
		actScale = actMetaVar.Scale	
		self.assertEqual(self.index, actIdx)
		self.assertEqual(self.scale, actScale)

#TODO: Alternative constructor which takes 3 atom points for a plane + the unit cell may be good
#Would stop me having to define planePoints explicitly...
class TestDistancePointPlaneFixedPointForPlane_getColvar(unittest.TestCase):

	def setUp(self):
		self.planeIndices = [0,1,2]
		self.atomIdx = 3
		self.cellA = _loadTestCellA()
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.DistancePointPlaneColVar_fixedPointsForPlane.fromInpCellAndAtomIndices(self.cellA, self.atomIdx, self.planeIndices)

	def _runTestFunct(self):
		return self.testObjA.getColvarValueFromInpGeom(self.cellA)

	def testExpectedCaseA_pbcsUnimportant(self):
		expVal = -2 #Maybe other sign...will see
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testExpected_pbcsMatter(self):
		self.cellA.cartCoords = [ [4,4,9,"X"],
	                              [6,6,9,"X"],
	                              [7,6,9,"X"],
	                              [8,8,1,"Y"] ]
		self.createTestObjs()
		expVal = -2
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testExpectedForInpXyz(self):
		inpXyz = [0,0,7]
		expVal = -3 #Plane is at 4, direction seems to point downwards from CP2K point of view
		actVal = self.testObjA.getColVarValueForInpGeomAndXyzPointCoord(self.cellA, inpXyz)
		self.assertAlmostEqual(expVal, actVal)

class TestDistancePointPlaneFixedPointForPlane(unittest.TestCase):

	def setUp(self):
		self.atomPointIndex = 3
		self.planePoints = [ [0,0,1], [2,2,1], [4,2,1] ] #integers make this simpler to deal with
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.DistancePointPlaneColVar_fixedPointsForPlane(self.atomPointIndex, self.planePoints)
		self.pyCp2kObj = methReg.createCP2KObjFromMethodStr("spe_standard")

	def testExpectedCaseA(self):
		self.testObjA.addColVarToSubsys(self.pyCp2kObj)
		actColVar = self.pyCp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS.COLVAR_list[-1].DISTANCE_POINT_PLANE
		
		actColSection = self.pyCp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS.COLVAR_list[-1].DISTANCE_POINT_PLANE
		self.assertEqual( [self.atomPointIndex+1], actColSection.POINT_list[3].Atoms )

		# Deal with the plane
		for idx in range(3):
			self.assertEqual( "FIX_POINT", actColSection.POINT_list[idx].Type )
			self.assertEqual( self.planePoints[idx], actColSection.POINT_list[idx].Xyz )

		#Combine all
		self.assertEqual( 4, actColSection.Atom_point )
		self.assertEqual( [1,2,3], actColSection.Atoms_plane )

	def testConstructorFromInpCellAndAtomIndices(self):
		inpCell = _loadTestCellA()
		#Mod cart coords a bit without changing the average plane equation
		newCartCoords = inpCell.cartCoords
		newCartCoords[0][2] += 1
		newCartCoords[1][2] -= 1
		inpCell.cartCoords = newCartCoords

		#A bit implermentation dependent on the xyz values...
		atomPointIdx, planeIndices = 3, [0,1,2]
		expAtomPointIdx, expPlaneXyz = atomPointIdx, [ [0,0,4], [1,0,4], [0,1,4] ] 
		actObj = tCode.DistancePointPlaneColVar_fixedPointsForPlane.fromInpCellAndAtomIndices(inpCell, atomPointIdx, planeIndices)

		self.assertEqual( expAtomPointIdx, actObj.atomPointIndex )
		for exp,act in it.zip_longest(expPlaneXyz, actObj.planeXyzVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffpointIdx(self):
		objA = copy.deepcopy(self.testObjA)
		self.atomPointIndex += 1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffPlanePoints(self):
		objA = copy.deepcopy(self.testObjA)
		self.planePoints[0][1] += 2
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testToAndFromDictConsistent(self):
		objA = copy.deepcopy(self.testObjA)
		inpDict = self.testObjA.toDict()
		objB = tCode.DistancePointPlaneColVar_fixedPointsForPlane.fromDict(inpDict)
		self.assertEqual(objA, objB)


def _loadTestCellA():
	lattParams, lattAngles = [10,10,10], [90,90,90]
	cartCoords = [ [4,4,4,"X"],
	               [6,6,4,"X"],
	               [7,6,4,"X"],
	               [8,8,6,"Y"] ]
	outCell = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	outCell.cartCoords = cartCoords
	return outCell


class TestDistancePointPlaneColVarAtt2(unittest.TestCase):

	def setUp(self):
		self.atomPlaneIndices = [4,5,6]
		self.atomPointIndex = 2

		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.DistancePointPlaneColVar_att2(self.atomPlaneIndices, self.atomPointIndex)
		self.pyCp2kObj = methReg.createCP2KObjFromMethodStr("spe_standard")

	def testExpectedCaseA(self):
		self.testObjA.addColVarToSubsys(self.pyCp2kObj)
		actColVar = self.pyCp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS.COLVAR_list[-1].DISTANCE_POINT_PLANE
		
		actColSection = self.pyCp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS.COLVAR_list[-1].DISTANCE_POINT_PLANE
		self.assertEqual( [self.atomPointIndex+1], actColSection.POINT_list[3].Atoms )
		self.assertEqual( [self.atomPlaneIndices[0]+1], actColSection.POINT_list[0].Atoms )
		self.assertEqual( [self.atomPlaneIndices[1]+1], actColSection.POINT_list[1].Atoms )
		self.assertEqual( [self.atomPlaneIndices[2]+1], actColSection.POINT_list[2].Atoms )

		self.assertEqual( 4, actColSection.Atom_point )
		self.assertEqual( [1,2,3], actColSection.Atoms_plane )


class TestDistancePointPlaneColVar(unittest.TestCase):

	def setUp(self):
		self.atomPlaneIndices = [4,5,6]
		self.atomPointIndex = 2

		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.DistancePointPlaneColVar(self.atomPlaneIndices, self.atomPointIndex)
		self.pyCp2kObj = methReg.createCP2KObjFromMethodStr("spe_standard")

	def testExpectedCaseA(self):
		self.testObjA.addColVarToSubsys(self.pyCp2kObj)
		actColVar = self.pyCp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS.COLVAR_list[-1].DISTANCE_POINT_PLANE
		expPointIdx = self.atomPointIndex + 1
		expPlaneIndices = [x+1 for x in self.atomPlaneIndices]
		self.assertEqual( expPointIdx, actColVar.Atom_point )
		self.assertEqual( expPlaneIndices, actColVar.Atoms_plane )


class TestGetXyzValsFromSurfPlaneEquation(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [8,9,10], [90,90,90]
		self.d = 3
		self.unitLength = True
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,self.d) #Our normal points along z

	def _runTestFunct(self):
		args = [self.planeEqn, self.cellA]
		kwargs = {"unitLength":self.unitLength}
		return tCode.getXyzValsFromSurfacePlaneEquationAndInpCellStandard(*args, **kwargs)

	def testExpectedValSimpleCaseA(self):
		expVals = [ [0,0,self.d], [1,0,self.d], [0,1,self.d] ]
		actVals = self._runTestFunct()
		self.checkExpAndActValsMatch(expVals, actVals)

	def testExpectedNonUnitLength(self):
		expVals = [ [0,0,self.d], [8,0,self.d], [0,9,self.d] ]
		self.unitLength = False
		actVals = self._runTestFunct()
		self.checkExpAndActValsMatch(expVals, actVals)

	def checkExpAndActValsMatch(self, expVals, actVals):
		self.assertEqual( len(expVals), len(actVals) )

		for expCoord, actCoord in zip(expVals, actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expCoord, actCoord)]





