
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





