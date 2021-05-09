
import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.method_register as methReg
import gen_basis_helpers.cp2k.collective_vars as tCode


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



