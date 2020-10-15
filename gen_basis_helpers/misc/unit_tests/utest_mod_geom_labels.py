
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.misc.mod_geom_labels as tCode

class TestModSurfaceGeomLabels(unittest.TestCase):

	def setUp(self):
		self.lattParamsA = [3,3,3]
		self.lattAnglesA = [90,90,90]

		self.fractCoordsA = [ [0, 0, 0.5, "Mg"],
		                      [1, 1, 0.5, "Mg"],
		                      [0, 0, 1.0, "Mg"],
		                      [1, 1, 1.0, "Mg"],
		                      [0, 0, 1.5, "Mg"],
		                      [0, 0, 1.5, "Mg"] ]

		self.createTestObjs()

	def createTestObjs(self):
		self.testCellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.testCellA.fractCoords = self.fractCoordsA

	def testExpectedLabelsForCellA(self):
		expLabels = ["Mg_surface", "Mg_surface", "Mg", "Mg", "Mg_surface", "Mg_surface"]
		tCode.modSurfaceAtomLabels(self.testCellA)
		actLabels = [x[-1] for x in self.testCellA.fractCoords]
		self.assertEqual(expLabels, actLabels)

	def testExpectedLabelsWithOneSurfaceAtomExcluded(self):
		expLabels = ["O","Mg_surface", "Mg", "Mg", "Mg_surface", "O"]
		self.fractCoordsA[0][-1] = "O"
		self.fractCoordsA[-1][-1] = "O"
		self.createTestObjs()
		tCode.modSurfaceAtomLabels(self.testCellA, excludeElements=["O"])
		actLabels = [x[-1] for x in self.testCellA.fractCoords]
		self.assertEqual(expLabels, actLabels)

	def testExpectedLablsWithEntireSurfacePlaneExcluded(self):
		expLabels = ["O", "O", "Mg_surface", "Mg_surface", "Mg_surface", "Mg_surface"]
		self.fractCoordsA [0][-1] = "O"
		self.fractCoordsA [1][-1] = "O"
		self.createTestObjs()
		tCode.modSurfaceAtomLabels(self.testCellA, excludeElements=["O"])
		actLabels = [x[-1] for x in self.testCellA.fractCoords]
		self.assertEqual(expLabels, actLabels)

	def testCustomModFunct(self):
		modFunct = lambda x: "fake_ele"
		expLabels = ["fake_ele", "fake_ele", "Mg", "Mg", "fake_ele", "fake_ele"]
		tCode.modSurfaceAtomLabels(self.testCellA, modFunct=modFunct)
		actLabels = [x[-1] for x in self.testCellA.fractCoords]
		self.assertEqual(expLabels, actLabels)



class TestGetKindIdx(unittest.TestCase):

	def setUp(self):
		self.eleListA  = ["Mg","Mg","Mg","Ne","H"]
		self.createTestObjs()

	def createTestObjs(self):
		self.coordsA = [ [0,0,0] + [x] for x in self.eleListA ]	

	def testExpectedGivenA(self):
		testKinds = ["Mg","Ne","H"]
		expIndices = [1,2,3]
		actIndices = [tCode.getKindIdxForLabelFromCoords(self.coordsA,x) for x in testKinds]
		self.assertEqual(expIndices, actIndices)

	def testRaisesExpectedWhenKeyMissing(self):
		testKind = "Not_Present"
		with self.assertRaises(KeyError):
			tCode.getKindIdxForLabelFromCoords(self.coordsA, testKind)

