

import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.get_indices_from_geom_core as tCode

class TestGetSpecialIndicesTemplate(unittest.TestCase):

	def setUp(self):
		self.inpCoordsA = [ [1,2,3,"Mg"], [4,5,6,"Mg"], [5,6,7,"Mg"] ]
		self.retIndicesA = [1,2,3]
		self.retIndicesB = [2]
		self.createTestObjs()

	def createTestObjs(self):
		self.inpGeomA = self._createUCell()
		self.filterFunctA = mock.Mock()
		self.filterFunctB = mock.Mock()
		self.filterFunctA.side_effect = lambda instance, inpGeom, indices: self.retIndicesA 
		self.filterFunctB.side_effect = lambda instance, inpGeom, indices: self.retIndicesB
		currFuncts = [self.filterFunctA, self.filterFunctB]
		self.testObjA = tCode.GetSpecialIndicesFromInpGeomTemplate(currFuncts)

	def _createUCell(self):
		outCell = uCellHelp.UnitCell(lattParams=[8,8,8], lattAngles=[90,90,90])
		outCell.cartCoords = self.inpCoordsA
		return outCell

	def testForSimpleFuncts_noScratchSpace(self):
		expVal = self.retIndicesB
		actVal = self.testObjA.getIndicesFromInpGeom(self.inpGeomA)
		self.filterFunctA.assert_called_with(self.testObjA, self.inpGeomA, [0,1,2])
		self.filterFunctB.assert_called_with(self.testObjA, self.inpGeomA, self.retIndicesA)
		self.assertEqual(expVal,actVal)


class TestFilterToAtomsOfCertainElements(unittest.TestCase):

	def setUp(self):
		self.instance = mock.Mock()
		self.coordsA = [ [0,0,0,"X"], [1,2,3,"Y"], [4,5,6,"Z"], [7,8,9,"Y"] ]
		self.eleTypes = ["X","Y"]
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = self._createUCellA()
		self.testFunctA = tCode.FilterToExcludeElesNotInList(self.eleTypes)

	def _createUCellA(self):
		outCell = uCellHelp.UnitCell(lattParams=[10,10,10], lattAngles=[90,90,90])
		outCell.cartCoords = self.coordsA
		return outCell

	def testExpectedResult_allInputIndices(self):
		inpIndices = [x for x in range(len(self.coordsA))]
		expIndices = [0,1,3]
		actIndices = self.testFunctA(self.instance, self.cellA, inpIndices)
		self.assertEqual(expIndices,actIndices)

	def testExpectedResult_onlySomeIndicesInput(self):
		inpIndices = [0,2,3]
		expIndices = [0,3]
		actIndices = self.testFunctA(self.instance, self.cellA, inpIndices)
		self.assertEqual(expIndices, actIndices)

class TestFilterToOuterPlanes(unittest.TestCase):

	def setUp(self):
		self.instance = mock.Mock()
		self.coordsA =  [ [0,0,1,"X"], [0,0,1.1], [1,2,3,"Y"], [4,5,6,"Z"], [7,8,9,"Y"] ]
		self.distTol = 1e-2
		self.top = True
		self.bottom = True
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = self._createUCellA()
		self.testFunctA = tCode.FilterToOuterSurfaceAtoms(top=self.top, bottom=self.bottom, distTol=self.distTol)

	def _createUCellA(self):
		outCell = uCellHelp.UnitCell(lattParams=[10,10,10], lattAngles=[90,90,90])
		outCell.cartCoords = self.coordsA
		return outCell

	def testExpectedForBothSurfaces(self):
		inpIndices = [x for x in range(len(self.coordsA))]
		expIndices = [0,4]
		actIndices = self.testFunctA(self.instance, self.cellA, inpIndices)
		self.assertEqual(expIndices,sorted(actIndices))

	def testExpectedForBottomSurfaceOnly(self):
		self.top, self.bottom = False, True
		self.createTestObjs()
		inpIndices = [x for x in range(len(self.coordsA))]
		expIndices = [0]
		actIndices = self.testFunctA(self.instance, self.cellA, inpIndices)
		self.assertEqual(expIndices,sorted(actIndices))

	def testExpectedForTopSurfaceOnly(self):
		self.top, self.bottom = True, False
		self.createTestObjs()
		inpIndices = [x for x in range(len(self.coordsA))]
		expIndices = [4]
		actIndices = self.testFunctA(self.instance, self.cellA, inpIndices)
		self.assertEqual(expIndices,sorted(actIndices))

	def testBothSurface_onlySomeIndicesInput(self):
		inpIndices = [0,1,2,3]
		expIndices = [0,3]
		actIndices = self.testFunctA(self.instance, self.cellA, inpIndices)
		self.assertEqual(expIndices,sorted(actIndices))

	def testBottomSurface_higherTolerance(self):
		self.bottom, self.top = True, False
		self.distTol = 5e-1
		self.createTestObjs()
		inpIndices = [x for x in range(len(self.coordsA))]
		expIndices = [0,1]
		actIndices = self.testFunctA(self.instance, self.cellA, inpIndices)
		self.assertEqual(expIndices,sorted(actIndices))

	def testSensibleReturnedForSingleAtomCase(self):
		inpIndices = [1]
		expIndices = [1]
		actIndices = self.testFunctA(self.instance, self.cellA, inpIndices)
		self.assertEqual(expIndices, actIndices)


















