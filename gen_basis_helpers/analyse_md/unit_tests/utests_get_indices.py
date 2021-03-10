

import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.shared.plane_equations as planeEqnHelp
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
		self.nLayers = 1
		self.coordsA =  [ [0,0,1,"X"], [0,0,1.1], [1,2,3,"Y"], [4,5,6,"Z"], [7,8,9,"Y"] ]
		self.distTol = 1e-2
		self.top = True
		self.bottom = True
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = self._createUCellA()
		self.testFunctA = tCode.FilterToOuterSurfaceAtoms(top=self.top, bottom=self.bottom, distTol=self.distTol, nLayers=self.nLayers)

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

	def testFilterToOuterTwoSurfacePlanes_topOnly(self):
		self.nLayers = 2
		self.bottom = False
		self.createTestObjs()
		inpIndices = [x for x in range(len(self.coordsA))]
		expIndices = [3,4]
		actIndices = self.testFunctA(self.instance, self.cellA, inpIndices)
		self.assertEqual(expIndices, actIndices)

	def testFilterToOuterTwoSurface_botOnly(self):
		self.nLayers = 2
		self.top = False
		self.createTestObjs()
		inpIndices =  [x for x in range(len(self.coordsA))]
		expIndices = [0,1]
		actIndices = self.testFunctA(self.instance, self.cellA, inpIndices)
		self.assertEqual(expIndices, actIndices)


class TestFilterToAtomsWithinDistanceOfSurfacePlane(unittest.TestCase):

	def setUp(self):
		self.top = True
		self.bot = True
		self.maxDist = 3
		self.distTol = 1e-1
		self.planeEqnA = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,5)
		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.coordsA = [ [4,4,4,"X"], [5,5,5,"X"], [6,6,6,"X"] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coordsA
		self.inpIndices = [x for x in range(len(self.coordsA))]
		self.filterFunctA = tCode.FilterToAtomsWithinDistanceOfSurfacePlane(self.planeEqnA, self.maxDist, top=self.top, bottom=self.bot, distTol=self.distTol)

	def _runTestFunct(self):
		return self.filterFunctA(mock.Mock(), self.cellA, self.inpIndices)

	def testExpectedResultsTopAndBottom(self):
		expIndices = [0,1,2]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)

	def testExpectedResultsTopAndBottom_lowerMaxDist(self):
		self.maxDist = 0.5
		self.createTestObjs()
		expIndices = [1]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)

	def testExpectedResultsTopOnly(self):
		self.bot = False
		self.createTestObjs()
		expIndices = [1,2]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testDistTol_returnsBelowWhenTopOnly(self):
		self.planeEqnA = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,5.05)
		self.bot = False
		self.createTestObjs()
		expIndices = [1,2]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)


class TestFilterToExcludeWithoutNebsAmongstRemaining(unittest.TestCase):

	def setUp(self):
		self.maxDist = 1.5
		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.coordsA = [ [5,5,4,"X"], [5,5,5,"X"], [5,5,6,"X"], [5,5,7,"Y"] ]
		self.restrictToPairs = None
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coordsA
		self.inpIndices = [x for x in range(len(self.coordsA))]
		self.filterFunctA = tCode.FilterToExcludeIndicesWithoutNebsAmongstRemaning(self.maxDist, restrictToPairs=self.restrictToPairs)

	def _runTestFunct(self):
		return self.filterFunctA(mock.Mock(), self.cellA, self.inpIndices)

	def testExpected_noneFilteredOut(self):
		expIndices = self.inpIndices
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)

	def testExpected_lastAtomFiltered(self):
		self.inpIndices = [0,1,3]
		expIndices = [0,1]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)

	def testWithRestrictedNebPairs(self):
		self.restrictToPairs = [ ["X","Y"] ]
		self.createTestObjs()
		expIndices = [2,3]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)


class TestFilterToExcludeOutsideOutOfPlaneDistanceFromPoints(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coordsA = [ [5,5,5,"X"],  [5,6,8,"X"], [5,7,8,"X"], [5,8.5,8,"X"] ]
		self.inpPointsA = [ [5,6,2,"X"], [5,7,1,"X"] ]

		self.maxDist = 2
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,5)

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coordsA
		self.inpIndices = [x for x in range(len(self.coordsA))]
		self.filterFunctA = tCode.FilterToExcludeIndicesFurtherOutOfPlaneThanCutoff(self.maxDist, self.planeEqn, self.inpPointsA)

	def _runTestFunct(self):
		return self.filterFunctA(mock.Mock(), self.cellA, self.inpIndices)

	def testExpectedCase_allReturned(self):
		expIndices = [0,1,2,3]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)

	def testExpectedCase_oneMissing(self):
		self.maxDist = 1.1
		self.createTestObjs()
		expIndices = [0,1,2]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)


class TestFilterToExcludeIndicesBasedonNumberOfAtomsInSurfacePlanes(unittest.TestCase):

	def setUp(self):
		#Used to initialise filter function
		self.minAtomsInPlane = 1
		self.maxAtomsInPlane = 1
		self.planeTol = 5e-1
		self.restrictNebsToInpIndices = True

		#geometry used to test
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoordsA = [  [4,4,4,"X"],
		                      [2,2,4,"X"],
		                      [1,1,4,"X"],
		                      [2,2,5,"X"], 
		                      [3,3,6,"X"],
		                      [5,5,6,"X"] ]		
		#Indices to start with
		self.inpIndices = [x for x in range(len(self.cartCoordsA))]

		self.createTestObjs()

	def createTestObjs(self):
		currArgs = [self.minAtomsInPlane, self.maxAtomsInPlane]
		currKwargs = {"planeTol":self.planeTol, "restrictNebsToInpIndices":self.restrictNebsToInpIndices}
		self.testObjA = tCode.FilterToExcludeIndicesBasedOnNumberOfAtomsInSurfacePlane(*currArgs, **currKwargs)

		self.geomA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.geomA.cartCoords = self.cartCoordsA

	def _runTestFunct(self):
		dudInstanceArg = mock.Mock() #filter function takes a getIndicesInstance as first arg; but not needed for this (or most) case
		return self.testObjA(dudInstanceArg, self.geomA, self.inpIndices)

	def testExpectedCase_limitToSingleAtomInPlane(self):
		expIndices = [3]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testExpectedCase_NeighboursStraddlePBCs(self):
		#Only last 2 indices SHOULD be closer than the planeTol (onyl with 2 atoms in a plane)
		self.cartCoordsA = [  [4,4,2,"X"],
		                      [2,2,3,"X"],
		                      [1,1,4,"X"],
		                      [2,2,5,"X"], 
		                      [3,3,0.1,"X"],
		                      [5,5,9.9,"X"] ]
		self.minAtomsInPlane, self.maxAtomsInPlane = 2, 2
		self.createTestObjs()
		expIndices = [4,5]
		actIndices = self._runTestFunct()

		self.assertEqual(expIndices, actIndices)

	def testExpected_nebsRestrictedToInpIndices(self):
		self.inpIndices = [0,3,4,5]
		self.minAtomsInPlane, self.maxAtomsInPlane = 2,5
		self.restrictNebsToInpIndices = True
		self.createTestObjs()
		expIndices = [4,5]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testExpected_nebsNotRestrictedToInpIndices(self):
		self.inpIndices = [0,3,4,5]
		self.minAtomsInPlane, self.maxAtomsInPlane = 2,5
		self.restrictNebsToInpIndices = False
		self.createTestObjs()
		expIndices = [0,4,5]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testThrowsIfMaxAtomsInPlaneLowerThanMin(self):
		self.minAtomsInPlane, self.maxAtomsInPlane = 4, 2
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self._runTestFunct()



