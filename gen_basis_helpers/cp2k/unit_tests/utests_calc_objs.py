
import copy

import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.cp2k_calc_objs as tCode
import gen_basis_helpers.cp2k.method_register as methReg


class DudCP2KObj(tCode.CP2KCalcObj):
	""" Point of this is to allow me to test the actual class without interfering with it (e.g. through class decorators) for other test modules """
	pass

class TestCP2KCalcObjGridDescriptors(unittest.TestCase):

	def setUp(self):
		tCode.addRelGridDescriptorToCP2KCalcObjCLASS(DudCP2KObj)
		tCode.addAbsGridDescriptorToCP2KCalcObjCLASS(DudCP2KObj)
		self.pyCP2KObjA = methReg.createCP2KObjFromMethodStr("cp2k_test_object")
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = DudCP2KObj( self.pyCP2KObjA )

	def testAddingRelGridDescriptorGivesConsistentGetterAndSetter(self):
		startVal = copy.deepcopy( self.testObjA.relGridCutoff )
		testVal = 500
		self.assertTrue ( abs(testVal-startVal) > 0.01 ) #Test that our starting value is far from our test value
		self.testObjA.relGridCutoff = testVal
		actVal = self.testObjA.relGridCutoff
		self.assertAlmostEqual(testVal, actVal)

	def testAssertErrorThrownIfRelGridInitUnitsNotInEV(self):
		self.pyCP2KObjA.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.MGRID.Rel_cutoff = "[Ry] 500"
		with self.assertRaises(AssertionError):
			unused = self.testObjA.relGridCutoff

	def testAddingAbsGridDescriptorGivesConsistentGetterAndSetter(self):
		startVal = copy.deepcopy( self.testObjA.absGridCutoff )
		testVal = 500
		self.assertTrue ( abs(testVal-startVal) > 0.01 ) #Test that our starting value is far from our test value
		self.testObjA.absGridCutoff = testVal
		actVal = self.testObjA.absGridCutoff
		self.assertAlmostEqual(testVal, actVal)


class TestCP2KMaxScfDescriptor(unittest.TestCase):

	def setUp(self):
		self.pyCP2KObjA = methReg.createCP2KObjFromMethodStr("cp2k_test_object")
		tCode.addMaxScfDescriptorToCP2KCalcObjCLASS(DudCP2KObj)
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = DudCP2KObj( self.pyCP2KObjA )

	def testGetterAndSetterConsistent(self):
		origMaxCycles = self.testObjA.maxScf
		testMaxCycles = 24
		self.assertNotEqual(origMaxCycles, testMaxCycles)
		self.testObjA.maxScf = testMaxCycles
		actMaxCycles = self.testObjA.maxScf
		self.assertEqual(testMaxCycles, actMaxCycles)


#TODO: These are getting sorta repetitive, can probably use inheritance or composition (supply function to test) to solve it
class TestCP2KAddedMoDescriptor(unittest.TestCase):

	def setUp(self):
		self.pyCP2KObjA = methReg.createCP2KObjFromMethodStr("cp2k_test_object")
		tCode.addAddedMOsDescriptorToCP2KCalcObjCLASS(DudCP2KObj)
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = DudCP2KObj( self.pyCP2KObjA )

	def testGetterAndSetterConsistent(self):
		origAddedMos = self.testObjA.addedMOs
		testAddedMos = 8
		self.assertNotEqual( testAddedMos,origAddedMos )
		self.testObjA.addedMOs = testAddedMos
		self.assertEqual( testAddedMos, self.testObjA.addedMOs )


class TestCP2KCalcObjParsing(unittest.TestCase):

	def setUp(self):
		self.pyCP2KObjA = methReg.createCP2KObjFromMethodStr("cp2k_test_object")
		self.basePath = "fake_base_path_a.inp"
		self.createTestObjs()

	def createTestObjs(self):
		self.testCalcObjA = tCode.CP2KCalcObj(self.pyCP2KObjA, self.basePath)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_calc_objs.os.getcwd")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_calc_objs.os.path.abspath")
	def testExpPathToGeomOutput(self, mockedAbsPath, mockedGetCwd):
		mockedGetCwd.side_effect = lambda : ""
		mockedAbsPath.side_effect = lambda inpPath:inpPath
		self.createTestObjs()
		expPath = "fake_base_path_a-pos-1.xyz"
		actPath = self.testCalcObjA.outGeomPath
		self.assertEqual(expPath,actPath)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_calc_objs.CP2KCalcObj.outGeomPath", new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.cp2k.cp2k_calc_objs.parseCP2K.parseXyzFromGeomOpt")
	def testGetFinalCartCoordsFromOpt(self, mockedXyzParser, mockedPathGetter):
		expFinalGeom, expCartCoords = mock.Mock(), mock.Mock()
		expPath = "fake_path_a"
		expFinalGeom.cartCoords = expCartCoords
		mockedPathGetter.return_value = expPath
		mockedXyzParser.side_effect = lambda *args: {"all_geoms":[mock.Mock(),expFinalGeom]}

		actCartCoords = self.testCalcObjA._getFinalCartCoordsFromOpt()
		mockedXyzParser.assert_called_with(expPath)
		self.assertEqual(expCartCoords,actCartCoords)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_calc_objs.CP2KCalcObj._getFinalCartCoordsFromOpt")
	@mock.patch("gen_basis_helpers.cp2k.cp2k_calc_objs.parseCP2K.parseCpout")
	def testExpectedCartCoordsInUnitCellObj(self, mockedCPOutParser, mockedCoordGetter):
		expUCell, expCartCoords = mock.Mock(), mock.Mock()
		expCpoutDict = {"unitCell":expUCell}
		mockedCPOutParser.side_effect = lambda *args: expCpoutDict
		mockedCoordGetter.side_effect = lambda *args: expCartCoords
		actOutObj = self.testCalcObjA.parsedFile
		actCartCoords = actOutObj.unitCell.cartCoords
		self.assertEqual(expCartCoords,actCartCoords)	

	@mock.patch("gen_basis_helpers.cp2k.cp2k_calc_objs.parseCP2K.parseCpout")
	def testNoErrorRaisedIfOptXyzFileMissing(self, mockedCPOutParser):
		expUCell = mock.Mock()
		mockedCPOutParser.side_effect = lambda *args:{"unitCell":expUCell}
		parsedFile = self.testCalcObjA.parsedFile

