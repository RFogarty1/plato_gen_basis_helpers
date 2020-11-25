
import collections
import copy
import unittest
import unittest.mock as mock


import gen_basis_helpers.lammps_interface.misc_objs as miscObjHelp
import gen_basis_helpers.lammps_interface.potential_objs as tCode

class TestPriceTIP3PPotentialObj(unittest.TestCase):

	def setUp(self):
		self.ljOO = [2,3]
		self.ljOH = [3,4]
		self.ljHH = [5,6]
		self.eqmLengthOH = 2.5
		self.forceConstOH = 400
		self.eqmAngleHOH = 20
		self.forceConstHOH = 50
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"ljOO":self.ljOO, "ljOH":self.ljOH, "ljHH":self.ljHH, "eqmLengthOH":self.eqmLengthOH,
		              "forceConstOH":self.forceConstOH, "eqmAngleHOH":self.eqmAngleHOH,
		              "forceConstHOH":self.forceConstHOH}
		self.testObjA = tCode.PriceTIP3PPotential(**currKwargs)

	def testExpectedBondPotentials(self):
		expVals = [ tCode.InternalMoleculeHarmonicPotential(["O","H"], self.eqmLengthOH, self.forceConstOH) ]
		actVals = self.testObjA.bondPots
		self.assertEqual(expVals, actVals)

	def testExpectedAnglePotentials(self):
		expVals = [ tCode.InternalMoleculeHarmonicPotential(["H","O","H"], self.eqmAngleHOH, self.forceConstHOH) ]
		actVals = self.testObjA.anglePots
		self.assertEqual(expVals, actVals)

	def testExpectedLennardJonesPotentials(self):
		expVals = [ tCode.LennardJonesPotential(["O","O"], *self.ljOO),
		            tCode.LennardJonesPotential(["O","H"], *self.ljOH),
		            tCode.LennardJonesPotential(["H","H"], *self.ljHH) ]
		actVals = self.testObjA.lennardJonesPots
		self.assertEqual(expVals, actVals)


class TestInternalMolecularHarmonicPotential(unittest.TestCase):

	def setUp(self):
		self.elementsA = ["Mg","H"]
		self.eqmLengthA = 2.05
		self.forceConstantA = 100.05
		self.createTestObjs()

	def createTestObjs(self):
		args = [self.elementsA, self.eqmLengthA, self.forceConstantA]
		self.testObjA = tCode.InternalMoleculeHarmonicPotential(*args)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA, objB)

	def testUnequalCompareUnequal_diffEqmLength(self):
		objA = copy.deepcopy(self.testObjA)
		self.eqmLengthA += 0.1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)


class TestLennardJonesPotential(unittest.TestCase):

	def setUp(self):
		self.elementsA = ["X","Y"]
		self.epsilonA = 0.45
		self.sigmaA = 1.25
		self.createTestObjs()

	def createTestObjs(self):
		args = [self.elementsA, self.epsilonA, self.sigmaA]
		self.testObjA = tCode.LennardJonesPotential(*args)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffSigma(self):
		objA = copy.deepcopy(self.testObjA)
		self.sigmaA +=0.5
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)


class TestCommandsFromGlobalPotOpts(unittest.TestCase):

	def setUp(self):
		self.pairStyle = "pair style string"
		self.kSpaceStyle =  "kSpaceStyleString"
		self.bondStyle = "bond style string"
		self.angleStyle = "angle style string"
		self.dihedralStyle = None
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"pairStyle":self.pairStyle, "kSpaceStyle":self.kSpaceStyle, "bondStyle":self.bondStyle,
		             "angleStyle":self.angleStyle, "dihedralStyle":self.dihedralStyle}
		self.potOptsA = tCode.GlobalPotOptions(**kwargDict)
		self.testObjA = tCode.LammpsPotential(self.potOptsA, None, None)

	def testExpectedCommsA(self):
		expDict = collections.OrderedDict()
		expDict["pair_style"] = self.pairStyle
		expDict["kspace_style"] = self.kSpaceStyle
		expDict["bond_style"] = self.bondStyle
		expDict["angle_style"] = self.angleStyle
#		expDict["dihedral_style"] = self.dihedralStyle #Not present because its "None"
		actDict = self.testObjA._getCommandsFromGlobalPotOpts()
		self.assertEqual(expDict,actDict)


#TODO: Need to redo my input file parser; parsing the potential simply fails at the moment
class TestMapLJAndHarmonicToCommands(unittest.TestCase):

	def setUp(self):
		self.eleToTypeIdx = {"O":3,"H":4}
		self.bondToTypeIdx = {("O","H"):2}
		self.angleToTypeIdx = {("H","O","H"):3}
		self.numbFmt = "{:.2f}"
		#Parameters defininig a fake TIP3P-like potential
		self.ljOO = [2,3]
		self.ljOH = [3,4]
		self.ljHH = [5,6]
		self.eqmLengthOH = 2.5
		self.forceConstOH = 400
		self.eqmAngleHOH = 20
		self.forceConstHOH = 50
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"eleToTypeIdx":self.eleToTypeIdx, "bondToTypeIdx":self.bondToTypeIdx,
		              "angleToTypeIdx":self.angleToTypeIdx}
		typeMapA = miscObjHelp.TypeMaps(**currKwargs)
		self.testToMapA = self._getPotentialObj()
		self.testMapperA = tCode.MapLJAndHarmonicParamsToCommandDict(typeMapA,numbFmt=self.numbFmt)


	def _getPotentialObj(self):
		currKwargs = {"ljOO":self.ljOO, "ljOH":self.ljOH, "ljHH":self.ljHH, "eqmLengthOH":self.eqmLengthOH,
		              "forceConstOH":self.forceConstOH, "eqmAngleHOH":self.eqmAngleHOH,
		              "forceConstHOH":self.forceConstHOH}
		return tCode.PriceTIP3PPotential(**currKwargs)

	def testExpectedValsObtainedA(self):
		expDict = collections.OrderedDict()
		#pair_coeff info
		expDict["pair_coeff"]  = "3 3 {:.2f} {:.2f}".format(*self.ljOO) #OO LJ
		expDict["pair_coeff"] += "\npair_coeff " + "3 4 {:.2f} {:.2f}".format(*self.ljOH)  #OH LJ
		expDict["pair_coeff"] += "\npair_coeff " + "4 4 {:.2f} {:.2f}".format(*self.ljHH)  #HH LJ

		#bond_coeff and angle_coeff info
		expDict["bond_coeff"] =  "2 {:.2f} {:.2f}".format(self.forceConstOH , self.eqmLengthOH)
		expDict["angle_coeff"] = "3 {:.2f} {:.2f}".format(self.forceConstHOH, self.eqmAngleHOH)

		actDict = self.testMapperA(self.testToMapA)
		self.assertEqual(expDict, actDict)



class TestLammpsPotentialCommandDictCallsExpected(unittest.TestCase):

	def setUp(self):
		self.globalPotOpts = mock.Mock()
		self.potParams = mock.Mock()
		self.mapPotParams = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		args = [self.globalPotOpts, self.potParams, self.mapPotParams]
		self.testObjA = tCode.LammpsPotential(*args)

	@mock.patch("gen_basis_helpers.lammps_interface.potential_objs.LammpsPotential._getCommandsFromGlobalPotOpts")
	def testCommandDictCallsExpected(self, mockedCommsFromGlobalPotOpts):
		keyValA, keyValB = ["keyA","valA"], ["keyB","valB"]
		expDict = collections.OrderedDict( [keyValA,keyValB] )
		mockedCommsFromGlobalPotOpts.side_effect = lambda: collections.OrderedDict([keyValA])
		self.mapPotParams.side_effect = lambda x: collections.OrderedDict([keyValB])
		actDict = self.testObjA.commandDict
		mockedCommsFromGlobalPotOpts.assert_called_with()
		self.mapPotParams.assert_called_with(self.potParams)
		self.assertEqual(expDict,actDict)



