

import unittest
import gen_basis_helpers.cp2k.method_register as methods
import gen_basis_helpers.cp2k.basis_register as basis
import plato_pylib.shared.ucell_class as UCell
import gen_basis_helpers.cp2k.cp2k_file_helpers as tCode


#TODO: Figure out why this test seems to run relatively slow compared to others
class testModifyCp2kObj(unittest.TestCase):
	""" Testing that we modify the CP2K objects correctly by looking at output strings """

	def setUp(self):
		self.testLattVects = [ [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0] ]
		self.testFractCoords = [ [0.5, 0.5, 0.5, "Mg"] ]
		self.testBasisInfoA = basis.createCP2KBasisObjFromEleAndBasisStr("Mg","utests-basis")
		self.testBasisInfoB = basis.createCP2KBasisObjFromEleAndBasisStr("H","utests-basis")
		self.createTestObjs()

	def createTestObjs(self):
		self.startCP2KObj = methods.createCP2KObjFromMethodStr("cp2k_test_object")
		self.testUCellA = UCell.UnitCell.fromLattVects( self.testLattVects, self.testFractCoords )

	def testDefaultObjGivesExpectedInputStr(self):
		expStr = sorted(_getDefObjInputStr())
		actStr = sorted( self.startCP2KObj.get_input_string() )
		self.assertEqual(expStr, actStr)

	def testModBasedOnDictGivesExpectedObjTestA(self):
		""" Modify multiple options on the CP2K object and check we get expected output file """
		modKwargsDict = {"kpts"       : [20,20,20],
		                 "gridCutAbs" : 300,
		                 "gridCutRel" : 50,
		                 "maxscf"     : 5,
		                 "addedMOs": 10}
	
		expStr = sorted(_loadExpectedOutputStrTestA())
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, modKwargsDict)
		actStr = sorted( self.startCP2KObj.get_input_string() )
		self.assertEqual(expStr, actStr)

	def testAddingSubSystemGivesExpectedObjSingleAtomTestA(self):
		tCode.addSubSysSectionCp2kObj(self.startCP2KObj, self.testLattVects, self.testFractCoords, [self.testBasisInfoA])
		actStr = self.startCP2KObj.get_input_string()
		expStr = sorted(_loadExpectedOutputSubsSystTestSingleAtomA())
		self.assertEqual(expStr, sorted(actStr))

	def testAddUCellGeomAndAllBasisInfoSingleAtomTestA(self):
		tCode.addGeomAndBasisInfoToSimpleCP2KObj(self.startCP2KObj, self.testUCellA, [self.testBasisInfoA])
		actStr = self.startCP2KObj.get_input_string()
		expStr = _loadExpectedOutputAddBasisAndGeomSingleAtomTestA()
		self.assertEqual(sorted(expStr), sorted(actStr))

	def testAddUCellGeomAndAllBasisInfoTwoElementTestA(self):
		""" Check we can use basis sets from multiple files. NOTE: CP2K restricts the use of pseudopotentials to those from only ONE file """
		self.testBasisInfoB.basisFile= "fake_second_basis_file"
		tCode.addGeomAndBasisInfoToSimpleCP2KObj(self.startCP2KObj, self.testUCellA, [self.testBasisInfoA,self.testBasisInfoB])
		actStr = self.startCP2KObj.get_input_string()
		expStr = _loadExpectedOutputAddBasisAndGeomTwoAtomTestA()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testAddUCellGeomAndAllBasisInfoTwoElementsThrowsErrorForMultiplePotFile(self):
		self.testBasisInfoA.potFile = "something"
		self.testBasisInfoB.potFile = "different"
		with self.assertRaises(AssertionError):
			tCode.addGeomAndBasisInfoToSimpleCP2KObj(self.startCP2KObj, self.testUCellA, [self.testBasisInfoA,self.testBasisInfoB])

	def testPrintOrbitalMulliken(self):
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"printAOMullikenPop":True})
		expStr = _loadExpectedOutputStrPrintOrbitalMulliken()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testSetCharge(self):
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"charge":2})
		expStr = _loadExpectedOutputCharge2()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testCellOpt(self):
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"runType":"cell_opt"})
		expStr = _loadExpectedOutputCellOptA()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )



def _getDefObjInputStr():
	defStr = '&GLOBAL\n  PROJECT_NAME cp2k_file\n  PRINT_LEVEL MEDIUM\n  RUN_TYPE ENERGY\n&END GLOBAL\n&FORCE_EVAL\n  METHOD Quickstep\n  &DFT\n    POTENTIAL_FILE_NAME GTH_POTENTIALS\n    BASIS_SET_FILE_NAME BASIS_SET\n    &QS\n      EPS_DEFAULT 1.0E-10\n    &END QS\n    &XC\n      &XC_FUNCTIONAL PBE\n      &END XC_FUNCTIONAL\n    &END XC\n    &KPOINTS\n      SCHEME MONKHORST-PACK 1 1 1\n    &END KPOINTS\n    &MGRID\n      NGRIDS 4\n      REL_CUTOFF [eV] 50000\n      CUTOFF [eV] 5000\n    &END MGRID\n    &SCF\n      SCF_GUESS ATOMIC\n      ADDED_MOS 4\n      EPS_SCF 1.0E-7\n      MAX_SCF 300\n      &MIXING T\n        NBUFFER 8\n        ALPHA 0.4\n        METHOD BROYDEN_MIXING\n      &END MIXING\n      &SMEAR ON\n        ELECTRONIC_TEMPERATURE [K] 157.9\n        METHOD FERMI_DIRAC\n      &END SMEAR\n      &DIAGONALIZATION ON\n        ALGORITHM Standard\n      &END DIAGONALIZATION\n    &END SCF\n  &END DFT\n  &PRINT\n    &FORCES On\n    &END FORCES\n  &END PRINT\n&END FORCE_EVAL\n'
	return defStr

def _loadExpectedOutputStrTestA():
	defStr = '&GLOBAL\n  PROJECT_NAME cp2k_file\n  PRINT_LEVEL MEDIUM\n  RUN_TYPE ENERGY\n&END GLOBAL\n&FORCE_EVAL\n  METHOD Quickstep\n  &DFT\n    POTENTIAL_FILE_NAME GTH_POTENTIALS\n    BASIS_SET_FILE_NAME BASIS_SET\n    &QS\n      EPS_DEFAULT 1.0E-10\n    &END QS\n    &XC\n      &XC_FUNCTIONAL PBE\n      &END XC_FUNCTIONAL\n    &END XC\n    &KPOINTS\n      SCHEME MONKHORST-PACK 20 20 20\n    &END KPOINTS\n    &MGRID\n      NGRIDS 4\n      REL_CUTOFF [eV] 50\n      CUTOFF [eV] 300\n    &END MGRID\n    &SCF\n      SCF_GUESS ATOMIC\n      ADDED_MOS 10\n      EPS_SCF 1.0E-7\n      MAX_SCF 5\n      &MIXING T\n        NBUFFER 8\n        ALPHA 0.4\n        METHOD BROYDEN_MIXING\n      &END MIXING\n      &SMEAR ON\n        ELECTRONIC_TEMPERATURE [K] 157.9\n        METHOD FERMI_DIRAC\n      &END SMEAR\n      &DIAGONALIZATION ON\n        ALGORITHM Standard\n      &END DIAGONALIZATION\n    &END SCF\n  &END DFT\n  &PRINT\n    &FORCES On\n    &END FORCES\n  &END PRINT\n&END FORCE_EVAL\n'
	return defStr

def _loadExpectedOutputSubsSystTestSingleAtomA():
	startStr = _getDefObjInputStr()
	extraStr = '  &SUBSYS\n    &CELL\n      B [bohr] 0.00000000 2.00000000 0.00000000\n      C [bohr] 0.00000000 0.00000000 1.00000000\n      A [bohr] 1.00000000 0.00000000 0.00000000\n    &END CELL\n    &COORD\n      Mg 0.5 0.5 0.5\n      SCALED TRUE\n    &END COORD\n    &KIND Mg\n      POTENTIAL test_potential\n      ELEMENT Mg\n      BASIS_SET test_basis\n    &END KIND\n  &END SUBSYS\n'
	outStr = startStr + extraStr
	return outStr

def _loadExpectedOutputStrPrintOrbitalMulliken():
	outStr = '&GLOBAL\n  PROJECT_NAME cp2k_file\n  PRINT_LEVEL MEDIUM\n  RUN_TYPE ENERGY\n&END GLOBAL\n&FORCE_EVAL\n  METHOD Quickstep\n  &DFT\n    POTENTIAL_FILE_NAME GTH_POTENTIALS\n    BASIS_SET_FILE_NAME BASIS_SET\n    &QS\n      EPS_DEFAULT 1.0E-10\n    &END QS\n    &XC\n      &XC_FUNCTIONAL PBE\n      &END XC_FUNCTIONAL\n    &END XC\n    &KPOINTS\n      SCHEME MONKHORST-PACK 1 1 1\n    &END KPOINTS\n    &MGRID\n      NGRIDS 4\n      REL_CUTOFF [eV] 50000\n      CUTOFF [eV] 5000\n    &END MGRID\n    &SCF\n      SCF_GUESS ATOMIC\n      ADDED_MOS 4\n      EPS_SCF 1.0E-7\n      MAX_SCF 300\n      &MIXING T\n        NBUFFER 8\n        ALPHA 0.4\n        METHOD BROYDEN_MIXING\n      &END MIXING\n      &SMEAR ON\n        ELECTRONIC_TEMPERATURE [K] 157.9\n        METHOD FERMI_DIRAC\n      &END SMEAR\n      &DIAGONALIZATION ON\n        ALGORITHM Standard\n      &END DIAGONALIZATION\n    &END SCF\n    &PRINT\n      &MULLIKEN ON\n        PRINT_GOP TRUE\n      &END MULLIKEN\n    &END PRINT\n  &END DFT\n  &PRINT\n    &FORCES On\n    &END FORCES\n  &END PRINT\n&END FORCE_EVAL\n'
	return outStr


def _loadExpectedOutputAddBasisAndGeomSingleAtomTestA():
	outStr = '&GLOBAL\n  PRINT_LEVEL MEDIUM\n  RUN_TYPE ENERGY\n  PROJECT_NAME cp2k_file\n&END GLOBAL\n&FORCE_EVAL\n  METHOD Quickstep\n  &DFT\n    POTENTIAL_FILE_NAME test_potFile\n    BASIS_SET_FILE_NAME test_basFile\n    &XC\n      &XC_FUNCTIONAL PBE\n      &END XC_FUNCTIONAL\n    &END XC\n    &SCF\n      SCF_GUESS ATOMIC\n      EPS_SCF 1.0E-7\n      MAX_SCF 300\n      ADDED_MOS 4\n      &SMEAR ON\n        METHOD FERMI_DIRAC\n        ELECTRONIC_TEMPERATURE [K] 157.9\n      &END SMEAR\n      &DIAGONALIZATION ON\n        ALGORITHM Standard\n      &END DIAGONALIZATION\n      &MIXING T\n        ALPHA 0.4\n        METHOD BROYDEN_MIXING\n        NBUFFER 8\n      &END MIXING\n    &END SCF\n    &QS\n      EPS_DEFAULT 1.0E-10\n    &END QS\n    &KPOINTS\n      SCHEME MONKHORST-PACK 1 1 1\n    &END KPOINTS\n    &MGRID\n      REL_CUTOFF [eV] 50000\n      CUTOFF [eV] 5000\n      NGRIDS 4\n    &END MGRID\n  &END DFT\n  &PRINT\n    &FORCES On\n    &END FORCES\n  &END PRINT\n  &SUBSYS\n    &CELL\n      B [bohr] 0.00000000 2.00000000 0.00000000\n      A [bohr] 1.00000000 0.00000000 0.00000000\n      C [bohr] 0.00000000 0.00000000 1.00000000\n    &END CELL\n    &COORD\n      Mg 0.49999999999999994 0.5 0.5\n      SCALED TRUE\n    &END COORD\n    &KIND Mg\n      ELEMENT Mg\n      POTENTIAL test_potential\n      BASIS_SET test_basis\n    &END KIND\n  &END SUBSYS\n&END FORCE_EVAL\n'
	return outStr


def _loadExpectedOutputAddBasisAndGeomTwoAtomTestA():
	outStr = "&GLOBAL\n  PROJECT_NAME cp2k_file\n  PRINT_LEVEL MEDIUM\n  RUN_TYPE ENERGY\n&END GLOBAL\n&FORCE_EVAL\n  METHOD Quickstep\n  &DFT\n    POTENTIAL_FILE_NAME test_potFile\n    BASIS_SET_FILE_NAME test_basFile\n    BASIS_SET_FILE_NAME fake_second_basis_file\n    &MGRID\n      CUTOFF [eV] 5000\n      REL_CUTOFF [eV] 50000\n      NGRIDS 4\n    &END MGRID\n    &QS\n      EPS_DEFAULT 1.0E-10\n    &END QS\n    &XC\n      &XC_FUNCTIONAL PBE\n      &END XC_FUNCTIONAL\n    &END XC\n    &KPOINTS\n      SCHEME MONKHORST-PACK 1 1 1\n    &END KPOINTS\n    &SCF\n      ADDED_MOS 4\n      SCF_GUESS ATOMIC\n      MAX_SCF 300\n      EPS_SCF 1.0E-7\n      &DIAGONALIZATION ON\n        ALGORITHM Standard\n      &END DIAGONALIZATION\n      &MIXING T\n        METHOD BROYDEN_MIXING\n        NBUFFER 8\n        ALPHA 0.4\n      &END MIXING\n      &SMEAR ON\n        METHOD FERMI_DIRAC\n        ELECTRONIC_TEMPERATURE [K] 157.9\n      &END SMEAR\n    &END SCF\n  &END DFT\n  &SUBSYS\n    &COORD\n      Mg 0.49999999999999994 0.5 0.5\n      SCALED TRUE\n    &END COORD\n    &CELL\n      B [bohr] 0.00000000 2.00000000 0.00000000\n      C [bohr] 0.00000000 0.00000000 1.00000000\n      A [bohr] 1.00000000 0.00000000 0.00000000\n    &END CELL\n    &KIND Mg\n      POTENTIAL test_potential\n      ELEMENT Mg\n      BASIS_SET test_basis\n    &END KIND\n    &KIND H\n      POTENTIAL test_potential\n      ELEMENT H\n      BASIS_SET test_basis\n    &END KIND\n  &END SUBSYS\n  &PRINT\n    &FORCES On\n    &END FORCES\n  &END PRINT\n&END FORCE_EVAL\n"
	return outStr


def _loadExpectedOutputCharge2():
	defStr = '&GLOBAL\n  PROJECT_NAME cp2k_file\n  PRINT_LEVEL MEDIUM\n  RUN_TYPE ENERGY\n&END GLOBAL\n&FORCE_EVAL\n  METHOD Quickstep\n  &DFT\n    CHARGE 2\n    POTENTIAL_FILE_NAME GTH_POTENTIALS\n    BASIS_SET_FILE_NAME BASIS_SET\n    &QS\n      EPS_DEFAULT 1.0E-10\n    &END QS\n    &XC\n      &XC_FUNCTIONAL PBE\n      &END XC_FUNCTIONAL\n    &END XC\n    &KPOINTS\n      SCHEME MONKHORST-PACK 1 1 1\n    &END KPOINTS\n    &MGRID\n      NGRIDS 4\n      REL_CUTOFF [eV] 50000\n      CUTOFF [eV] 5000\n    &END MGRID\n    &SCF\n      SCF_GUESS ATOMIC\n      ADDED_MOS 4\n      EPS_SCF 1.0E-7\n      MAX_SCF 300\n      &MIXING T\n        NBUFFER 8\n        ALPHA 0.4\n        METHOD BROYDEN_MIXING\n      &END MIXING\n      &SMEAR ON\n        ELECTRONIC_TEMPERATURE [K] 157.9\n        METHOD FERMI_DIRAC\n      &END SMEAR\n      &DIAGONALIZATION ON\n        ALGORITHM Standard\n      &END DIAGONALIZATION\n    &END SCF\n  &END DFT\n  &PRINT\n    &FORCES On\n    &END FORCES\n  &END PRINT\n&END FORCE_EVAL\n'
	return defStr


def _loadExpectedOutputCellOptA():
	outStr = _getDefObjInputStr()
	outStr = outStr.replace("RUN_TYPE ENERGY", "RUN_TYPE CELL_OPT")
	outStr = outStr.replace("&FORCE_EVAL\n", "&FORCE_EVAL\n  STRESS_TENSOR ANALYTICAL\n")
	return outStr
