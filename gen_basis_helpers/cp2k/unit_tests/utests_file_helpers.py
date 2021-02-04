

import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.method_register as methods
import gen_basis_helpers.cp2k.basis_register as basis
import plato_pylib.shared.ucell_class as UCell
import gen_basis_helpers.cp2k.cp2k_file_helpers as tCode

import gen_basis_helpers.cp2k.cp2k_misc_objs as miscObjs


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

	def testAddUCellGeomAndAllBasisInfoSingleAtomTestA_ghostAtom(self):
		self.testBasisInfoA.ghost = True
		tCode.addGeomAndBasisInfoToSimpleCP2KObj(self.startCP2KObj, self.testUCellA, [self.testBasisInfoA])
		actStr = self.startCP2KObj.get_input_string()
		expStr = _loadExpectedOutputGhostAtom_a()
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
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"runType":"cell_opt", "geo_constrain_cell_angles":[True,True,True]})
		expStr = _loadExpectedOutputCellOptA()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testSmearingFalseWorks(self):
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"useSmearing":False})
		expStr = _loadExpectedOutputNoSmearing()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testReplaceEpsDef(self):
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"epsDef":5.6e-10})
		expStr = _loadExpectedOutputEpsDef_5pt6_e10()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testReplaceEpsGvgRspace(self):
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"epsGvgRSpace":5.2e-10})
		expStr = _loadExpectedOutputEpsGVG_5p2_e10()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testReplaceEpsPgfOrb(self):
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"epsPgfOrb":4.2e-10})
		expStr = _loadExpectedOutputEpsPgfOrb_4pt2_e10()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testReplaceEpsPPNL(self):
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"epsPPNL":3.2e-10})
		expStr = _loadExpectedOutputEpsPPNL_3pt2_e10()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testEpsRho(self):
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"epsRho":3.7e-10})
		expStr = _loadExpectedOutputEpsRho_3pt7_e10()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testEpsCoreCore(self):
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"epsCoreCharge":4.5e-10})
		expStr = _loadExpectedOutputEpsCoreCharge_4pt5_e10()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testRunTypeBSSE(self):
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"runtype":"bsse", "fragmentsBSSE":[[1,2],[3,4]]})
		expStr = _loadExpectedOutputBSSE_a()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testChangeXCFunctional(self):
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, {"xcFunctional":"blyp"})
		expStr = _loadExpectedOutputChangeXcFunctional()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testChangeGrimmeDispCorr(self):
		kwargDict = {"corrType": "DFTD3", "excludeKindsD3":[1],
		             "paramFile":"fake_file", "refFunctional": "BLYP", "printDFTD":True}
		grimmeObj = miscObjs.GrimmeDispersionCorrOptsCP2K(**kwargDict)
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, grimmeObj.modPyCP2KDict)
		expStr = _loadExpectedOutputChangeGrimmeObjA()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr), sorted(actStr) )

	def testOptB88SpecialFunctional(self):
		kwargDict = {"xcFunctional":"optb88_pw92"}
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, kwargDict)
		expStr = _loadExpectedOutputOptB88Functional()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr.replace(" ","")), sorted(actStr.replace(" ","")) )

	def testChangeDispNLCorr(self):
		kwargDict = {"corrType":"DRSLL", "cutoff":100, "kernelFileName":"fake_file.dat",
		             "verboseOutput":True}
		testOptObj = miscObjs.NonLocalDispersionsCorrOptsCP2K(**kwargDict)
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, testOptObj.modPyCP2KDict)
		expStr = _loadExpectedOutputChangeNonLocalDispCorrObjA()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr.replace(" ","")), sorted(actStr.replace(" ","")) )

	def testChangeScfConvParams(self):
		kwargDict = {"scfMixAlpha": "0.8", "scfMixMethod": "PULAY_MIXING"}
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, kwargDict)
		expStr = _loadExpectedOutputScfMixA()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr.replace(" ","")), sorted(actStr.replace(" ","")) )

	def testSurfDipoleCorr(self):
		kwargDict = {"useCorr":True,"surfDipoleDir":"x"}
		testOptObj = miscObjs.SurfaceDipoleCorrectionCP2K(**kwargDict)
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, testOptObj.modPyCP2KDict)
		expStr = _loadExpectedOutputSurfaceDipoleOptsA()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr.replace(" ","")), sorted(actStr.replace(" ","")) )

	def testMDSection(self):
		kwargDict = {"mdEnsemble":"NVT", "mdSteps":"100",
		             "mdTimeStep":"0.5", "mdTemperature":"300",
		             "mdThermostatType":"nose"}
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, kwargDict)
		expStr = _loadExpectedOutputSimpleMDOptionsA()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr.replace(" ","")), sorted(actStr.replace(" ","")) )

	#Random couple of kwargs to do with MD generally
	def testMiscOptsA(self):
		kwargDict = {"scfPrintRestart":False, "qsExtrapolationMethod":"LINEAR_P", "nGrids":5,
		             "walltime":2500, "prefDiagLib":"sl"}
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, kwargDict)
		expStr = _loadExpectedOutputMiscOptsA()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr.replace(" ","")), sorted(actStr.replace(" ","")) )

	def testTrajPrintOptsA(self):
		kwargDict = {"trajPrintEachMd":50, "trajPrintEachScf":100}
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, kwargDict)
		expStr = _loadExpectedOutputTrajPrintOptsA()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr.replace(" ","")), sorted(actStr.replace(" ","")) )

	def testRestartPrintOptsA(self):
		kwargDict = {"restartPrintEachMd":3}
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, kwargDict)
		expStr = _loadExpectedOutputRestartOptsA()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr.replace(" ","")), sorted(actStr.replace(" ","")) )

	def testAgainstAtomicPosConstraintsA(self):
		kwargDict = {"atPosConstraint_fixIdxPositions":  [2,4,1] , #Note: This is one-indexed in cp2k
		             "atPosConstraint_fixComponents":    ["XY","XY", "XYZ"] } 
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, kwargDict)
		expStr = _loadExpectedOutputAtomicPosConstraintsA()
		actStr = self.startCP2KObj.get_input_string()
		self.assertEqual( sorted(expStr.replace(" ","")), sorted(actStr.replace(" ","")) )

	def testAddColVars(self):
		colVarA, colVarB = mock.Mock(), mock.Mock()
		kwargDict = {"colVars": [colVarA, colVarB]}
		tCode.modCp2kObjBasedOnDict(self.startCP2KObj, kwargDict)
		colVarA.addColVarToSubsys.assert_called_with(self.startCP2KObj)
		colVarB.addColVarToSubsys.assert_called_with(self.startCP2KObj)



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
	constraintStr = "&MOTION\n  &CELL_OPT\n    KEEP_ANGLES TRUE\n  &END CELL_OPT\n&END MOTION\n"
	outStr = constraintStr + outStr
	return outStr


def _loadExpectedOutputNoSmearing():
	outStr = _getDefObjInputStr()
	outStr = outStr.replace("      &SMEAR ON\n        ELECTRONIC_TEMPERATURE [K] 157.9\n        METHOD FERMI_DIRAC\n      &END SMEAR",
	                        "      &SMEAR FALSE\n        ELECTRONIC_TEMPERATURE [K] 157.9\n        METHOD FERMI_DIRAC\n      &END SMEAR")
	return outStr


def _loadExpectedOutputEpsDef_5pt6_e10():
	outStr = _getDefObjInputStr()
	outStr = outStr.replace("EPS_DEFAULT 1.0E-10", "EPS_DEFAULT 5.6E-10")
	return outStr

def _loadExpectedOutputEpsGVG_5p2_e10():
	outStr = _getDefObjInputStr()
	outStr = outStr.replace("EPS_DEFAULT 1.0E-10", "EPS_DEFAULT 1.0E-10\n      EPS_GVG_RSPACE 5.2E-10")
	return outStr

def _loadExpectedOutputEpsPgfOrb_4pt2_e10():
	outStr = _getDefObjInputStr()
	outStr = outStr.replace("EPS_DEFAULT 1.0E-10", "EPS_DEFAULT 1.0E-10\n      EPS_PGF_ORB 4.2E-10")
	return outStr

def _loadExpectedOutputEpsPPNL_3pt2_e10():
	outStr = _getDefObjInputStr()
	outStr = outStr.replace("EPS_DEFAULT 1.0E-10", "EPS_DEFAULT 1.0E-10\n      EPS_PPNL 3.2E-10")
	return outStr

def _loadExpectedOutputEpsRho_3pt7_e10():
	outStr = _getDefObjInputStr()
	outStr = outStr.replace("EPS_DEFAULT 1.0E-10", "EPS_DEFAULT 1.0E-10\n      EPS_RHO 3.7E-10")
	return outStr

def _loadExpectedOutputEpsCoreCharge_4pt5_e10():
	outStr = _getDefObjInputStr()
	outStr = outStr.replace("EPS_DEFAULT 1.0E-10", "EPS_DEFAULT 1.0E-10\n      EPS_CORE_CHARGE 4.5E-10")
	return outStr

def _loadExpectedOutputBSSE_a():
	outStr = _getDefObjInputStr()
	outStr = outStr.replace("RUN_TYPE ENERGY", "RUN_TYPE BSSE")
	outStr = outStr.replace("Quickstep\n", "Quickstep\n  &BSSE\n    &FRAGMENT\n      LIST 1 2\n    &END FRAGMENT\n    &FRAGMENT\n      LIST 3 4\n    &END FRAGMENT\n  &END BSSE\n")
	return outStr

def _loadExpectedOutputGhostAtom_a():
	outStr = _loadExpectedOutputAddBasisAndGeomSingleAtomTestA()
	outStr = outStr.replace("ELEMENT Mg\n      POTENTIAL test_potential\n      BASIS_SET test_basis", "BASIS_SET test_basis\n      GHOST TRUE")
	return outStr

def _loadExpectedOutputChangeXcFunctional():
	outStr = _getDefObjInputStr()
	outStr = outStr.replace("XC_FUNCTIONAL PBE\n","XC_FUNCTIONAL BLYP\n")
	return outStr


def _loadExpectedOutputChangeGrimmeObjA():
	outStr = _getDefObjInputStr()
	newStr = "      &END XC_FUNCTIONAL\n      &VDW_POTENTIAL\n        POTENTIAL_TYPE PAIR_POTENTIAL\n        &PAIR_POTENTIAL\n"
	newStr += "          TYPE DFTD3\n          PARAMETER_FILE_NAME fake_file\n          REFERENCE_FUNCTIONAL BLYP\n"
	newStr += "          D3_EXCLUDE_KIND 1\n          &PRINT_DFTD ON\n          &END PRINT_DFTD\n        &END PAIR_POTENTIAL\n"
	newStr += "      &END VDW_POTENTIAL\n"
	outStr = outStr.replace("      &END XC_FUNCTIONAL\n", newStr) 
	return outStr


def _loadExpectedOutputOptB88Functional():
	outStr = _getDefObjInputStr()
	newStr = "&XC_FUNCTIONAL\n"
	newStr += "        &LIBXC\n"
	newStr += "          FUNCTIONAL XC_GGA_X_OPTB88_VDW\n"
	newStr += "        &END LIBXC\n"
	newStr += "        &PW92 TRUE\n"
	newStr += "        &END PW92\n"
	outStr = outStr.replace("&XC_FUNCTIONAL PBE\n", newStr)
	return outStr


def _loadExpectedOutputChangeNonLocalDispCorrObjA():
	outStr = _getDefObjInputStr()
	newStr = "      &END XC_FUNCTIONAL\n"
	newStr += "      &VDW_POTENTIAL\n"
	newStr += "         POTENTIAL_TYPE NON_LOCAL\n"
	newStr += "         &NON_LOCAL\n"
	newStr += "           TYPE DRSLL\n"
	newStr += "           VERBOSE_OUTPUT TRUE\n"
	newStr += "           KERNEL_FILE_NAME fake_file.dat\n"
	newStr += "           CUTOFF [eV] 100\n"
	newStr += "         &END NON_LOCAL\n"
	newStr += "      &END VDW_POTENTIAL\n"
	outStr = outStr.replace("      &END XC_FUNCTIONAL\n", newStr) 
	return outStr

def _loadExpectedOutputScfMixA():
	outStr = _getDefObjInputStr()
	outStr = outStr.replace("METHOD BROYDEN_MIXING", "METHOD PULAY_MIXING")
	outStr = outStr.replace("ALPHA 0.4", "ALPHA 0.8")
	return outStr

def _loadExpectedOutputSurfaceDipoleOptsA():
	outStr = _getDefObjInputStr()
	newStr = "&DFT\n    SURFACE_DIPOLE_CORRECTION TRUE\n"
	newStr += "    SURF_DIP_DIR X\n"
	outStr = outStr.replace("&DFT\n",newStr)
	return outStr

def _loadExpectedOutputSimpleMDOptionsA():
	outStr = _getDefObjInputStr()
	newStr  = "&MOTION\n &MD\n"
	newStr += "    ENSEMBLE NVT\n"
	newStr += "    STEPS 100\n"
	newStr += "    TIMESTEP 0.5\n"
	newStr += "    TEMPERATURE 300\n"
	newStr += "    &THERMOSTAT\n"
	newStr += "      TYPE NOSE\n"
	newStr += "    &END THERMOSTAT\n"
	newStr += "  &END MD\n"
	newStr += "&END MOTION\n"
	return newStr + outStr 

def _loadExpectedOutputMiscOptsA():
	outStr = _getDefObjInputStr()
	newScfPart = "    &SCF\n      &PRINT\n        &RESTART OFF\n        &END RESTART\n      &END PRINT\n"
	newQsPart = "    &QS\n      EXTRAPOLATION LINEAR_P\n"
	newGlobalPart = "&GLOBAL\n  WALLTIME 2500\n  PREFERRED_DIAG_LIBRARY SL\n"
	outStr = outStr.replace("    &QS\n" ,  newQsPart)
	outStr = outStr.replace("    &SCF\n", newScfPart)
	outStr = outStr.replace("&GLOBAL\n", newGlobalPart)
	outStr = outStr.replace("NGRIDS 4","NGRIDS 5")
	return outStr

def _loadExpectedOutputTrajPrintOptsA():
	outStr = _getDefObjInputStr()
	newStr = "&MOTION\n"
	newStr += "  &PRINT\n"
	newStr += "    &TRAJECTORY\n"
	newStr += "      &EACH\n"
	newStr += "        MD 50\n"
	newStr += "        QS_SCF 100\n"
	newStr += "      &END EACH\n"
	newStr += "    &END TRAJECTORY\n"
	newStr += "  &END PRINT\n"
	newStr += "&END MOTION\n"
	return newStr + outStr

def _loadExpectedOutputRestartOptsA():
	outStr = _getDefObjInputStr()
	newStr = "&MOTION\n"
	newStr += "  &PRINT\n"
	newStr += "    &RESTART\n"
	newStr += "      &EACH\n"
	newStr += "        MD 3\n"
	newStr += "      &END EACH\n"
	newStr += "    &END RESTART\n"
	newStr += "  &END PRINT\n"
	newStr += "&END MOTION\n"
	return newStr + outStr

def _loadExpectedOutputAtomicPosConstraintsA():
	outStr = _getDefObjInputStr()
	newStr =  "&MOTION\n"
	newStr += "  &CONSTRAINT\n"
	newStr += "    &FIXED_ATOMS\n"
	newStr += "      LIST 2\n"
	newStr += "      COMPONENTS_TO_FIX XY\n"
	newStr += "    &END FIXED_ATOMS\n"
	newStr += "    &FIXED_ATOMS\n"
	newStr += "      LIST 4\n"
	newStr += "      COMPONENTS_TO_FIX XY\n"
	newStr += "    &END FIXED_ATOMS\n"
	newStr += "    &FIXED_ATOMS\n"
	newStr += "      LIST 1\n"
	newStr += "      COMPONENTS_TO_FIX XYZ\n"
	newStr += "    &END FIXED_ATOMS\n"
	newStr += "  &END CONSTRAINT\n"
	newStr += "&END MOTION\n"
	return newStr + outStr







