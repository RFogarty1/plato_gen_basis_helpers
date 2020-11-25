
import collections
import copy
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.lammps_interface.lammps_geom as tCode

class TestEleToTypeIdxMapping(unittest.TestCase):

	def setUp(self):

		self.lattParamsA = [10,10,10]
		self.lattAnglesA = [90,90,90]
		self.cartCoordsA = [ [1,2,3,"Mg"],
		                     [4,5,6,"Mg"],
		                     [3,2,1,"H"],
		                     [4,2,1,"Mg"],
		                     [3,3,3,"O"] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.testCellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.testCellA.cartCoords = self.cartCoordsA

	def testExpectedOutputSimpleCaseA(self):
		expDict = {"Mg":1, "H":2, "O":3}
		actDict = tCode.getEleToTypeIdxMapFromUnitCell(self.testCellA)
		self.assertEqual(expDict,actDict)


class TestSimulationBoxObj(unittest.TestCase):

	def setUp(self):
		self.xRangeA = [0,10]
		self.yRangeA = [0,10]
		self.zRangeA = [0,10]
		self.tiltFactors = [0,0,0]
		self.createTestObjs()

	def createTestObjs(self):
		argsA = [self.xRangeA, self.yRangeA, self.zRangeA]
		self.testObjA = tCode.LammpsSimulationBox(*argsA, tiltFactors=self.tiltFactors)

	def testEquality_equalObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA,objB)

	def testEquality_unequalObjsCompareUnequal_rangeVals(self):
		objA = copy.deepcopy(self.testObjA)
		self.yRangeA[-1] += 0.1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)

	def testEquality_unqualObjsCompareUnequal_tiltFactors(self):
		objA = copy.deepcopy(self.testObjA)
		self.tiltFactors[-1] += 0.5
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testExpectedCreatedFromCubicCell(self):
		testCell = uCellHelp.UnitCell(lattParams=[10,10,10], lattAngles=[90,90,90])
		expObj = self.testObjA
		actObj = tCode.LammpsSimulationBox.fromUnitCell(testCell)
		self.assertEqual(expObj, actObj)

	def testExpectedCreatedFromHexagonalCell(self):
		xHi,yHi,zHi, = 10, 8.7, 10
		xy, xz, yz = -5, 0, 0
		expObj = tCode.LammpsSimulationBox( [0,xHi], [0,yHi], [0,zHi], [xy,xz,yz] )
		lattVects = [ [xHi, 0  , 0  ],
		              [xy , yHi, 0  ],
		              [xz , yz , zHi] ]
		testCell = uCellHelp.UnitCell.fromLattVects(lattVects)
		actObj = tCode.LammpsSimulationBox.fromUnitCell(testCell)
		self.assertEqual(expObj,actObj)


class TestGetDataDictAtomTypeFull(unittest.TestCase):

	def setUp(self):
		self.numbFmtCoords = "{:.5f}"
		self.numbFmtCharge = "{:.4f}"
		self.numbFmtMass = "{:.3f}"
		self.lattParamsA = [10,10,10]
		self.lattAnglesA = [90,90,90]
		self.cartCoordsA = [ [ 2,2,2,"O" ],
		                     [ 4,4,4,"H" ],
		                     [ 5,5,5, "H" ] ]
		self.bondInfoA = [ [1, 1, 1, 2],
		                   [2, 1, 1, 3] ]
		self.angleInfoA = [ [1, 1, 2, 1, 3] ]
		self.eleToMassDict = {"O":16, "H":1}
		self.eleToTypeDict = {"O":1,"H":2}
		self.eleToChargeDict = {"O":-2, "H":1}
		self.moleculeIds = [1,2,3]
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"lattParams":self.lattParamsA, "lattAngles":self.lattAnglesA}
		self.testCellA = uCellHelp.UnitCell(**currKwargs)
		self.testCellA.cartCoords = self.cartCoordsA

		self.testGeomObjA = self._getLammpsGeomObj()
#		self.testGeomObjA = tCode.LammpsGeom(self.testCellA)
		self.testFunctObjA = tCode.GetDataDictFromLammpsGeomAtomStyleFull(numbFmtCoords=self.numbFmtCoords,
		                                                              numbFmtCharge=self.numbFmtCharge, numbFmtMasses=self.numbFmtMass)
	#TODO: FINISH THIS
	def _getLammpsGeomObj(self):
		kwargDict = {"eleToTypeIdx":self.eleToTypeDict, "eleToCharge":self.eleToChargeDict, "eleToMass":self.eleToMassDict,
		             "geomToBondInfo":lambda *args:self.bondInfoA, "geomToAngleInfo": lambda *args:self.angleInfoA,
		             "geomToMoleculeIDs":lambda *args:self.moleculeIds}
		return tCode.LammpsGeom(self.testCellA, **kwargDict)

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_geom._MapGeomInfoToDataDict_atomTypeFull")
	def testExpectedCallsMadeForTotalDict(self, mockedGeomInfoToDataDictCls):
		#define a mock for our data map instance
		mockMapInstance = mock.Mock()
		mockedGeomInfoToDataDictCls.side_effect = lambda *args,**kwargs: mockMapInstance

		#Fake strings
		atomsSectionStr, connectSectionStr = "atom_section_str", "conn_section_str"
		headerSectionStr, massesSectionStr = "header_section_str", "masses_section_str"

		mockMapInstance.getAtomsStrFromUnitCell.side_effect = lambda *args,**kwargs: atomsSectionStr
		mockMapInstance.getConnectivityStr.side_effect = lambda *args,**kwargs: connectSectionStr
		mockMapInstance.getHeaderStrFromUnitCellAndKwargs.side_effect = lambda *args,**kwargs: headerSectionStr
		mockMapInstance.getMassesDictStrFromMassDict.side_effect = lambda *args,**kwargs: massesSectionStr

		#Get the expected ordered dict
		expVals = [ ["LAMMPS Atom File", headerSectionStr],
		            ["Masses", massesSectionStr],
		            ["Atoms", atomsSectionStr],
		            ["Bonds", connectSectionStr],
		            ["Angles", connectSectionStr] ]
		expOrderedDict = collections.OrderedDict(expVals)
		expMassDict = {1:16, 2:1}
		actOrderedDict = self.testFunctObjA( self.testGeomObjA )

		#Check the expected calls are made
		currKwargDict = {"numbFmt":self.numbFmtCoords, "bondInfo":self.bondInfoA, "angleInfo":self.angleInfoA}
		mockMapInstance.getHeaderStrFromUnitCellAndKwargs.assert_called_with(self.testCellA,**currKwargDict)
		mockMapInstance.getMassesDictStrFromMassDict.assert_called_with(expMassDict,numbFmt=self.numbFmtMass)
		currKwargDict = {"eleToTypeIdx":self.eleToTypeDict, "eleToChargeIdx": self.eleToChargeDict,
		                 "moleculeIds":self.moleculeIds, "coordNumbFmt":self.numbFmtCoords,
		                 "chargeNumbFmt":self.numbFmtCharge}
		mockMapInstance.getAtomsStrFromUnitCell.assert_called_with(self.testCellA, **currKwargDict)
		mockMapInstance.getConnectivityStr.assert_any_call(self.bondInfoA)
		mockMapInstance.getConnectivityStr.assert_any_call(self.angleInfoA)

		self.assertEqual(expOrderedDict,actOrderedDict)


class TestMapGeomInfoToData_fullAtomStyle_masses(unittest.TestCase):

	def setUp(self):
		self.typeMassDictA = {1:20, 2:14, 3:33}
		self.createTestObjs()

	def createTestObjs(self):
		self.testMapperA = tCode.createMapGeomInfoToDataDictObjFromAtomStyle("full")

	def testMassStrAsExpected(self):
		expStr = ""
		expStr += "\t1\t{:.4f}\n".format(20)
		expStr += "\t2\t{:.4f}\n".format(14)
		expStr += "\t3\t{:.4f}\n".format(33)
		actStr = self.testMapperA.getMassesDictStrFromMassDict(self.typeMassDictA)
		self.assertEqual(expStr,actStr)

class TestMapGeomInfoToData_fullAtomStyle_headerBox(unittest.TestCase):

	def setUp(self):
		self.xRangeA = [0,5]
		self.yRangeA = [0,7]
		self.zRangeA = [0,11]
		self.tiltFactorsA = [2,3,1]
		self.lattParamsA = [10,10,10]
		self.lattAnglesA = [90,90,90]
		self.cartCoordsA = [ [1,2,3,"X"], [4,5,6,"Y"], [6,7,8,"X"] ] #Just need the length of them
		self.bondInfoA = [ [1,1,1,2], [1,1,1,3] ]
		self.angleInfoA = [ [1,1,1,2,3] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.testMapperA = tCode.createMapGeomInfoToDataDictObjFromAtomStyle("full")
		currArgs = [self.xRangeA, self.yRangeA, self.zRangeA]
		self.testSimBoxA = tCode.LammpsSimulationBox(*currArgs,self.tiltFactorsA)
		currKwargs = {"lattParams":self.lattParamsA, "lattAngles":self.lattAnglesA}
		self.testCellA = uCellHelp.UnitCell(**currKwargs)
		self.testCellA.cartCoords= self.cartCoordsA

	def testGetSimBoxHeaderSectionAsExpected(self):
		numbFmt = "{:.4f}"
		expStr = ""
		expStr += "\t{:.4f}\t{:.4f}\txlo xhi\n".format(*self.xRangeA)
		expStr += "\t{:.4f}\t{:.4f}\tylo yhi\n".format(*self.yRangeA)
		expStr += "\t{:.4f}\t{:.4f}\tzlo zhi\n".format(*self.zRangeA)
		expStr += "\t{:.4f}\t{:.4f}\t{:.4f}\txy xz yz".format(*self.tiltFactorsA)
		actStr = self.testMapperA._getHeaderCellDimsPartFromLammpsSimBox(self.testSimBoxA)
		self.assertEqual(expStr, actStr)

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_geom._MapGeomInfoToDataDict_atomTypeFull._getHeaderCellDimsPartFromLammpsSimBox")
	@mock.patch("gen_basis_helpers.lammps_interface.lammps_geom.LammpsSimulationBox")
	def testExpectedCallsWhenGettingHeaderCellDimsFromUnitCell(self, mockedSimBoxCls, mockedGetHeader):
		numbFmt = "{:.3f}"
		expOutput, mockSimBoxInstance = mock.Mock(), mock.Mock()
		fakeCell = mock.Mock()
		mockedSimBoxCls.fromUnitCell.side_effect = lambda *args:mockSimBoxInstance
		mockedGetHeader.side_effect = lambda *args,**kwargs: expOutput

		actOutput = self.testMapperA._getHeaderCellDimsPartFromUnitCell(fakeCell, numbFmt=numbFmt)
		mockedSimBoxCls.fromUnitCell.assert_called_with(fakeCell)
		mockedGetHeader.assert_called_with(mockSimBoxInstance, numbFmt=numbFmt)
		self.assertEqual(expOutput, actOutput)

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_geom._MapGeomInfoToDataDict_atomTypeFull._getHeaderCellDimsPartFromUnitCell")
	def testGetHeaderSectionStr(self, mockGetCellDimsStr):
		mockOutStr, numbFmt = "fake_cell_dims_section_str", "fake_fmt_str"
		mockGetCellDimsStr.side_effect = lambda *args, **kwargs: mockOutStr
		expStr = ""
		expStr += "\t3 atoms\n"
		expStr += "\t2 bonds\n"
		expStr += "\t1 angles\n"
		expStr += "\t0 dihedrals\n"
		expStr += "\t0 impropers\n"
		expStr += "\t2 atom types\n"
		expStr += "\t1 bond types\n"
		expStr += "\t1 angle types"
		expStr += "\n" + mockOutStr
		currKwargs = {"numbFmt":numbFmt, "bondInfo":self.bondInfoA, "angleInfo":self.angleInfoA}
		actStr = self.testMapperA.getHeaderStrFromUnitCellAndKwargs(self.testCellA, **currKwargs)
		mockGetCellDimsStr.assert_called_with(self.testCellA, numbFmt=numbFmt)
		self.assertEqual(expStr, actStr)


class TestMapGeomInfoToData_fullAtomStyle_connectivityData(unittest.TestCase):

	def setUp(self):
		self.bondInfoA = [ [1,1,1,2], [1,1,1,3] ]
		self.angleInfoA = [ [1,1,1,2,3] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode._MapGeomInfoToDataDict_atomTypeFull()

	def testBondInfoStrAsExpected(self):
		expStr = ""
		expStr += "\t1\t1\t1\t2\n"
		expStr += "\t1\t1\t1\t3\n"
		actStr = self.testObjA.getConnectivityStr(self.bondInfoA)
		self.assertEqual(expStr, actStr)

	def testAngleInfoStrAsExpected(self):
		expStr = ""
		expStr += "\t1\t1\t1\t2\t3\n"
		actStr = self.testObjA.getConnectivityStr(self.angleInfoA)
		self.assertEqual(expStr,actStr)


class TestMapGeomInfoToData_fullAtomStyle_atomSection(unittest.TestCase):

	def setUp(self):
		self.lattParamsA, self.lattAnglesA = [10,10,10], [90,90,90]
		self.eleToTypeIdx = {"O":1, "H":2}
		self.eleToCharge = {"O":-1.3, "H": 1.5}
		self.moleculeIDs = [1,2,4]
		self.cartCoordsA = [ [1,2,3,"O"],
		                     [4,5,6,"H"],
		                     [7,8,9,"H"] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode._MapGeomInfoToDataDict_atomTypeFull()
		self.testCellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.testCellA.cartCoords = self.cartCoordsA

	def testExpectedAtomsSection(self):
		#NOTE: Format is atomID,moleculeID,atomType,charge,x,y,z
		numbFmtCharge, numbFmtCoord = "{:.4f}", "{:.5f}"
		expStr  = ""
		expStr += "\t1 1 1 {:.4f} {:.5f} {:.5f} {:.5f}\n".format(-1.3, 1, 2, 3) #First oxygen
		expStr += "\t2 2 2 {:.4f} {:.5f} {:.5f} {:.5f}\n".format(1.5, 4, 5, 6) #First hydrogen
		expStr += "\t3 4 2 {:.4f} {:.5f} {:.5f} {:.5f}\n".format(1.5, 7, 8, 9) #Second hydrogen

		currKwargs = {"eleToTypeIdx":self.eleToTypeIdx, "eleToChargeIdx":self.eleToCharge, "moleculeIds":self.moleculeIDs,
		              "coordNumbFmt":numbFmtCoord, "chargeNumbFmt":numbFmtCharge}
		actStr = self.testObjA.getAtomsStrFromUnitCell(self.testCellA, **currKwargs)
		self.assertEqual(expStr,actStr)

	def testRaisesIfReqKwargsUnset(self):
		#Molecule IDs are also needed for this to work, hence should throw if not present
		currKwargs = {"eleToTypeIdx":self.eleToTypeIdx, "eleToChargeIdx":self.eleToCharge}
		with self.assertRaises(ValueError):
			self.testObjA.getAtomsStrFromUnitCell(self.testCellA, **currKwargs)



