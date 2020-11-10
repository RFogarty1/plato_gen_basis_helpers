
import copy
import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.cp2k_misc_objs as tCode

class TestGrimmeDispObj(unittest.TestCase):

	def setUp(self):
		self.corrType = "dftd3"
		self.excludeKindsD3 = [1,2]
		self.paramFile = "fake_file.dat"
		self.refFunctional = "blyp"
		self.printDFTD = True
		self.createTestObjs()

	def createTestObjs(self):
		self.kwargDictA = {"corrType": self.corrType, "excludeKindsD3": self.excludeKindsD3,
		                   "paramFile": self.paramFile, "refFunctional": self.refFunctional,
		                   "printDFTD": self.printDFTD}
		self.testObjA = tCode.GrimmeDispersionCorrOptsCP2K(**self.kwargDictA)

	def testEqualObjsCompareEqualA(self):
		objA = copy.deepcopy( self.testObjA )
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA, objB)	

	def testUnequalObjsCompareUnequal_excludeKinds(self):
		objA = copy.deepcopy(self.testObjA)
		self.excludeKindsD3.append(2)
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testToAndFromDictConsistentA(self):
		objA = copy.deepcopy(self.testObjA)
		outDict = objA.toDict()
		objB = tCode.GrimmeDispersionCorrOptsCP2K.fromDict(outDict)
		self.assertEqual(objA,objB)

class testNonLocalDispersionObj(unittest.TestCase):

	def setUp(self):
		self.corrType = "DRSLL"
		self.cutoff = 100
		self.kernelFileName = "fake_kernel"
		self.verboseOutput = True
		self.createTestObjs()

	def createTestObjs(self):
		self.kwargDictA = {"corrType":self.corrType, "cutoff":self.cutoff,
		                   "kernelFileName":self.kernelFileName,
		                   "verboseOutput":self.verboseOutput}
		self.testObjA = tCode.NonLocalDispersionsCorrOptsCP2K(**self.kwargDictA)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA,objB)

	def testUnequalObjsCompareUnequal_cutoffDiff(self):
		objA = copy.deepcopy(self.testObjA)
		self.cutoff += 10
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testToAndFromDictConsistentA(self):
		objA = copy.deepcopy(self.testObjA)
		outDict = objA.toDict()
		objB = tCode.NonLocalDispersionsCorrOptsCP2K.fromDict(outDict)
		self.assertEqual(objA, objB)



class testSurfaceDipoleCorrectionObj(unittest.TestCase):

	def setUp(self):
		self.useCorr = False
		self.surfDipoleDir = "x"
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.SurfaceDipoleCorrectionCP2K(useCorr=self.useCorr, surfDipoleDir=self.surfDipoleDir)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA, objB)

	def testUnequalObjsCompareUnequal(self):
		objA = copy.deepcopy( self.testObjA )
		self.surfDipoleDir = "y"
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual( objA.surfDipoleDir, objB.surfDipoleDir )
		self.assertNotEqual( objA, objB )

	def testToAndFromDictConsistentA(self):
		objA = copy.deepcopy(self.testObjA)
		outDict = objA.toDict()
		objB = tCode.SurfaceDipoleCorrectionCP2K.fromDict(outDict)
		self.assertEqual(objA, objB)





