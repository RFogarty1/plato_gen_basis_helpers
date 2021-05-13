
import copy
import unittest
import unittest.mock as mock



import gen_basis_helpers.analyse_md.ads_sites_core as tCode


class StubClass(tCode.FixedIndicesAdsSiteBase):

	def __init__(self, **kwargs):
		for key in kwargs:
			setattr(self, key, kwargs[key])

class TestFixedIndicesAdsSitesEquality(unittest.TestCase):

	def setUp(self):
		self.siteName = "fake_name_a"
		self.atomIndices = [3,4,5]
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"atomIndices":self.atomIndices, "siteName":self.siteName}
		self.testObjA = StubClass(**kwargDict)

	def testEqualObjsCompareEqualA(self):
		objA = copy.deepcopy( self.testObjA )
		self.createTestObjs()
		self.assertEqual(objA, self.testObjA)

	def testEqualObjsCompareEqual_diffOrderAtomIndices(self):
		objA = copy.deepcopy( self.testObjA )
		self.atomIndices = reversed(self.atomIndices)
		self.createTestObjs()
		self.assertEqual(objA, self.testObjA)

	def testUnequalObjsCompareUnequal_diffIndices(self):
		objA = copy.deepcopy(self.testObjA)
		self.atomIndices [1] += 2
		self.createTestObjs()
		self.assertNotEqual(objA, self.testObjA)

	def testUnequalObjsCompareUnequal_diffSiteName(self):
		objA = copy.deepcopy(self.testObjA)
		self.siteName += "_modded"
		self.createTestObjs()
		self.assertNotEqual(objA, self.testObjA)

