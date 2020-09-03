
import unittest
import unittest.mock as mock

import gen_basis_helpers.fit_cp2k_basis.cp2k_tables_updater as tCode

class TestBasisCoeffUpdaterStandard(unittest.TestCase):

	def setUp(self):
		self.filePathA = "fake_path"
		self.eleName = "Mg"
		self.basisNames = ["basisA"]
		self.coeffToBasisMapperA = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.BasisCoeffUpdaterCP2KStandard(self.filePathA, self.eleName, self.basisNames, self.coeffToBasisMapperA)

	@mock.patch("gen_basis_helpers.fit_cp2k_basis.cp2k_tables_updater.BasisCoeffUpdaterCP2KStandard._updateUsingCurrentCoeffs")
	def testCoeffsChangedOnUpdate(self, mockedUpdateFromSelf):
		testCoeffs = [1,2,3]
		self.testObjA.coeffs = [1]
		self.assertNotEqual(testCoeffs,self.testObjA.coeffs)
		self.testObjA.updateCoeffs(testCoeffs)
		self.assertEqual(testCoeffs, self.testObjA.coeffs)

	@mock.patch("gen_basis_helpers.fit_cp2k_basis.cp2k_tables_updater.parseCP2KBasis.getCP2KBasisFromPlatoOrbitalGauPolyBasisExpansion")
	def testGetBasisInCP2kFormatHasCorrectCalls(self, mockedBasisConverter):
		testCoeffs = [1,2]
		self.testObjA.coeffs = testCoeffs
		expBasisPlato, expBasisCP2K, expAngMom = [mock.Mock()], mock.Mock(), 4
		expBasisPlato[0].label = expAngMom

		self.coeffToBasisMapperA.side_effect = lambda *args:expBasisPlato
		mockedBasisConverter.side_effect = lambda *args, **kwargs: expBasisCP2K

		actBasisCP2K = self.testObjA._getBasisSetCP2KFormat()
		self.coeffToBasisMapperA.assert_called_with(testCoeffs)
		mockedBasisConverter.assert_called_with(expBasisPlato, [expAngMom], self.eleName, basisNames=self.basisNames)
		self.assertEqual(expBasisCP2K, actBasisCP2K)

	@mock.patch("gen_basis_helpers.fit_cp2k_basis.cp2k_tables_updater.parseCP2KBasis.ParsedBasisFileCP2K")
	@mock.patch("gen_basis_helpers.fit_cp2k_basis.cp2k_tables_updater.BasisCoeffUpdaterCP2KStandard._getBasisSetCP2KFormat")
	def testGetCP2KBasisFullFileObj(self, mockedGetCP2KBasis, mockedFullFileCls):
		ourBasis, expOutput = mock.Mock(), mock.Mock()

		mockedGetCP2KBasis.side_effect = lambda: ourBasis
		mockedFullFileCls.side_effect = lambda *args,**kwargs: expOutput

		actOutput = self.testObjA._getCP2KBasisFullFileObj()
		mockedGetCP2KBasis.assert_called_with()
		mockedFullFileCls.assert_called_with(self.filePathA, [ourBasis])
		self.assertEqual(expOutput, actOutput)

	@mock.patch("gen_basis_helpers.fit_cp2k_basis.cp2k_tables_updater.parseCP2KBasis.writeBasisFileFromParsedBasisFileObj")
	@mock.patch("gen_basis_helpers.fit_cp2k_basis.cp2k_tables_updater.BasisCoeffUpdaterCP2KStandard._getCP2KBasisFullFileObj")
	def testUpdateUsingCurrentCoeffs(self, mockedGetParsedFileObj, mockedWriter):
		expParsedFileObj = mock.Mock()
		mockedGetParsedFileObj.side_effect = lambda : expParsedFileObj
		self.testObjA._updateUsingCurrentCoeffs()
		mockedWriter.assert_called_with(self.filePathA, expParsedFileObj)

