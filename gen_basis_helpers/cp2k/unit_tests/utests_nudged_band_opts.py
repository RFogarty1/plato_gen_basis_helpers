
import copy
import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.cp2k_nudged_band_options as tCode

class TestNudgedBandReplicas(unittest.TestCase):

	def setUp(self):
		self.coordsStart = [ [1,2,3],[4,5,6] ]
		self.coordsEnd   = [ [6,7,8],[7,8,9] ] 
		self.coordsInter = None
#		self.lenConvFactor = 1
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.NudgedBandReplicas(self.coordsStart, self.coordsEnd, interCoords=self.coordsInter)

#		self.testObjA = tCode.NudgedBandReplicas(self.coordsStart, self.coordsEnd, interCoords=self.coordsInter,
#		                                         lenConvFactor=self.lenConvFactor)

	def testExpectedOptDict_startAndEndDefined(self):
		expReplicaCoords = [self.coordsStart,self.coordsEnd]
		expOutDict = {"nudgedband_replica_coords":expReplicaCoords}
		actOutDict = self.testObjA.optDict
		for key in expOutDict.keys():
			self.assertEqual(expOutDict[key], actOutDict[key])

	def testExpectedOptDict_startAndEndDefined_withSuperfluosAtomicSymbols(self):
		expReplicaCoords = [copy.deepcopy(x) for x in [self.coordsStart,self.coordsEnd]]
		expOutDict = {"nudgedband_replica_coords":expReplicaCoords}

		self.coordsStart = [x+["H"] for x in self.coordsStart]
		self.coordsEnd = [x+["O"] for x in self.coordsEnd]
		self.createTestObjs()

		actOutDict = self.testObjA.optDict
		self._checkExpAndActDictsEqual(expOutDict, actOutDict)

	def testExpectedOptDict_twoInters(self):
		interCoordsA = [ [3,2,1], [2,3,1] ]
		interCoordsB = [ [5,1,3], [4,6,2] ]
		self.coordsInter = [interCoordsA, interCoordsB]
		self.createTestObjs()

		expOutDict = {"nudgedband_replica_coords":[self.coordsStart, interCoordsA, interCoordsB, self.coordsEnd]}
		actOutDict = self.testObjA.optDict

		self._checkExpAndActDictsEqual(expOutDict, actOutDict)

	@unittest.skip("Will try and force the cp2k side to use bohr...")
	def testExpectedOptDict_convFactorApplied(self):
		self.lenConvFactor = 2

		#startCoords
		startCoords = list()
		for sCoords in self.coordsStart:
			currCoords = [self.lenConvFactor*x for x in sCoords]
			startCoords.append( currCoords )

		#endCoords
		endCoords = list()
		for eCoords in self.coordsEnd:
			currCoords = [self.lenConvFactor*x for x in eCoords]
			endCoords.append( currCoords)


#		startCoords = [ self.lenConvFactor*self.coordsStart ]

		self.assertTrue(False)

	def _checkExpAndActDictsEqual(self, expDict, actDict):
		for key in expDict.keys():
			self.assertEqual(expDict[key], actDict[key])


class TestNudgedBandOpts(unittest.TestCase):

	def setUp(self):
		self.numbReplicas = 4
		self.procsPerReplicaEnv = 8
		self.springConstant = None
		self.nebType = "CI-NEB"
		self.alignFrames = False
		self.rotateFrames = False
		self.printInitConfigInfo = True
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"numbReplicas":self.numbReplicas, "procsPerReplicaEnv":self.procsPerReplicaEnv,
		             "springConstant":self.springConstant, "nebType":self.nebType, "alignFrames":self.alignFrames,
		             "rotateFrames":self.rotateFrames, "printInitConfigInfo":self.printInitConfigInfo}
		self.testObjA = tCode.NudgedBandOptsStd(**kwargDict)

	def testExpectedOptDictA(self):
		expOptDict = {"nudgedband_numbReplicas":self.numbReplicas, "nudgedband_procsPerReplicaEnv":self.procsPerReplicaEnv,
		              "nudgedband_type":self.nebType, "nudgedband_alignFrames":self.alignFrames,
		              "nudgedband_rotateFrames":self.rotateFrames, "nudgedband_printInitConfigInfo":self.printInitConfigInfo}
		actOptDict = self.testObjA.optDict

		for key in expOptDict.keys():
			self.assertEqual( expOptDict[key], actOptDict[key] )




