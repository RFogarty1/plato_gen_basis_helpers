#!/usr/bin/python3



import sys
sys.path.append('..')
import unittest

import numpy as np
import inv_sk_helpers as tCode

sys.path.append('/media/ssd1/rf614/Work/usr_scripts/coding/Plato_Analysis_Lib_Functions')
import parse_inv_sk as parseInvSk
import plato_pylib.plato.parse_tbint_files as parseTbint

class testGetDataInvSK(unittest.TestCase):

	def setUp(self):
		self.invSKObjA = createMockInvSKObjA()
		self.paresdTbintFileA = createMockParsedTbintA()
	def tearDown(self):
		pass


	def testDistanceVsHopVal(self):
		testShellA, testShellB = 0, 0
		testXProp, testYProp = "distance", "hval"
		testBondType = "sigma"


		expectedVals = np.array(( [2.0,-2.5],[4.0,-0.5] ))
		actualVals = tCode.getInvSKDataTwoShells(self.invSKObjA, testShellA, testShellB, testBondType,
		                                         testXProp, testYProp)

		self.assertTrue( np.allclose(expectedVals, actualVals) )

	def testScreenFunctVsHopVal(self):
		testShellA, testShellB = 1,1
		testXProp, testYProp = "screenfunct", "hval"
		testBondType = "sigma"
		expScreenVals = [8.5,10.2]
		expHVals = [-1.5,0.2]

		expectedVals = np.array(( [expScreenVals[0],expHVals[0]],[expScreenVals[1],expHVals[1]] ))
		actualVals = tCode.getInvSKDataTwoShells(self.invSKObjA, testShellA, testShellB, testBondType,
		                                         testXProp, testYProp)

		self.assertTrue( np.allclose(expectedVals, actualVals) )

	def testDistanceVsHerror(self):
		testShellA, testShellB = 0,0
		testXProp, testYProp = "distance", "herror"
		testBondType = "sigma"
		expX, expY = [2.0,4.0], [3.0,1.7]

		expectedVals = np.array(( [expX[0],expY[0]], [expX[1], expY[1]] ))
		actualVals = tCode.getInvSKDataTwoShells(self.invSKObjA, testShellA, testShellB, testBondType,
		                                         testXProp, testYProp, tbintDataObjs=self.paresdTbintFileA)

		self.assertTrue( np.allclose(expectedVals, actualVals) )


	def testScreenFunctTimesHr(self):
		testShellA, testShellB = 0,0
		testXProp, testYProp = "screenfunct_times_hr", "herror"
		testBondType = "sigma"

		expX, expY = [8.5*0.5,10.2*1.2], [3.0,1.7]
		expectedVals = np.array(( [expX[0],expY[0]], [expX[1], expY[1]] ))
		actualVals = tCode.getInvSKDataTwoShells(self.invSKObjA, testShellA, testShellB, testBondType,
		                                         testXProp, testYProp, tbintDataObjs=self.paresdTbintFileA)

		self.assertTrue( np.allclose(expectedVals, actualVals) )


def createMockInvSKObjA():
	fakePos = [0.0,0.0,0.0]
	fakeLVal = 0
	distVals = [2.0,2.0,4.0,4.0]
	shellVals = [(0,0),(1,1),(0,0),(1,1)]
	hVals = [ -2.5, -1.5, -0.5, 0.2 ]
	screenFunctVals = [8.5,8.5, 10.2, 10.2]
	allObjs = list()	

	for idx,unused in enumerate(distVals):
		currObj = parseInvSk.InvSKField(posA=fakePos,posB=fakePos,shellA=shellVals[idx][0],
		                                shellB=shellVals[idx][1],lA=fakeLVal,lB=fakeLVal,
		                                dist=distVals[idx],hValSigma=hVals[idx],
		                                screenFunct=screenFunctVals[idx])
		allObjs.append(currObj)
	return parseInvSk.InvSKAllData(allObjs)


def createMockParsedTbintA():
	outDict = {}

	allHopInts = list()
	rVals = [1.9,4.1]
	intVals_ss_sigma = np.array(( [rVals[0], 0.5],[rVals[1],1.2] ))
	intVals_pp_sigma = np.array(( [rVals[0], 0.7],[rVals[1],1.5] ))

	ssSigmaHop = parseTbint.TbintIntegrals(shellA=0,shellB=0,orbSubIdx=1,
	                            integrals=intVals_ss_sigma, intType="hopping")
	ppSigmaHop = parseTbint.TbintIntegrals(shellA=1,shellB=1,orbSubIdx=1,
	                            integrals=intVals_pp_sigma, intType="hopping")

	allHopInts.append(ssSigmaHop)
	allHopInts.append(ppSigmaHop)

	outDict = {"hopping":allHopInts}

	return outDict

if __name__ == '__main__':
	unittest.main()




