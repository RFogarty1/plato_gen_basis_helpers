

import unittest

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.atom_combo_opts_obj_maps as optObjMaps
import gen_basis_helpers.analyse_md.classification_distr_opt_objs as classDistrOptObjHelp



class TestWaterCountTypesMinDistAndHBond(unittest.TestCase):

	#Taken from TestDiscHBondCounterBetweenGroupsOxyDistFilter mostly
	def setUp(self):
		#1) All geometric parameters for testing
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#roll 90,pitch=45, OH len~1, HOH angle 104.5. Then just added translation vectors
		self.waterACoords = [ [0,0,0,"O"], [-0.13,0,0.99,"H"], [0.99,0,-0.13,"H"] ]
		#translation = [2,0,0]; no rotations
		self.waterBCoords = [ [2,0,0,"O"], [2.61, 0.79, 0, "H"], [2.61,-0.79,0,"H"] ]
		#translation = [2+(2*0.61), 2*0.79, 0]
		self.waterCCoords = [ [3.22, 1.58, 0, "O"],  [3.83, 2.37, 0, "H"], [3.83, 0.79, 0, "H"] ]
		#translation = [2+(2*0.61), -2*0.79,0]
		self.waterDCoords = [ [3.22, -1.58, 0, "O"], [3.83, -0.79, 0, "H"], [3.83, -2.37, 0, "H"] ]

		self.xCoord = [[0,0,0,"X"]]
		self.coords = self.waterACoords + self.waterBCoords + self.waterCCoords + self.waterDCoords + self.xCoord

		#Options
		self.oxyIndices = [0,3,6,9]
		self.hyIndices = [ [1,2], [4,5], [7,8], [10,11] ]
		self.distFilterIndices = [12]
		self.distFilterVals = [ [0,3], [3,5] ]#AB should be one group, with CD as the other. + easy to flip this
		self.maxOO = 3 #AC h-bond would be possible iff this was set high enough i suspect
		self.maxAngle = 35

		self.binResObjs = [None,None] #Irrelevent to these tests so....

		self.nDonorFilterRanges = [ [-0.5,1.5], [0.5,2.5] ]
		self.nAcceptorFilterRanges = [ [-0.5,0.5], [0.5,2.5] ]
		self.nTotalFilterRanges = [ [-0.5,1.5], [0.5,2.5] ]
		self.checkInputConsistent = True

		self.createTestObjs()

	def createTestObjs(self):
		#Sort the geometry out
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create an options object
		currArgs = [self.binResObjs, self.oxyIndices, self.hyIndices, self.distFilterIndices, self.distFilterVals]
		currKwargs = {"nDonorFilterRanges": self.nDonorFilterRanges, "nAcceptorFilterRanges": self.nAcceptorFilterRanges,
		              "nTotalFilterRanges": self.nTotalFilterRanges, "maxOOHBond": self.maxOO, "maxAngleHBond": self.maxAngle,
		              "checkInputConsistent": self.checkInputConsistent}
		self.optObj = classDistrOptObjHelp.WaterCountTypesMinDistAndHBondSimpleOpts(*currArgs,**currKwargs)

		#Get sparse matrix populator + populate it
		self.sparseMatrixCalculator = optObjMaps.getSparseMatrixCalculatorFromOptsObjIter([self.optObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binner object
		self.testObj = optObjMaps.getMultiDimBinValGetterFromOptsObjs([self.optObj])

	def testExpectedA(self):
		expBinVals = [ (1,0) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

	def testExpectedB_diffRanges(self):
		#Put first 2 in first bin, second 2 in 2nd bin
		self.nDonorFilterRanges    =  [  [-0.5,3.5], [-0.5,0.5] ]
		self.nAcceptorFilterRanges =  [  [-0.5,1.5], [0.5,1.5 ] ]
		self.nTotalFilterRanges    =  [  [0.5,3.5 ], [0.5,1.5 ] ]
		self.createTestObjs()

		expBinVals = [ (2,2) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals,actBinVals)



