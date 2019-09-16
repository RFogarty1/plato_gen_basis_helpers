#!/usr/bin/python3

"""Test that we correctly convert plato_pylib unitCell class into ASE atoms class


"""

import itertools as it
import math
import sys
import unittest

sys.path.append("..")

import plato_pylib.shared.ucell_class as UCell
from ase.atoms import Atoms
import ase.build.bulk as aseBuilder

import ase_conversions as tCode

class TestForMgHcp(unittest.TestCase):

	def setUp(self):
		lattVects = ( [6.06,0.0,0.0], [-3.0299495954, 5.2480557244,0.0], [0.0,0.0,9.84] )
		fractPos = ( [0.0, 0.0, 0.0], [0.33333300, 0.66666700, 0.50000000] )
		atomLabels = ["Mg","Mg"]
		fractCoords = [ list(x)+[y] for x,y in it.zip_longest(fractPos,atomLabels) ]
		self.testUCell = UCell.UnitCell.fromLattVects(lattVects, fractCoords=fractCoords)
		self.aseCell = Atoms(pbc=True, cell=lattVects, symbols=atomLabels, scaled_positions = fractPos)

	def testAseAtomsObjCorrect(self):
		convertedAseCell = tCode.getAseAtomsObjFromPylibUCell(self.testUCell)
		self.assertTrue( checkAtomsObjsEqual(self.aseCell,convertedAseCell,floatTol=1e-7) )


def checkAtomsObjsEqual(objA, objB, floatTol=1e-8):
	""" Atoms Objects dont account for floating point imprecision, hence they need overriding """

	if (objA==objB):
		return True

	if (not isinstance(objA,Atoms) or (not isinstance(objB,Atoms))):
		return False

	a = objA.arrays
	b = objB.arrays

	#Mostly duplicating tests, but adding a float tolerance to the positions case
	if len(objA)!=len(objB):
		return False

	posDists = [getDistTwoVects(x,y) for (x,y) in it.zip_longest(a["positions"],b["positions"])]

	if not all( [x<floatTol for x in posDists] ):
		return False

	if not (a["numbers"] == b["numbers"]).all():
		return False


	#TODO: Likely need to do a try:Except to handle the case where there is no cell for full generality
	for vectA,vectB in it.zip_longest(objA.cell,objB.cell):
		diffs = [x-y for (x,y) in it.zip_longest(vectA,vectB)]
		if not all ( [x<floatTol for x in diffs] ):
			return False

	#Imperfect test. [True,True,True] != True at the moment. Need to isinstance(iter) likely
	if objA.pbc.all() != objB.pbc.all():
		return False

	return True


def getDistTwoVects(x,y):
	diffSqrVect = [(a-b)**2 for a,b in it.zip_longest(x,y)]
	return math.sqrt( sum(diffSqrVect) )




if __name__ == '__main__':
	unittest.main()


