#!/usr/bin/python3

"""Purpose of this module is to convert my types into 
relevant ASE types
"""

from ase.atoms import Atoms


def getAseAtomsObjFromPylibUCell(inpUCell):
	scaledPositions = [x[:-1] for x in inpUCell.fractCoords] #We remove the atom symbol from the end
	symbols = inpUCell._elementList #Will be None if never set, which is fine
	lattVects = inpUCell.lattVects
	pbcs = True
	outObj = Atoms(scaled_positions=scaledPositions, symbols=symbols, cell=lattVects, pbc=pbcs) 

	print("outObj = {}".format(outObj))
	return outObj

