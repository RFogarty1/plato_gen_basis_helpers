
import math
import numpy as np

import MDAnalysis.lib.distances as distLib

from . import mdanalysis_interface as mdAnalysisInter

#Will be too slow; we need to convert ALL the cart-coords.
#def calcSingleDistBetweenCentralIndices_minImageConv(inpCell, idxA, idxB):
#	cartCoords = inpCell.cartCoords


def calcSingleDistBetweenCoords_minImageConv(inpCell, coordA, coordB):
	""" Gets the distance between two co-ordinates in inpCell using the minimum image convention; thus, this will return the smallest distance possible with the PBCs
	
	Args:
		inpCell: (plato_pylib UnitCell object) We use this to get the lattice vectors needed to account for PBCs
		coordA: (len 3 iter) Co-ordinates of first atom
		coordB: (len 3 iter) Co-ordinates of second atom
 
	Returns
		dist: (float) The distance between the co-ordinates after taking PBCs into account
 
	"""
	dims = mdAnalysisInter.getMDAnalysisDimsFromUCellObj(inpCell)
	dists = distLib.calc_bonds(np.array([coordA[:3]]),np.array([coordB[:3]]),box=dims)
	assert len(dists)==1
	return dists[0]


def calcSingleAngleBetweenCoords_minImageConv(inpCell, coordA, coordB, coordC):
	""" Gets the angle (ABC) between three co-ordinates in inpCell using the minimum image convention; thus, this will use the closest possible A,B,C versions
	
	Args:
		inpCell: (plato_pylib UnitCell object) We use this to get the lattice vectors needed to account for PBCs
		coordA: (len 3 iter) Co-ordinates of first atom
		coordB: (len 3 iter) Co-ordinates of second atom
		coordC: (len 3 iter) Co-ordinates of third atom
 
	Returns
		angle: (float) The angle (in degrees) between the co-ordinates after taking PBCs into account
 
	"""
	dims = mdAnalysisInter.getMDAnalysisDimsFromUCellObj(inpCell)
	args = [ np.array([coordA[:3]]), np.array([coordB[:3]]), np.array([coordC[:3]]) ]
	angles = distLib.calc_angles(*args, box=dims)
	assert len(angles)==1
	return math.degrees(angles[0])

