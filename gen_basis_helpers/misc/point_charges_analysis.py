
import itertools as it

import numpy as np

from ..analyse_md import calc_dists as calcDistsHelp

def getCoulombEnergyBetweenIndicesForPointCharges(inpGeom, charges, indicesA, indicesB, lenConv=1):
	""" Gets the interaction energy between groups of point charges by using Coulombs law. By default we assume geometry using angstrom while charges are in units of elementary charge
	
	Args:
		inpGeom: (plato_pylib UnitCell Object) Contains the geometry, used to get distances between point charges
		indicesA: (iter of ints)
		indicesB: (iter of ints)
		charges: (iter of floats) Elementart charges for ALL atoms in the geometry
		lenConv: (float) Conversion factor for length; default is angstrom
 
	Returns
		outEnergy: (float) Estimated Coulombic interaction energies between indicesA and indicesB. Default units are eV
 
	Raises:
		ValueError: If theres any overlap between indicesA and indicesB

	"""

	#Check indicesA and indicesB are mutually exclusive;
	if len( set(indicesA).intersection( set(indicesB) ) )!=0:
		raise ValueError("Overlapping indices found; which isnt allowed") 

	#Calcylate the interactions
	coulombMatrix = getCoulombEnergyInteractionMatrix(inpGeom, charges, indicesA, indicesB, lenConv=lenConv)

	totalInts = 0
	for idxA,valA in enumerate(indicesA):
		for idxB,valB in enumerate(indicesB):
			totalInts += coulombMatrix[idxA][idxB]

	return totalInts


#Tested indirectly with getCoulombEnergyBetweenIndicesForPointCharges
def getCoulombEnergyInteractionMatrix(inpGeom, charges, indicesA, indicesB, lenConv=1):
	""" Gets a matrix with interaction energies between atom indices
	
	Args:
		inpGeom: (plato_pylib UnitCell Object) Contains the geometry, used to get distances between point charges
		indicesA: (iter of ints)
		indicesB: (iter of ints)
		charges: (iter of floats) Elementart charges for ALL atoms in the geometry
		lenConv: (float) Conversion factor for length; default is angstrom
 
	Returns
		outMatrix; (NxM matrix) Contains Coulombic interaction energies between indicesA and indicesB
 
	Notes:
		a)Theres a minimum distance (not settable at the mo) due to division by zero errors
		b) Default units are angstrom and eV
	"""
	#Get distance and charge matrices needed
	distMatrix = calcDistsHelp.calcDistanceMatrixForCell_minImageConv(inpGeom, indicesA=indicesA, indicesB=indicesB)
	chargesA = [ charges[idx] for idx in indicesA ]
	chargesB = [ charges[idx] for idx in indicesB ]

	chargesMatrix = np.zeros( list(distMatrix.shape) + [2] )
	for idxA, idxB in np.ndindex( distMatrix.shape[:2] ):
		chargesMatrix[idxA][idxB] = [chargesA[idxA], chargesB[idxB]]

	#
	minVal = 1e-9
	outMatrix = np.zeros( np.array(distMatrix).shape )
	for idx in np.ndindex(outMatrix.shape):
		qA, qB = chargesA[idx[0]], chargesB[idx[1]]
		if distMatrix[idx] < minVal:
			outMatrix[idx] = np.inf
		else:
			outMatrix[idx] = getCoulombEnergyTwoPointsStandard(qA, qB, distMatrix[idx], lenConv=lenConv)

	return outMatrix


def getCoulombEnergyTwoPointsStandard(qA, qB, dist, lenConv=1, energyConv=1):
	""" Gets the energy between charges at two points separated by dist, using Coulombs law
	
	Args:
		qA: (float) The charge on point A. Units=elementary charge
		qB: (float) The charge on point B. Units=elementary charge
		dist: (float) The distance between points A and B
		lenConv: (float) Conversion factor for length; default is angstrom
		energyConv: (float) Conversion factor for energy; default is eV
 
	Returns
		outE: (float) Energy from Coulombs law. Standard units are eV.
 
	"""
	coulombsConstant = 14.3996*lenConv # eV Angstrom e^{-2} where e is elementary charge
	outVal = coulombsConstant * ((qA*qB)/ dist)
	return outVal*energyConv






















