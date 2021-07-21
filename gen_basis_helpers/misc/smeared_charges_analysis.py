

import math

import numpy as np

import plato_pylib.shared.unit_convs as uConvHelp

import gen_basis_helpers.analyse_md.calc_dists as calcDistHelp



def getStandardExponentDict_Cardenas2016_neutralVals():
	""" Get dictionary of {element:exponent} using ATOMIC UNITS
	
	"""
	startDict = _getStandardHubbardVals_Cardenas2016_neutralVals()
	outDict = dict()
	prefactor = math.pi / 8
	for key in startDict.keys():
		outDict[key] = prefactor * ( (startDict[key]*uConvHelp.EV_TO_HA)**2 )

	return outDict

def _getStandardHubbardVals_Cardenas2016_neutralVals():
	""" Gets values in eV from DOI:10.1039/c6cp04533b Phys. Chem. Chem. Phys., 2016,18, 25721-25734 
	
	This uses the neutral atoms
	"""
	outDict = {"H":2*12.84, "O":2*12.16, "Mg":2*7.65}
	return outDict


def getCoulombInteractionEnergyStandard(inpGeom, exponentDict, charges, indicesA, indicesB, distLenConv=1, eConv=1):
	""" Gets Coulombic interaction energy between indicesA and indicesB
	
	Args:
		inpGeom: (plato_pylib UnitCell object) Contains the geometry
		exponentDict: (dict) Keys are element/kind symbols. Values are what we put in (exponent/pi)^{3/2} * e^{-exponent*r^2} charge distribution
		charges: (iter of floats) Atomic charges; units should be elementary charge (which is the easiest way anyway)
		indicesA: (iter of ints) Indices in group A, we calcualte interaction energy between group A and B (but not WITHIN each group)
		indicesB: (iter of ints) Indices in group B, we calcualte interaction energy between group A and B (but not WITHIN each group)
		distLenConv: (float) We multiply the distance matrix by this value. If input geometry is in angstrom this should usually be set to the angstrom->bohr conversion		
		eConv: (float) Default units are generally going to be Hartree
 
	Returns
		outVal: (float) The total interaction energy between atoms in indicesA and indicesB
 
	"""
	cartCoords = inpGeom.cartCoords
	exponents = [ exponentDict[coord[-1]] for coord in cartCoords ]

	coulombInteractionMatrix = getCoulombInteractionMatrix(inpGeom, exponents, charges, indicesA, indicesB, distLenConv=distLenConv)

	interactionMatrix = np.where( np.isnan(coulombInteractionMatrix), 0, coulombInteractionMatrix)
	outSum = np.sum( interactionMatrix )
	outVal = (0.5*outSum)*eConv

	return outVal


#0.5*sum(Matrix) is the interaction energy
def getCoulombInteractionMatrix(inpGeom, exponents, charges, indicesA, indicesB, distLenConv=1, minDist=1e-6):
	""" Calculates Coulombic energy between relevant atoms. Generally expect this to give results in atomic units (Hartree). Uses q_A*q_B*getIntegralTwoSphericalGaussianSmearedDensities(dist, alphaA, alphaB) to get the pairwise energies
	
	Args:
		inpGeom: (plato_pylib UnitCell object) Contains the geometry
		exponents: (iter of floats) Charges distrib. are described by (exponent/pi)^{3/2} * e^{-exponent*r^2}. These should generally be in atomic units (bohr)
		charges: (iter of floats) Atomic charges; units should be elementary charge (which is the easiest way anyway)
		indicesA: (iter of ints) Indices in group A, we calcualte interaction energy between group A and B (but not WITHIN each group)
		indicesB: (iter of ints) Indices in group B, we calcualte interaction energy between group A and B (but not WITHIN each group)
		distLenConv: (float) We multiply the distance matrix by this value. If input geometry is in angstrom this should usually be set to the angstrom->bohr conversion
		minDist: (float) Minimum distance to consider; at time of writing shorter distances are treated as having an interaction of zero. I may fix this to use the proper r->0 expression later. REGARDLESS: This function shouldnt be used for distances which are very near zero anyway.

	NOTE:
		I expect both exponents and charges to be given for ALL atoms in the cell, not just ones that are required to calculate interaction energies

	Returns
		outMatrix: (NxN np array) This contains the interaction energies between atoms in indicesA and indicesB
 
	"""
	#
	distMatrix = calcDistHelp.calcDistanceMatrixForCell_minImageConv(inpGeom, indicesA=indicesA, indicesB=indicesB)
	exponentsA = [exponents[idx] for idx in indicesA]
	exponentsB = [exponents[idx] for idx in indicesB]
	chargesA = [charges[idx] for idx in indicesA]
	chargesB = [charges[idx] for idx in indicesB]

	#
	distMatrixIndices =  [idx for idx in np.ndindex(distMatrix.shape)]
	outMatrix = np.zeros( (len(indicesA),len(indicesB)) )

	for arrayIdx in distMatrixIndices:
		idxA, idxB = arrayIdx
		currDist = distMatrix[idxA][idxB]*distLenConv
		alphaA, alphaB = exponentsA[idxA], exponentsB[idxB]
		chargeA, chargeB = chargesA[idxA], chargesB[idxB]
		currIntegral = getIntegralTwoSphericalGaussianSmearedDensities(currDist, alphaA, alphaB)
		outMatrix[idxA][idxB] = chargeA*chargeB*currIntegral

	return outMatrix


def getIntegralTwoSphericalGaussianSmearedDensities(dist, alphaA, alphaB):
	""" Calculates the Hartree integral (  \int \\fract{\\rho(r)\\rho(r')}{|r-r'|} dr dr' ) for two charge densities described by single spherical Gaussian functions. Used to get hartree energy by multiplying by q_{A}q_{B}.

	Note the formula is taken from the plato manual
	
	Args:
		dist: (float) Distance between the Gaussians
		alphaA: (float) Width parameter (i.e. exponent) of the first Gaussian
		alphaB: (float) Width parameter (i.e. exponent) of the second Gaussian
			 
	Returns
		outVal: (float) Value of the integral
 
	"""
	prefactor = (math.pi**3) / ( ((alphaA*alphaB)**(3/2))*dist )
	erfArg = math.sqrt( (alphaA*alphaB)/(alphaA+alphaB) ) * dist
	erfTerm = math.erf(erfArg)
	return prefactor*erfTerm


