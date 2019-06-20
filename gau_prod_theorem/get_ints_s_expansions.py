#!/usr/bin/python3

''' Purpose of this is simply to get the overlap of mcweda weight functions at rc*2 '''

import gau_prod_theorem as gProd
import itertools as it
import plato_pylib.plato.parse_gau_files as parseGau


def getSelfOverlapMcWedaWeightFromGauPolyBasis(gPolyBasis, dist:float):
	#Step 1 is to convert the basis into a specialised class
	allGauPrimsA = list()
	allGauPrimsB = list()
	for exp,coeff in it.zip_longest(gPolyBasis.exponents, gPolyBasis.r0Coeffs):
		allGauPrimsA.append( gProd.GauPrim.fromZDirOnly(exp,coeff,0) )
		allGauPrimsB.append( gProd.GauPrim.fromZDirOnly(exp,coeff,dist) )

	#Step two is to use these expansions to get the integral
	return _getIntegralTwoSOrbitalExpansions(allGauPrimsA,allGauPrimsB)


def _getIntegralTwoSOrbitalExpansions(expA:list, expB:list):
	sum = 0.0
	for gA, gB in it.product(expA,expB):
		gAB = gProd.combineTwoGauPrims(gA,gB)
		sum += gAB.getIntegralAllSpace()
	return sum

