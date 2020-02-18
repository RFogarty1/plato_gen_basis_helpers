
import itertools as it
import math
import numpy as np

def getVHaFromSphericDensityOnOneDimGrid(inpGrid, dNu):
	""" Calculate the Hartree potential for a spherical density when given a grid
	
	Args:
		inpGrid:(nx2) tuple list (e.g. [(1,2),(2,4)]). x-values are distances from origin, y values are densities at those points
		dNu: A float found in the *.bas output file. I have no idea what it means	
	Returns
		outGrid:(nx2) tuple list; x-values are distances from origin, y values are densities at those points

	Note: Algorithm just stolen out of plato basis program
	
	"""

	#I have no idea what these parameter names mean
	xVals = [x[0] for x in inpGrid]
	yVals = [x[1] for x in inpGrid]

	workA = [4*math.pi*x[0]*x[1] for x in inpGrid]
	workB = [4*math.pi*x[0]*x[0]*x[1] for x in inpGrid]

	derivWorkA = _diff(xVals, workA, dNu)
	derivWorkB = _diff(xVals, workB, dNu)

	#Calculate q values, which are used to get the hartree potential
	qVals, derivQVals = [0.0], [0.0]

	for idx in range(1,len(xVals)):
		currQ = qVals[idx-1] + _integrateCubic(xVals, workB, derivWorkB, idx-1, idx)
		currDerivQ = derivQVals[idx-1] + _integrateCubic(xVals, workA, derivWorkA, idx-1, idx) 
		qVals.append(currQ), derivQVals.append(currDerivQ)

	#Now get the actual hartree potential values (hopefully)
	vHaVals = list()
	for idx,y in enumerate(yVals):
		currVal = 2 * ( derivQVals[len(yVals)-1] - derivQVals[idx] + (qVals[idx]/xVals[idx]) )
		vHaVals.append(currVal)

	return [(x,y) for x,y in it.zip_longest(xVals,vHaVals)]


#Probably an anoying finite difference thing as the source of the magic numbers (e.g. 25)
def _diff(xVals, yVals, dNu):

	firstVal = (-25.0* yVals[0]/12) + (4*yVals[1]) - (3*yVals[2]) + (4*yVals[3]/3) - (yVals[4]/4)
	secondVal = (-yVals[0]/4) - (5*yVals[1]/6) + (3*yVals[2]/2) - (yVals[3]/2.0) + (yVals[4]/12)

	outVals = [firstVal, secondVal]
	
	for idx,y in enumerate(yVals[2:],2):
		if idx < len(yVals) - 2:
			currVal = (yVals[idx-2] / 12) - (2*yVals[idx-1]/3) + (2*yVals[idx+1]/3) - (yVals[idx+2]/12)
			outVals.append(currVal)

	n = len(yVals)
	nMinus2 = (-yVals[n-5]/12) + (yVals[n-4]/2) - (3*yVals[n-3]/2) + (5*yVals[n-2]/6) + (yVals[n-1]/4)
	nMinus1 = (yVals[n-5]/4) - (4*yVals[n-4]/3) + (3*yVals[n-3]) - (4*yVals[n-2]) + (25*yVals[n-1]/12)
	outVals.append(nMinus2)
	outVals.append(nMinus1)


	for idx,y in enumerate(yVals):
		outVals[idx] = outVals[idx] / (xVals[idx]*dNu)

	return outVals


def _integrateCubic(xVals, yVals, dyVals, firstIdx, lastIdx):
	total = 0.0

	for idx in range(firstIdx, lastIdx):
		deltaX = xVals[idx+1] - xVals[idx]
		deltaYPositive = yVals[idx] + yVals[idx+1]
		deltaDerivY = dyVals[idx] - dyVals[idx+1]
		total += deltaX  *  ( (6*deltaYPositive) + (deltaX*deltaDerivY)  ) / 12 

	return total


