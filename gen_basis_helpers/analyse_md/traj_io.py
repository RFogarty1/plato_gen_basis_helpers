
import re

import plato_pylib.shared.ucell_class as uCellHelp

from . import traj_core as trajCore


def writeTrajToSimpleExtendedXyzFormat(traj, outPath):
	outList = list()
	for trajStep in traj:
		currSection = _getFileAsListForSingleTraj(trajStep)
		outList.extend(currSection)

	outStr = "\n".join(outList)
	with open(outPath,"wt") as f:
		f.write(outStr)

def _getFileAsListForSingleTraj(traj):
	outList = [ str( len(traj.unitCell.cartCoords) ) ]

	#Get comment line
	lattVectStr = ""
	for lVect in traj.unitCell.lattVects:
		lattVectStr += " ".join([str(x) for x in lVect]) + " "

	stepStr = " step={}".format(str(traj.step)) if traj.step is not None else ""
	timeStr = " time={:.2f}".format(traj.time) if traj.time is not None else ""
	commentLine = "Lattice=\"{}\"".format(lattVectStr)
	commentLine += " Properties=species:S:1:pos:R:3"
	commentLine += stepStr + timeStr
	outList.append(commentLine)

	#Get the cart coords
	for coord in traj.unitCell.cartCoords:
		currStr = "{} {:.8f} {:.8f} {:.8f}".format( coord[-1], *coord[0:3] )
		outList.append(currStr)

	return outList

def readTrajFromSimpleExtendedXyzFormat(inpPath):
	fileAsList = _readFileIntoList(inpPath)

	lineIdx = 0
	outSteps = list()
	while lineIdx<len(fileAsList):
		lineIdx, currStep = _parseNextExtendedXyzSection(fileAsList, lineIdx)
		outSteps.append(currStep)

	outObj = trajCore.TrajectoryInMemory(outSteps)
	return outObj

def _parseNextExtendedXyzSection(fileAsList, lineIdx):
	nAtoms = int( fileAsList[lineIdx].strip() )
	lineIdx += 1

	#Parser the comment line; which holds latt params and possibly step/time
	currTrajStep = _parseExtendedXyzCommentLine(fileAsList[lineIdx])
	lineIdx += 1

	#parse the cart coords
	cartCoords = list()
	for idx in range(nAtoms):
		currLine = fileAsList[lineIdx].strip().split()
		currCoords = [float(x) for x in currLine[1:4]] + [currLine[0]]
		cartCoords.append(currCoords)
		lineIdx+=1

	currTrajStep.unitCell.cartCoords = cartCoords
	return lineIdx, currTrajStep

def _parseExtendedXyzCommentLine(commentLine):
	#1) Get the unit cell parameters
	pattern = "Lattice=\".*\""
	lattStr = re.findall(pattern, commentLine)[0].replace("Lattice=","").replace("\"","")
	lattVals = [float(x) for x in lattStr.split()]
	assert len(lattVals)==9
	lattVects = [ lattVals[:3], lattVals[3:6], lattVals[6:9] ]
	unitCell = uCellHelp.UnitCell.fromLattVects(lattVects)

	#Get a step if possible
	pattern = "step=[0-9]+"
	currMatches = re.findall(pattern, commentLine)
	if len(currMatches)==0:
		step = None
	elif len(currMatches)==1:
		step = int( currMatches[0].replace("step=","") )
	else:
		raise ValueError("Found {} matches for step on line {}".format(len(currMatches), commentLine))

	#Get a time if possible
	pattern = "time=[0-9\.]+"
	currMatches = re.findall(pattern, commentLine)
	if len(currMatches)==0:
		time = None
	elif len(currMatches) == 1:
		time = float( currMatches[0].replace("time=","") )

	#Merge all into a trajStep
	outObj = trajCore.TrajStepBase(unitCell=unitCell,step=step,time=time)
	return outObj


def _readFileIntoList(inpPath):
	with open(inpPath,"rt") as f:
		fileAsList = f.readlines()
	return fileAsList

