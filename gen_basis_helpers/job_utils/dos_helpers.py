

import math
import os
import subprocess
import sys

import numpy as np

import analyse_plato_bandstructures as bStructs
import plato_pylib.parse_plato_out_files as platoOut
from ..shared import ch_dir as chDir 



def getDosPlotData(inpFile:"Path to input file", smearing:"eV"=0.1, stepSize:"eV"=0.02, runDos=True):
	if inpFile.endswith(".occ"):
		return _getDosInfoPlato(inpFile, smearing=smearing, runDos=runDos)
	else:
		raise ValueError("{} is an invalid path for getDosPlotData".format(inpFile))




def _getDosInfoPlato(inpFile, smearing=0.1, stepSize=0.01, runDos=True):
	#Figure out number of points needed to get ~stepSize
	extraSigmaWidth = smearing*5.0 #determines how far we plot energies relative to smearing and min/max eigenvals; val taken from plato src code
	parsedFile = platoOut.parseOccFile(inpFile)
	minE,maxE =  np.min(parsedFile["eigen_vals"]) - extraSigmaWidth , np.max(parsedFile["eigen_vals"]) + extraSigmaWidth
	nPts = ( abs(maxE - minE) ) / stepSize

	# Figure out the run command + directory to change to[TODO: call function below instead once sure it works]
	folder, fileNameWithExt = os.path.split(inpFile)
	if folder == "":
		folder = os.getcwd()
	fileName, unused = os.path.splitext(fileNameWithExt)
	runComm = "plotdos -n {} -s {} {}".format(nPts, smearing, fileName)

	# Run the command + get the output data
	with chDir.ChDir(os.getcwd(), folder):
		if runDos:
			subprocess.check_call(runComm,shell=True)
		outDict = _parsePlatoPlotDosOutput(fileName + "-dos.csv")
	return outDict


def getDosRunComm_plato(inpFile, smearing:"ev", stepSize:"eV"):
	inpFile = os.path.splitext(inpFile)[0] + ".occ"

	extraSigmaWidth = smearing*5.0 #determines how far we plot energies relative to smearing and min/max eigenvals; val taken from plato src code
	parsedFile = platoOut.parseOccFile(inpFile)
	minE,maxE =  np.min(parsedFile["eigen_vals"]) - extraSigmaWidth , np.max(parsedFile["eigen_vals"]) + extraSigmaWidth
	nPts = ( abs(maxE - minE) ) / stepSize

	folder, fileNameWithExt = os.path.split(inpFile)
	folder = os.path.abspath(folder)
	fileName, unused = os.path.splitext(fileNameWithExt)
	runComm = "cd {};plotdos -n {} -s {} {}".format(folder, math.ceil(nPts), smearing, fileName)

	return runComm


def _parsePlatoPlotDosOutput(filePath):
	outDict = dict()
	with open(filePath,"rt") as f:
		fileAsList = f.readlines()

	energies, dos = list(), list()

	#Note, 1st line is a header
	for idx,line in enumerate(fileAsList[1:]): 
		currLine = line.strip().split(",")
		if len(currLine) == 1:
			infoIdx = idx+2 #Line after a blank gives fermi energy (+2 is because we skipped the 0th line)
			break 
		else:
			energies.append( float(currLine[0]) )
			dos.append( float(currLine[1]) )

	eFermi = float( fileAsList[infoIdx].strip().split(",")[0] )	
	dosData = list()
	for energy,dosVal in zip(energies,dos):
		dosData.append( (energy,dosVal) )

	#Put vals in outDict
	outDict["dosData".lower()] = dosData
	outDict["eFermi".lower()] = eFermi

	return outDict



def parseOptaDosDatFile(filePath):
	outDict = dict()
	energies, dos = list(), list()

	with open(filePath,"rt") as f:
		fileAsList = f.readlines()

	for line in fileAsList:
		currLine = line.strip().split()
		if not(currLine[0].startswith("#")):
			energies.append( float(currLine[0]) )
			dos.append( float(currLine[1]) )

	dosData = list()
	for eVal,dosVal in zip(energies,dos):
		dosData.append( (eVal,dosVal) )

	outDict["dosData".lower()] = dosData
	return outDict


def getShiftedDosData(dosData, shiftVal):
	outDosData = list()
	for x,y in dosData:
		outDosData.append( (x+shiftVal,y) )
	return outDosData

