#!/usr/bin/python3

import os
import shutil
import subprocess
import sys
import plato_pylib.plato.mod_plato_inp_files as platoInp
import plato_pylib.plato.parse_gau_files as parseGau

from ..shared import ch_dir as chDir

def getDefaultFitGaussDict():
	''' Returns dict of default keyword:string used to make *.gin file. '''
	''' Using a function (vs global) so that it doesnt get altered between calls '''

	defDict = {"Orbitals":[0],
	           "Density":[0, 1],
	           "NeutralAtom":[0, 1],
	           "OverlapTolerance": 10.0,
	           "MinimumRatio": 1.2,
	           "MinimumResolution":0.05,
	           "NonLocalPotential":[[0, 0]],
	           "weightfuncts":[0]}
	defDict = {k.lower():v for k,v in defDict.items()}
	return defDict


def getPlatoStrDictFromInpDict(inpDict):
	platDict = {k.lower():v for k,v in inpDict.items()}
	
	platDict["orbitals"] = " ".join( [str(x) for x in inpDict["orbitals"]] )
	platDict["density"] = " ".join( [str(x) for x in inpDict["density"]] )
	platDict["NeutralAtom".lower()] = " ".join( [str(x) for x in inpDict["NeutralAtom".lower()]] )
	platDict["weightfuncts"] = " ".join( [str(x) for x in inpDict["weightfuncts"]] )

	nlStr = ""
	for proj in inpDict["NonLocalPotential".lower()]:
		nlStr += " ".join([str(x) for x in proj]) + "\n"
	platDict["NonLocalPotential".lower()] = nlStr

	for key in platDict.keys():
		platDict[key] = str(platDict[key])

	return platDict

def writeGinFileFromDict(outPath,inpDict):
	#set default values
	optDict = getDefaultFitGaussDict()

	#Update + get in nice format for plato
	inpDict = {k.lower():v for k,v in inpDict.items()}
	optDict.update(inpDict)
	outDict = getPlatoStrDictFromInpDict(optDict)

	#write output
	platoInp.writePlatoOutFileFromDict(outPath,outDict)

def fitFromGinDict(basePath:"Path to *.bas file w/o extension", appendStr, optDict, maxFitRepeats=0)->"gauFilePath":

	#Create files
	newGauPath,newGinPath =  _createFilesForFitFromGinDict(basePath, appendStr, optDict)

    #Run fitgauss
	runFitGauss(newGinPath,maxFitRepeats)

	return newGauPath


def _createFilesForFitFromGinDict(basePath:"Path to *.bas file w/o extension", appendStr, optDict):
    basePath = os.path.splitext(basePath)[0] #Just incase path included extension
    oldBasPath = basePath + ".bas"    
    newGinPath = basePath + appendStr + ".gin"
    newBasPath = newGinPath.replace(".gin",".bas")
    newGauPath = newGinPath.replace(".gin",".gau")

    shutil.copy2(basePath+".bas", newBasPath)
    writeGinFileFromDict(newGinPath, optDict)

    return newGauPath,newGinPath


def runFitGauss(ginPath, nYes=0):
	basePath = os.path.splitext(ginPath)[0]
	workFolder, fileName = os.path.split(basePath)

	with chDir.ChDir(os.getcwd(),workFolder): #Context manager to temporarily step into a directory
		yesPath = createYesFile(nYes)
		subprocess.check_call( "fitgauss {} < {}".format(fileName,yesPath) , shell=True)
	outPath = os.path.join(workFolder,fileName + ".gau")

	return outPath




def createYesFile(nYes):
	yesStr = "\n".join(["y" for x in range(nYes)]) + "\nn\n" 
	filePath = "yes.txt"
	f = open("yes.txt",'w')
	f.write(yesStr)
	f.close()
	return filePath


def getFitFromPyFitProgOut(filePath):
	with open(filePath,"r") as f:
		fileAsListWithComments = f.readlines()
	fileAsList = list()
	for line in fileAsListWithComments:
		if line.find('#') == -1:
			fileAsList.append(line)
	parsedFit,lIdx = parseGau.parseFit(fileAsList,0)
	return parseGau.GauPolyBasis.fromIterable(parsedFit)


