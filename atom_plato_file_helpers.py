#!/usr/bin/python3

import os
import subprocess
import sys
sys.path.append('/media/ssd1/rf614/Work/usr_scripts/coding/Plato_Analysis_Lib_Functions')
import plato_pylib.plato.mod_plato_inp_files as platoInp
from tbint_helpers import ChDir

''' Code to help create *.atm files easily in python3'''



def writeAtmFileFromDict(outPath, inpDict):
	#First set some boring terms to default values
	optDict = {"xc":"LDA",
	           "nGrid":1000,
	           "nLoops":1000,
	           "potMix":0.7,
	           "orbMix":0.7,
	           "verbosity":2,
	           "job":0
				}
	optDict = {k.lower():v for k,v in optDict.items()}

	#Update based on usr options
	dCopy = {k.lower():v for k,v in inpDict.items()}
	optDict.update(dCopy)

	#Translate to plato keywords
	functToFlag = {"lda":0,"pbe":1}
	mapKeyToAtm = {"nGrid".lower():"NumberOfGridPoints",
	               "nLoops".lower():"NumberOfLoops",
	               "potMix".lower():"PotentialMix",
	               "orbMix".lower():"OrbitalMix",
	               "PP".lower():"PseudoPotential",
	               "nOrbs".lower():"NumberOfOrbitals",
	               "rc": "CutoffRadius"}
	#Where only the keyword needs changing
	for key in mapKeyToAtm.keys():
		if key in optDict.keys():
			optDict[ mapKeyToAtm[key] ] = optDict.pop(key)

	#Set relevant values where key AND value need changing
	if "xc" in optDict.keys():
		optDict["XCFlag".lower()] = functToFlag[ optDict["xc"].lower() ]
		optDict.pop("xc")
	if "extPot".lower() in optDict.keys():
		optDict["ExternalPotential"] = optDict["extpot"].toPlatoStr()
		optDict.pop("extpot")
	if "orbInfo".lower() in optDict.keys():
		optDict["OrbitalList"] = "\n".join( [x.toPlatoStr() for x in optDict["orbinfo"]] )
		optDict.pop("orbinfo")
	if "coreInfo".lower() in optDict.keys():
		optDict["Core"] = optDict["coreinfo"].toPlatoStr()
		optDict.pop("coreinfo")

	#Convert everything to a str
	optDict = {k:str(v) for k,v in optDict.items()}

	platoInp.writePlatoOutFileFromDict(outPath, optDict)


def runBasis(basePath, **kwargs):
	kwargs = {k.lower():v for k,v in kwargs.items()}
	splitIndices = kwargs.get("splitIndices".lower(),None)
	nOrbs = kwargs.get("nOrbs".lower(),None)

	#Deal with filePath
	baseFolder, baseFile = os.path.split(basePath)
	baseFile = os.path.splitext(baseFile)[0]
	if baseFolder == "":
		baseFolder = os.getcwd()

	#Create the relevant string to run basis
	extraFlags = list()
	if splitIndices is not None:
		extraFlags.append("-sp")
	commStr = "basis " + " ".join(extraFlags) + " " + baseFile

	#Create file to split basis if needed
	if splitIndices is not None:
		splitOrbRespFile = "yes.txt"
		_createYesNoResponseFile(baseFolder, nOrbs, splitIndices, fName=splitOrbRespFile)

	#Run the program
	with ChDir(os.getcwd(),baseFolder):
		if splitIndices is not None:
			subprocess.check_call(commStr + "<{}".format(splitOrbRespFile) ,shell=True)
		else:
			subprocess.check_call(commStr,shell=True)
	


def _createYesNoResponseFile(folder, numbResponses:int, yesIndices:"list", fName=None):
	if fName is None:
		fName = "yes.txt"
	outVals = ["y" if x in yesIndices else "n" for x in range(numbResponses) ]

	outStr = "\n".join(outVals)
	outPath = os.path.join(folder , fName)

	with open(outPath,"wt") as f:
		f.write(outStr)
	return outPath

class OrbInfoAtmCode():
	def __init__(self, n:"int", l:"int", spinUp:"float, between 0 and 1 for occ. of spin-up", spinDown, write:"int 1 to write, 0 else"):
		self.n = int(n)
		self.l = int(l)
		self.spinUp = float(spinUp)
		self.spinDown = float(spinDown)
		self.write = int(write)
	def toPlatoStr(self):
		return "{} {} {:.4f} {:.4f} {}".format(self.n, self.l, self.spinUp, self.spinDown, self.write)

class ExtPotObj():
	def __init__(self, flag, v0=None, r0=None, potPower=None):
		self.flag = flag
		self.v0 = v0
		self.r0 = r0
		self.potPower = potPower
		#Input checking
		if flag > 0:
			assert v0 is not None, "if flag={} then v0 needs to be set".format(flag)
			assert r0 is not None, "if flag={} then r0 needs to be set".format(flag)
			assert potPower is not None, "if flag={} then r0 needs to be set".format(flag)

	def toPlatoStr(self):
		if self.flag==0:
			return "{}".format(self.flag)
		else:
			return "{}\n{:.4f} {:.4f} {:.4f}".format(self.flag, self.v0, self.r0, self.potPower)

class CoreInfo():
	def __init__(self, maxL, nMaxList):
		self.maxL = maxL
		self.nMaxList = list(nMaxList)
	def toPlatoStr(self):
		formStr = "{} " + " ".join(["{}" for x in range(len(self.nMaxList))])
		return formStr.format(self.maxL,*self.nMaxList)



