
import os
import pathlib

import plato_pylib.plato.mod_plato_inp_files as platoInp



#Method class should encapsulate all the things that vary between different methods in plato
class PlatoMethod():

	def __init__(self, optDict, runCommFunction, strDictFromOptDictFunction, gridKwarg=None):
		self.optDict = optDict
		self.runCommFunction = runCommFunction #interface is getRunComm(fileName)
		self.getStrDictFromOptDict = strDictFromOptDictFunction #interface is getStrDictFromOptDict
		self.gridKwarg = "IntegralMeshSpacing".lower() if gridKwarg is None else gridKwarg #Different for dft case annoyingly
#		self.gridKwarg = "IntegralMeshSpacing".lower()

	@property
	def kpts(self):
		return self.optDict["BlochStates".lower()]

	@kpts.setter
	def kpts(self,value):
		self.optDict["BlochStates".lower()] = value

	@property
	def dataSet(self):
		return self.optDict["dataset"]

	@dataSet.setter
	def dataSet(self,value:"str, plato relative path"):
		self.optDict["dataset"] = value

	@property
	def integGrid(self):
		return self.optDict[self.gridKwarg]

	@integGrid.setter
	def integGrid(self,value):
		self.optDict[self.gridKwarg] = value

	def getStrDictWithStruct(self, struct:"ucell obj"):
		strDict = self.getStrDictFromOptDict(self.optDict)
		geomDict = platoInp.getPlatoGeomDictFromUnitCell(struct)
		strDict.update(geomDict)
		return strDict






class CalcObj():
	def __init__(self, filePath, strDict, strDictWriteFunction, fileParser, runCommFunction, inpFileExt, outFileExt):
		self.filePath = os.path.splitext(filePath)[0]
		self.strDict = strDict
		self.inpFileExt = inpFileExt
		self.outFileExt = outFileExt
		self._strDictWriteFunction = strDictWriteFunction
		self._fileParser = fileParser
		self._runCommFunction = runCommFunction

	@classmethod
	def fromEnforcedKwargs(cls, **kwargs):
		kwargs = {k.lower():v for k,v in kwargs.items()}
		reqArgsOrder = ["filePath", "strDict", "strDictWriteFunction", "fileParser",
		                "runCommFunction", "inpFileExt", "outFileExt"]
		argsList = list()
		for x in reqArgsOrder:
			argsList.append( kwargs[x.lower()] )
	
		return cls(*argsList)

	def getRunComm(self):
		return self._runCommFunction(self.filePath + self.inpFileExt)

	def writeFile(self):
		writeFolder = os.path.split(self.filePath)[0]
		pathlib.Path(writeFolder).mkdir(exist_ok=True, parents=True)
		self._strDictWriteFunction(self.filePath + self.inpFileExt, self.strDict)

	def parseOutFile(self):
		return self._fileParser(self.filePath + self.outFileExt)

