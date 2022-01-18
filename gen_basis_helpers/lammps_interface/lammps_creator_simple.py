
import collections

from . import lammps_calc_objs as calcObjHelp
from ..shared import method_objs as baseObjs



class LammpsCalcObjFactorySimple(baseObjs.CalcMethodFactoryBase):

	registeredKwargs = set(baseObjs.CalcMethodFactoryBase.registeredKwargs)

	registeredKwargs.add("units") #TODO: Set default value for units to real
	registeredKwargs.add("atomStyle")
	registeredKwargs.add("timeStep")
	registeredKwargs.add("thermoStyle")
	registeredKwargs.add("printThermoEveryNSteps")
	registeredKwargs.add("velocityObj")
	registeredKwargs.add("lammpsGeomObj")
	registeredKwargs.add("potentialObj")
	registeredKwargs.add("ensembleObj")
	registeredKwargs.add("nSteps")
	registeredKwargs.add("dumpOptions")
	registeredKwargs.add("boundaries")
	registeredKwargs.add("walls")
	registeredKwargs.add("chargeEqualisation")
	registeredKwargs.add("groups")
	registeredKwargs.add("velocityRescaling")

	#KEY FUNCTION
	def _createFromSelf(self):
		scriptFileDict = self._getScriptFileDict()
		dataFileDict = self._getDataFileDict()
		typeIdxToEle = {v:k for k,v in self.lammpsGeomObj.eleToTypeIdx.items()}
		outObj = calcObjHelp.LammpsCalcObjStandard(self.workFolder, self.fileName, dataFileDict, scriptFileDict, typeIdxToEle=typeIdxToEle)
		return outObj

	def _getDataFileDict(self):
		dataDict = self.lammpsGeomObj.getDataDictFunct(self.lammpsGeomObj)
		return dataDict

	def _getScriptFileDict(self):
		outOptions = calcObjHelp.ScriptFileOptionsStandard()
		outOptions.initOpts = self._getInitCommandDict()
		outOptions.forceFieldOpts = self.potentialObj.commandDict if self.potentialObj is not None else collections.OrderedDict()
		outOptions.fixSection = self._getFixCommandDict()
		outOptions.settingsSection = self._getSettingsCommandDict()
		outOptions.runSection = collections.OrderedDict([["run",str(self.nSteps)]]) if self.nSteps is not None else collections.OrderedDict()
		outOptions.outputSection = self.dumpOptions.commandDict if self.dumpOptions is not None else collections.OrderedDict()
		return outOptions.getOutputDict()

	def _getInitCommandDict(self):
		outComms = list()
		if self.units is not None:
			outComms.append( ["units", self.units] )

		if self.atomStyle is not None:
			outComms.append( ["atom_style", self.atomStyle] )

		#Needs defining BEFORE we define the box
		if self.boundaries is not None:
			outComms.append( ["boundary", " ".join([str(x) for x in self.boundaries]) ] )

		outComms.append(["read_data", self.fileName+".data"])

		return collections.OrderedDict(outComms)

	def _getFixCommandDict(self):
		outStr = ""
		currFixIdx = 1

		def _appendStdFixCommandForAttr(attr, outStr, currFixIdx):
			currVal = getattr(self,attr)
			if currVal is not None:
				if outStr!="":
					outStr += "\nfix "

				try:
					iter(currVal)
				except TypeError:
					outStr += str(currFixIdx) + " " + currVal.fixStr
					currFixIdx += 1
				else:
					for idx,currObj in enumerate(currVal):
						if idx==0:
							outStr += str(currFixIdx) + " " + currObj.fixStr
						else:
							outStr += "\nfix " + str(currFixIdx) + " " + currObj.fixStr
						currFixIdx += 1

			return outStr, currFixIdx

		outStr, currFixIdx = _appendStdFixCommandForAttr("ensembleObj", outStr, currFixIdx)
		outStr, currFixIdx = _appendStdFixCommandForAttr("chargeEqualisation", outStr, currFixIdx)
		outStr, currFixIdx = _appendStdFixCommandForAttr("velocityRescaling", outStr, currFixIdx)


		if self.walls is not None:
			for wall in self.walls:
				if outStr=="":
					outStr += str(currFixIdx) + " " + wall.fixStr
				else:
					outStr += "\nfix {} {}".format( str(currFixIdx), wall.fixStr )
				currFixIdx += 1

		if outStr == "":
			return collections.OrderedDict()
		
		return collections.OrderedDict([["fix",outStr]])
	
	def _getSettingsCommandDict(self):
		outDict = collections.OrderedDict()

		if self.groups is not None:
			outDict["group"] = self._getAllGroupsCommand()

		if self.velocityObj is not None:
			outDict["velocity"] = self.velocityObj.commandStr

		if self.timeStep is not None:
			outDict["timestep"] = "{:.4f}".format(self.timeStep)

		if self.thermoStyle is not None:
			outDict["thermo_style"] = self.thermoStyle

		if self.printThermoEveryNSteps is not None:
			outDict["thermo"] = str(self.printThermoEveryNSteps)

		return outDict

	def _getAllGroupsCommand(self):
		
		try:
			iter(self.groups)
		except TypeError:
			outStr = self.groups.groupStr
		else:
			outStr = ""
			for idx,obj in enumerate(self.groups):
				if idx==0:
					outStr += obj.groupStr
				else:
					outStr += "\ngroup " + obj.groupStr

		return outStr

