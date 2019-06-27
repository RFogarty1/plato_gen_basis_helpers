import os

from method_objs import CalcMethod
import plato_pylib.parseOther.parse_cp2k_files as parseCP2K



class CP2KCalcObj(CalcMethod):
	def __init__(self, pycp2kObj, basePath=None):
		self.cp2kObj = pycp2kObj
		if basePath is None:
			self.basePath = os.path.abspath( os.path.join( os.getcwd(), "cp2k_file" ) )
		else:
			self.basePath = os.path.abspath( os.path.splitext(basePath)[0] )
		
	
	def writeFile(self):
		self.cp2kObj.project_name = os.path.split(self.basePath)[1]
		self.cp2kObj.working_directory = os.path.split(self.basePath)[0]
		self.cp2kObj.write_input_file()
		
	@property
	def outFilePath(self):
		return self.basePath + ".cpout"
	
	@property
	def nCores(self):
		return 1 #I only implement the serial version. Error should occur if trying to set this
	
	@property
	def runComm(self):
		inpPath = self.basePath + ".inp"
		inpFolder = os.path.abspath(os.path.split(inpPath)[0] )
		inpFName = os.path.split(inpPath)[1]
		outFName = os.path.split(self.outFilePath)[1]
		commFmt = "cd {};cp2k.sopt {}>{}"
		return commFmt.format(inpFolder, inpFName, outFName)
	
	@property
	def parsedFile(self):
		parsedDict = parseCP2K.parseCpout(self.outFilePath)
		return SimpleNameSpace(**parsedDict)





#Optional descriptors that can be added

def addInpPathDescriptorToCP2KCalcObjCLASS(inpCls):
	attrName = "inpPath"
	setattr(inpCls, attrName, InpPathDescriptorCP2K(attrName))


#Defining an inpPath property for CP2KCalcObj
class InpPathDescriptorCP2K(object): 
	def __init__(self, attrName):
		self.name = attrName
	def __get__(self, instance, owner):
		return instance.basePath + ".inp"
  
	def __set__(self, instance, value):
		raise NotImplementedError("Cannot set attribute {}".format(self.name))

