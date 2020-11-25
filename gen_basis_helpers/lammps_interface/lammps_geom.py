
import plato_pylib.shared.ucell_class as uCellHelp

import collections
import itertools as it


class LammpsGeom():

	def __init__(self, unitCell, eleToTypeIdx=None, eleToCharge=None, eleToMass=None, geomToBondInfo=None, 
	             geomToAngleInfo=None, geomToMoleculeIDs=None, getDataDictFunct=None):
		""" Initializer
		
		Args:
			unitCell: (plato_pylib UnitCell object)
			eleToTypeIdx: (dict) Keys are eleKeys while values are type indices (integers to identify atom types)
			eleToCharge: (dict) Keys are eleKeys while values are atomic charges (not always needed, potential dependent)
			eleToMass: (dict) Keys are eleKeys while 
			geomToBondInfo: f(unitCell)-> list of bond info. Not needed for some atom_style values
			geomToAngleInfo: f(unitCell)->list of angle info. Not needed for some atom_style values
			geomToMoleculeIDs: f(unitCell)->list of ints. Each is a molecule ID for the corresponding atom
			getDataDictFunct: f(LammpsGeom)->dict where the dict contains headers/values for writing a data file
 
		"""
		self.unitCell = unitCell
		self.eleToTypeIdx = getEleToTypeIdxMapFromUnitCell(unitCell) if eleToTypeIdx is None else eleToTypeIdx
		self.eleToCharge = eleToCharge
		self.eleToMass = uCellHelp.getEleKeyToMassDictStandard() if eleToMass is None else eleToMass
		self.geomToBondInfo = geomToBondInfo
		self.geomToAngleInfo = geomToAngleInfo
		self.geomToMoleculeIDs = geomToMoleculeIDs
		self.getDataDictFunct = getDataDictFunct



class LammpsSimulationBox():

	def __init__(self, xRange, yRange, zRange, tiltFactors=None):
		""" Initializer
		
		Args:
			xRange: [xlo, xhi] defines x-boundaries of box
			yRange: [ylo, yhi] defines y-boundaries of box
			zRange: [zlo, zhi] defines z-boundaries of box
			tiltFactors: (len-3 iter) [xy,xz,yz] 
			OVERALL the lattice vectors that get defined from these values are:
			A = (xhi-xlo,0,0); B = (xy,yhi-ylo,0); C = (xz,yz,zhi-zlo)
				 
		"""
		self._eqTol = 1e-6
		self.xRange = list(xRange)
		self.yRange = list(yRange)
		self.zRange = list(zRange)
		self.tiltFactors = [0,0,0] if tiltFactors is None else tiltFactors

	@classmethod
	def fromUnitCell(cls, unitCell):
		""" Alternative initializer. Note this relies on UnitCell class using the same arbitrary choice of lattice vectors as LAMMPS
		
		Args:
			unitCell: (plato_pylib UnitCell object)
				 
		"""
		lattVects = unitCell.lattVects
		xRange = [0, lattVects[0][0]]
		yRange = [0, lattVects[1][1]]
		zRange = [0, lattVects[2][2]]
		xy = lattVects[1][0]
		xz = lattVects[2][0]
		yz = lattVects[2][1]
		tiltFactors = [xy, xz, yz]
		return cls(xRange, yRange, zRange, tiltFactors=tiltFactors)

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)
		numListAttrs = ["xRange", "yRange", "zRange", "tiltFactors"]
		for attr in numListAttrs:
			selfVals, otherVals = getattr(self,attr), getattr(other,attr)
			if len(selfVals) != len(otherVals):
				return False
			else:
				diffs = [abs(x-y) for x,y in zip(selfVals,otherVals)]
				if any([x>eqTol for x in diffs]):
					return False
		return True


class GetDataDictFromLammpsGeomAtomStyleFull():
	""" Callable class. f(LammpsGeom instance)->OrderedDict is the interface, where the OrderedDict contains keys/values for the data file

	"""

	def __init__(self, numbFmtCoords="{:.5f}", numbFmtCharge="{:.4f}", numbFmtMasses="{:.4f}", modTiltFactors=None):
		""" Initializer
		
		Args:
			numbFmtCoords: (str) format string for atomic co-ordinates
			numbFmtCharge: (str) format string for atomic charges
			numbFmtMasses: (str) format string for atomic masses
			modTiltFactors: f([xy,xz,yz])->[xy,xz,yz] This is useful to deal with small float errors in xy for hexagonal cells. A perfect hexagonal cell has xy=0.5*xhi which is the MAXIMUM value allowed. Thus if a float error is in the wrong direction it can lead to an error. Setting modTiltFactors can deal with this by (for example) multiplying xy by 0.999 or similar
 
		"""
		self.numbFmtCoords = numbFmtCoords
		self.numbFmtCharge = numbFmtCharge
		self.numbFmtMasses = numbFmtMasses
		self.modTiltFactors = modTiltFactors

	def getDataDict(self, lammpsGeomInstance):
		mapFuncts = createMapGeomInfoToDataDictObjFromAtomStyle("full", modTiltFactors=self.modTiltFactors)
		outDict = collections.OrderedDict()

		#Get the various info dict
		bondInfo = lammpsGeomInstance.geomToBondInfo(lammpsGeomInstance.unitCell)
		angleInfo = lammpsGeomInstance.geomToAngleInfo(lammpsGeomInstance.unitCell)
		moleculeIDs = lammpsGeomInstance.geomToMoleculeIDs(lammpsGeomInstance.unitCell)
		massDict = dict()
		for key in lammpsGeomInstance.eleToTypeIdx:
			currType = lammpsGeomInstance.eleToTypeIdx[key]
			currMass = lammpsGeomInstance.eleToMass[key]
			massDict[currType] = currMass


		currKwargDict = {"numbFmt":self.numbFmtCoords, "bondInfo":bondInfo, "angleInfo":angleInfo}
		outDict["LAMMPS Atom File"] = mapFuncts.getHeaderStrFromUnitCellAndKwargs(lammpsGeomInstance.unitCell, **currKwargDict)
		outDict["Masses"] = mapFuncts.getMassesDictStrFromMassDict(massDict, numbFmt=self.numbFmtMasses)

		currKwargDict = {"eleToTypeIdx":lammpsGeomInstance.eleToTypeIdx, "eleToChargeIdx":lammpsGeomInstance.eleToCharge,
		                 "moleculeIds":moleculeIDs, "coordNumbFmt":self.numbFmtCoords, "chargeNumbFmt":self.numbFmtCharge}
		outDict["Atoms"] = mapFuncts.getAtomsStrFromUnitCell(lammpsGeomInstance.unitCell, **currKwargDict)
		outDict["Bonds"] = mapFuncts.getConnectivityStr(bondInfo)
		outDict["Angles"] = mapFuncts.getConnectivityStr(angleInfo)

		return outDict



	def __call__(self, lammpsGeomInstance):
		return self.getDataDict(lammpsGeomInstance)



class MapGeomInfoToDataDictBase():
	
	def getMassesDictStrFromMassDict(self, massDict, numbFmt="{:.4f}"):
		raise NotImplementedError("")


def createMapGeomInfoToDataDictObjFromAtomStyle(atomStyle, modTiltFactors=None):
	""" Description of function
	
	Args:
		atomType: (str) The value of the atom_style keyword. The value of this keyword affects the format of the output file (e.g. if type was atom then bonds dont need defining)
			 
	Returns
		outMapper: (MapGeomInfoToDataDictBase object). Contains various functions for mapping geometry to data file strings for a given atomType
		modTiltFactors: f([xy,xz,yz])->[xy,xz,yz] This is useful to deal with small float errors in xy for hexagonal cells. A perfect hexagonal cell has xy=0.5*xhi which is the MAXIMUM value allowed. Thus if a float error is in the wrong direction it can lead to an error. Setting modTiltFactors can deal with this by (for example) multiplying xy by 0.999 or similar
 
	Raises:
		 KeyError: If atomStyle is not a valid keyword
	"""
	if atomStyle.lower()=="full":
		return _MapGeomInfoToDataDict_atomTypeFull(modTiltFactors=modTiltFactors)
	else:
		raise KeyError(atomStyle)


class _MapGeomInfoToDataDict_atomTypeFull(MapGeomInfoToDataDictBase):

	def __init__(self,modTiltFactors=None):
		self.modTiltFactors = modTiltFactors

	def getAtomsStrFromUnitCell(self, unitCell, eleToTypeIdx=None, eleToChargeIdx=None, moleculeIds=None,
	                            coordNumbFmt="{:.5f}", chargeNumbFmt="{:.4f}"):

		if eleToTypeIdx is None:
			raise ValueError(eleToTypeIdx) 
		elif eleToChargeIdx is None:
			raise ValueError(eleToChargeIdx)
		elif moleculeIds is None:
			raise ValueError(moleculeIds)
	
		outFmt = "\t{} {} {} " + chargeNumbFmt + " " + " ".join([coordNumbFmt for x in range(3)]) + "\n"
		cartCoords = unitCell.cartCoords
		outStr = ""
		for idx,(coord,molId) in enumerate(it.zip_longest(cartCoords,moleculeIds), start=1):
			currEle = coord[-1]
			currType, currCharge = eleToTypeIdx[currEle], eleToChargeIdx[currEle]
			currLine = outFmt.format(idx, molId, currType, currCharge, *coord[:3])
			outStr += currLine
		return outStr

	def getConnectivityStr(self, conInfo):
		""" Get string for either bond or angles section.
		
		Args:
			conInfo: (iter of iters) Each element is one line, represented by a list of ints (4 values for bondInfo, 5 values for angleInfo)
				 
		Returns
			conStr: The string used to define either bond/angle sections in the output data file
	 
		"""
		outStr = ""
		for info in conInfo:
			outStr += "\t" + "\t".join([str(x) for x in info]) + "\n"
		return outStr

	#TODO: Can likely be factored into base class
	def getHeaderStrFromUnitCellAndKwargs(self, unitCell, numbFmt="{:.4f}", bondInfo=None, angleInfo=None):
		cellDimsStr = self._getHeaderCellDimsPartFromUnitCell(unitCell, numbFmt=numbFmt)
		nBonds = 0 if bondInfo is None else len(bondInfo)
		nAngles = 0 if angleInfo is None else len(angleInfo)
		nDihedrals = 0
		nImproper = 0
		nAtomTypes =len( getEleToTypeIdxMapFromUnitCell(unitCell).keys() )
		nBondTypes = max([x[1] for x in bondInfo]) if bondInfo is not None else 0
		nAngleTypes = max([x[1] for x in angleInfo]) if angleInfo is not None else 0
		outStr = ""
		outStr += "\t{} atoms\n".format(len(unitCell.cartCoords))
		outStr += "\t{} bonds\n".format(nBonds)
		outStr += "\t{} angles\n".format(nAngles)
		outStr += "\t{} dihedrals\n".format(nDihedrals)
		outStr += "\t{} impropers\n".format(nImproper)
		outStr += "\t{} atom types\n".format(nAtomTypes)
		outStr += "\t{} bond types\n".format(nBondTypes)
		outStr += "\t{} angle types\n".format(nAngleTypes)
		outStr += cellDimsStr
		return outStr

	def _getHeaderCellDimsPartFromUnitCell(self, unitCell, numbFmt="{:.4f}"):
		simBox = LammpsSimulationBox.fromUnitCell(unitCell)
		outStr = self._getHeaderCellDimsPartFromLammpsSimBox(simBox, numbFmt=numbFmt)
		return outStr

	def _getHeaderCellDimsPartFromLammpsSimBox(self, simBox, numbFmt="{:.4f}"):
		rangeFmtPrefix = "\t" + numbFmt + "\t" + numbFmt +"\t"
		tiltFactorFmt = "\t" + numbFmt + "\t" + numbFmt + "\t" + numbFmt +"\txy xz yz"
		outStr = ""
		outStr += rangeFmtPrefix.format(*simBox.xRange) + "xlo xhi\n"
		outStr += rangeFmtPrefix.format(*simBox.yRange) + "ylo yhi\n"
		outStr += rangeFmtPrefix.format(*simBox.zRange) + "zlo zhi\n"
		tiltFactors = simBox.tiltFactors if self.modTiltFactors is None else self.modTiltFactors(simBox.tiltFactors)
		outStr += tiltFactorFmt.format(*tiltFactors)
		return outStr



	#TODO: Can likely factor into the base class
	def getMassesDictStrFromMassDict(self, massDict, numbFmt="{:.4f}"):
		""" Get the body string for the "Masses" entry from a dict of masses
		
		Args:
			massDict: (dict) Keys are atomType indices (1,2,3..) while values are masses to assign to those types
			numbFmt: (str) Define how we format the mass into a string
 
		Returns
			massStr: (str) The string used to define masses in the output data file. Note these values will always be in the order of atom indices
	 
		"""
		lineFmtStr = "\t{}\t"+numbFmt+"\n"
		outStr = ""
		for key,val in massDict.items():
			outStr += lineFmtStr.format(key,val)
		return outStr


def getEleToTypeIdxMapFromUnitCell(unitCell):
	""" Gets a dictionary mapping element types (e.g. "Mg") to type indices (e.g. 1,2) used in lammps
	
	Args:
		unitCell: (plato_pylib UnitCell)
			 
	Returns
		 outDict: (dict) e.g. {"Mg":1, "O":2, "H":3}. Always starts at 1 and the order depends on the order in the unitCell (if Mg is the first atom in unitCell.cartCoords then it will have idx=1)
 
	"""
	cartCoords = unitCell.cartCoords
	eleTypeList = list()
	outDict = dict()

	currIdx = 1
	for coord in cartCoords:
		currEle = coord[-1]
		if currEle not in eleTypeList:
			eleTypeList.append(currEle)
			outDict[currEle] = currIdx
			currIdx += 1

	return outDict

#Taken from https://gist.github.com/lukasrichters14/c862644d4cbcf2d67252a484b7c6049c
def _getEleKeyToMassDictStandard():
	elements_dict = {'H' : 1.008,'HE' : 4.003, 'LI' : 6.941, 'BE' : 9.012,\
	             'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
	             'F' : 18.998, 'NE' : 20.180, 'NA' : 22.990, 'MG' : 24.305,\
	             'AL' : 26.982, 'SI' : 28.086, 'P' : 30.974, 'S' : 32.066,\
	             'CL' : 35.453, 'AR' : 39.948, 'K' : 39.098, 'CA' : 40.078,\
	             'SC' : 44.956, 'TI' : 47.867, 'V' : 50.942, 'CR' : 51.996,\
	             'MN' : 54.938, 'FE' : 55.845, 'CO' : 58.933, 'NI' : 58.693,\
	             'CU' : 63.546, 'ZN' : 65.38, 'GA' : 69.723, 'GE' : 72.631,\
	             'AS' : 74.922, 'SE' : 78.971, 'BR' : 79.904, 'KR' : 84.798,\
	             'RB' : 84.468, 'SR' : 87.62, 'Y' : 88.906, 'ZR' : 91.224,\
	             'NB' : 92.906, 'MO' : 95.95, 'TC' : 98.907, 'RU' : 101.07,\
	             'RH' : 102.906, 'PD' : 106.42, 'AG' : 107.868, 'CD' : 112.414,\
	             'IN' : 114.818, 'SN' : 118.711, 'SB' : 121.760, 'TE' : 126.7,\
	             'I' : 126.904, 'XE' : 131.294, 'CS' : 132.905, 'BA' : 137.328,\
	             'LA' : 138.905, 'CE' : 140.116, 'PR' : 140.908, 'ND' : 144.243,\
	             'PM' : 144.913, 'SM' : 150.36, 'EU' : 151.964, 'GD' : 157.25,\
	             'TB' : 158.925, 'DY': 162.500, 'HO' : 164.930, 'ER' : 167.259,\
	             'TM' : 168.934, 'YB' : 173.055, 'LU' : 174.967, 'HF' : 178.49,\
	             'TA' : 180.948, 'W' : 183.84, 'RE' : 186.207, 'OS' : 190.23,\
	             'IR' : 192.217, 'PT' : 195.085, 'AU' : 196.967, 'HG' : 200.592,\
	             'TL' : 204.383, 'PB' : 207.2, 'BI' : 208.980, 'PO' : 208.982,\
	             'AT' : 209.987, 'RN' : 222.081, 'FR' : 223.020, 'RA' : 226.025,\
	             'AC' : 227.028, 'TH' : 232.038, 'PA' : 231.036, 'U' : 238.029,\
	             'NP' : 237, 'PU' : 244, 'AM' : 243, 'CM' : 247, 'BK' : 247,\
	             'CT' : 251, 'ES' : 252, 'FM' : 257, 'MD' : 258, 'NO' : 259,\
	             'LR' : 262, 'RF' : 261, 'DB' : 262, 'SG' : 266, 'BH' : 264,\
	             'HS' : 269, 'MT' : 268, 'DS' : 271, 'RG' : 272, 'CN' : 285,\
	             'NH' : 284, 'FL' : 289, 'MC' : 288, 'LV' : 292, 'TS' : 294,\
	             'OG' : 294}
	outDict = dict()
	for key,val in elements_dict.items():
		outDict[key.capitalize()] = val
	return outDict




