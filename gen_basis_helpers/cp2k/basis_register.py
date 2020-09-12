

""" Module holds the objects needed to map user-strings to CP2K basis-set objects """

import os
import plato_pylib.parseOther.parse_cp2k_basis as parseCP2KBasis

import gen_basis_helpers.cp2k.private.get_default_basis_strs as getDefaultBasisStrs


#_METHOD_STRS_TO_BASIS_SETS =
_BASIS_STRS_TO_CP2K_BASIS_CREATORS = getDefaultBasisStrs.defaultBasisStrsToCreators
_BASIS_STRS = set([key for key in _BASIS_STRS_TO_CP2K_BASIS_CREATORS.keys()])



def getRegisteredCP2KObjCreatorStrings():
	return set(_BASIS_STRS)


def decoRegisterCP2KBasisCreatorToMethodStr(basisStr, overwrite=False):
	def registerFunct(funct):
		registerCP2KBasisCreatorToMethodStr(funct,basisStr,overwrite=overwrite)
		return funct
	return registerFunct


#Duplicated from method_register essentially
def registerCP2KBasisCreatorToMethodStr(creator, basisStr, overwrite=False):
	""" Registers function "creator" to the string basisStr. This means that createCP2KBasisObjFromEleAndBasisStr(eleStr, methodStr) will return the object defined by the creator function
	
	Args:
		creator: (function, returns ) Function that takes no arguments and returns 
		basisStr: (str) Key which you want to use to make this creation method accesible. Note this is case insensitive
		overwrite: (Bool, optional) Whether to overwrite methodStr if its already been registered. Default is False.

	Raises:
		AssertionError: If methodStr has already been registered and overwrite=False 

	"""

	if overwrite is False:
		assert basisStr.lower() not in getRegisteredCP2KObjCreatorStrings(), "{} must not already be registerd, unless overwrite is set to True".format(basisStr)
	_BASIS_STRS.add(basisStr.lower())
	_BASIS_STRS_TO_CP2K_BASIS_CREATORS[basisStr.lower()] = creator




def createCP2KBasisObjsFromStrDict(strDict):
	""" Creates an iter of CP2KBasis objects for a dictionary input 
	
	Args:
		strDict:(dict) keys are element strs (e.g. "H","Mg") while values are keys to specific basis set/pseudopot. combination
			
	Returns
		basisObjs: (iter of CP2KBasisObjStandard) Each entry represents the basis-set information for one element
	
	Raises:
		KeyError: If any values in strDict arent valid basis set identifiers. See getRegisteredCP2KObjCreatorStrings for available options
	"""
	outList = list()
	for key in strDict.keys():
		currBasisObj = createCP2KBasisObjFromEleAndBasisStr(key, strDict[key])
		outList.append( currBasisObj )
	return outList


def createCP2KBasisObjFromEleAndBasisStr(eleStr, basisStr):
	""" Creates a CP2KBasis object when given an element symbol and string represnting basis set/pseudopot to use. See getRegisteredCP2KObjCreatorStrings for available options
	
	Args:
		eleStr: (str, case-insensitive) The element used
		basisStr: (str, case-insensitive) The key to a specific basis set/pseudopot. combination
			
	Returns
		outObj: (CP2KBasisObjStandard) Object representation of information needed to add basis-set information to a CP2K file for one element
	
	Raises:
		KeyError: If basisStr is not present
	"""
	outObj = _BASIS_STRS_TO_CP2K_BASIS_CREATORS[basisStr.lower()]()
	outObj.element = eleStr.capitalize()

	return outObj

def createCP2KBasisSetInPlatoFmt(eleStr, basisStr, basisFolderPath):
	""" Gets a CP2K basis set in a format often used with plato (including reversing the normalisation of coefficients)
	
	Args:
		eleStr: (str, case-insensitive) The element used
		basisStr: (str, case-insensitive) The key to a specific basis set/pseudopot. combination
		basisFilePath: (str) Path to where the CP2K basis files are
	
	Returns
		gauPolyExpansions: (iter of GauPolyBasis objects) Each represents ONE basis function in the standard format for plato. Angular momentum stored on .label
	
	"""
	origBasisObj = createCP2KBasisObjFromEleAndBasisStr(eleStr, basisStr)
	origBasisPath = os.path.join( basisFolderPath, origBasisObj.basisFile )
	parsedFile = parseCP2KBasis.parseCP2KBasisFile(origBasisPath)
	cp2kFormatBasis = parsedFile.getUniqueBasisSet(eleStr, origBasisObj.basis)
	platoFmtBasis = parseCP2KBasis.getGauPolyBasisFunctionsFromCP2KBasisSet(cp2kFormatBasis)
	return platoFmtBasis

