

""" Module holds the objects needed to map user-strings to methods """

import gen_basis_helpers.cp2k.private.get_default_method_strs as get_default_method_strs

_METHOD_STRS_TO_CP2K_OBJ_CREATORS = get_default_method_strs.defaultMethodStrsToObjCreators
_METHOD_STRS = set([key for key in _METHOD_STRS_TO_CP2K_OBJ_CREATORS.keys()]) #Used to help track if we have duplicates

def decoRegisterCP2KObjCreatorToMethodStr(methodStr, overwrite=False):
	def registerFunct(funct):
		registerCP2KObjCreatorToMethodStr(funct,methodStr, overwrite=overwrite)
		return funct
	return registerFunct


def registerCP2KObjCreatorToMethodStr(creator, methodStr, overwrite=False):
	""" Registers function "creator" to the string methodStr. This means that createCP2KObjFromMethodStr(methodStr) will return the object defined by the creator function
	
	Args:
		creator: (function, returns pycp2k CP2K object) Function that takes no arguments and returns cp2k object
		methodStr: (str) Key which you want to use to make this creation method accesible. Note this is case insensitive
		overwrite: (Bool, optional) Whether to overwrite methodStr if its already been registered. Default is False.

	Raises:
		AssertionError: If methodStr has already been registered and overwrite=False 

	"""

	if overwrite is False:
		assert methodStr.lower() not in getRegisteredCP2KObjCreatorStrings(), "{} must not already be registerd, unless overwrite is set to True".format(methodStr)
	_METHOD_STRS.add(methodStr.lower())
	_METHOD_STRS_TO_CP2K_OBJ_CREATORS[methodStr.lower()] = creator


def getRegisteredCP2KObjCreatorStrings():
	return set(_METHOD_STRS)

def createCP2KObjFromMethodStr(methodStr):
	""" Creates a pycp2k CP2K object corresponding to the input methodStr. See getRegisteredCP2KObjCreatorStrings() for available strings.
	
	Args:
		methodStr: (str, case-insensitive) The key to a specific cp2k obj
			
	Returns
		outObj: (pycp2k CP2K object) Object representation of a CP2K file.
	
	Raises:
		KeyError: If methodStr is not present
	"""
	return _METHOD_STRS_TO_CP2K_OBJ_CREATORS[methodStr.lower()]()


