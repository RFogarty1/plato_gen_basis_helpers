

import gen_basis_helpers.castep.private.get_default_method_strs as get_default_method_strs


_METHOD_STRS_TO_PARAM_DICTS = get_default_method_strs.defaultMethodStrsToObjCreators
_METHOD_STRS = set([key for key in _METHOD_STRS_TO_PARAM_DICTS.keys()]) #Used to help track if we have duplicates

#Shamelessly duplicated everything from the CP2K version. TODO: Should factor out the registration functions from both really


def decoRegisterCreatorToMethodStr(methodStr, overwrite=False):
	def registerFunct(funct):
		registerCreatorToMethodStr(funct,methodStr, overwrite=overwrite)
		return funct
	return registerFunct


def registerCreatorToMethodStr(creator, methodStr, overwrite=False):
	""" Registers function "creator" to the string methodStr. This means that createCP2KObjFromMethodStr(methodStr) will return the object defined by the creator function
	
	Args:
		creator: (function, returns castep param dict) Function that takes no arguments and returns param dict
		methodStr: (str) Key which you want to use to make this creation method accesible. Note this is case insensitive
		overwrite: (Bool, optional) Whether to overwrite methodStr if its already been registered. Default is False.

	Raises:
		AssertionError: If methodStr has already been registered and overwrite=False 

	"""

	if overwrite is False:
		assert methodStr.lower() not in getRegisteredCreatorStrings(), "{} must not already be registerd, unless overwrite is set to True".format(methodStr)
	_METHOD_STRS.add(methodStr.lower())
	_METHOD_STRS_TO_PARAM_DICTS[methodStr.lower()] = creator


def getRegisteredCreatorStrings():
	return set(_METHOD_STRS)


def createParamDictFromMethodStr(methodStr):
	""" Creates a castep param dict corresponding to the input methodStr. See getRegisteredCreatorStrings() for available strings.
	
	Args:
		methodStr: (str, case-insensitive) The key to a specific *.param file representation
			
	Returns
		outObj: (param dict) Dictionary representation of a castep *.param file.
	
	Raises:
		KeyError: If methodStr is not present
	"""
	return _METHOD_STRS_TO_PARAM_DICTS[methodStr.lower()]()



