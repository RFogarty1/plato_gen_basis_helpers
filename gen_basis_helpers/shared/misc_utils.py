
import functools
import types

#Use as a function really (e.g. with fragile(open("file.txt","rt")) as f:)
#call raise fragile.Break to escape out of the context manager
class fragile(object):
    class Break(Exception):
      """Break out of the with statement"""

    def __init__(self, value):
        self.value = value

    def __enter__(self):
        return self.value.__enter__()

    def __exit__(self, etype, value, traceback):
        error = self.value.__exit__(etype, value, traceback)
        if etype == self.Break:
            return True
        return error


#Class decorator to apply a composite search method. We need a composite and non-composite version of this
def getObjectsWithComponentsInstanceWrapper(isComposite=None):
	""" Attach a function getObjectsWithComponents to an instance. Function searches through object
	
	Args:
		inpCls: Input instnace. This must have a .label attribute which returns a BaseLabel object (see shared). If isComposite the instance also needs to have a *.objs attribute which stores the leaf objects.
		 
		isComposite (bool): True means decorate a composite class, False means decorate a leaf class
			
	Returns
		The instance will have the "getObjectsWithComponents(self, components, caseSensitive=True)" function attached. Components
		is an iter of strings used to identify the object. e.g ["hcp","methodA"]. The function returns a list of objects which have
		all the requested components as part of their labels 
			

	"""

	if isComposite is None:
		raise ValueError("isComposite = {} is invalid. It must be True or False".format(isComposite))

	def clsDeco(inpInit):
		@functools.wraps(inpInit)
		def wrappedInit(instance,*args,**kwargs):
			inpInit(instance, *args, **kwargs)
			if isComposite:
				instance.getObjectsWithComponents = types.MethodType( _getObjectsWithComponentsForComposite, instance )
			else:
				instance.getObjectsWithComponents = types.MethodType( _getObjectsWithComponentsForNonComposites, instance )

		return wrappedInit


	return clsDeco


def _getObjectsWithComponentsForComposite(self, components, caseSensitive=True):
	""" Returns objects which match the components
	
	Args:
		components (str list): Strings which will be matched against BaseLabel (or inherited cls) components. Any object with 
	all of these components present will be returned 
		caseSensitive(bool): Whether we need to match the case of each component or not 
	Returns
		objList (iter): List of output objects where labels match components requested. 
 
	"""
	outObjs = list()
	for x in self.objs:
		outObjs.extend( x.getObjectsWithComponents(components, caseSensitive=caseSensitive) )
	return outObjs


def _getObjectsWithComponentsForNonComposites(self, components, caseSensitive=True):
	""" Returns analyser objects which match the components
	
	Args:
		components (str list): Strings which will be matched against self.label.components (or inherited cls). 
		caseSensitive(bool): Whether we need to match the case of each component or not 
	Returns
		objList (iter): List of objects. Empty if object components dont match, else length 1 with reference to self. The iter is required
		                since the same interface needs to be present for composite/non-composite objects
 
	"""

	inpComponents = list(components)
	assert len(self.label)==1
	allComps = list(self.label[0].components)
	if not caseSensitive:
		inpComponents = [x.lower() for x in inpComponents]
		allComps = [x.lower() for x in allComps]


	objConsistentWithComponents = True
	for x in inpComponents:
		if x not in allComps:
			objConsistentWithComponents = False

	if objConsistentWithComponents:
		return [self]
	else:
		return list()



 
