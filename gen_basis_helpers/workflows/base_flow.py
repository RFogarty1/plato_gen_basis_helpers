
import gen_basis_helpers.shared.misc_utils as misc

""" These workflows are not immediately compatable with those from the fitting code. I needed workflows that were composite-friendly, and those simply werent """

class BaseWorkflow():
	"""BaseWorkflow: An object used to encapsulate the calculation of some property (e.g. a defect energy)
	"""


	@property
	def preRunShellComms(self):
		""" List of string commands that will get run before run() is called. Higher-level functions can therefore 
		combine these for a set of WorkFlows, which should lead to more efficient parralelisation. Return None or
		empty list to not run any of these commands.
		
		"""
		return list()

	@property
	def namespaceAttrs(self):
		""" iter of Lists of attribute names for the properties this workspace calculates + places in a Namespace. See run()
		Note: Each leaf object has 1 list of these, but we always return an iter of lists to make it easier to use composite workflows
		
		"""
		raise NotImplementedError()

	@property
	def workFolder(self):
		""" List of workFolders used. Returning an empty list is allowed. This may (or may not) be used to check our workFlows dont overwrite each others work. (Note: This is the Base, abstract property docstring) """
		raise NotImplementedError()

	def run(self):
		""" Runs the workflow and populates the output attr with an iter of Namespaces
		containing calculated properties/values as fields/values. Field names should match those in namespaceAttrs
				
		Note: This should generally be a pretty fast analysis part, preRunShellComms are generally used to handle expensive parts (e.g. running a qunatum chemistry calculation).
		"""
		raise NotImplementedError()



class BaseLabelledWorkflow(BaseWorkflow):
	"""BaseLabelledWorkflow: An object used to encapsulate the calculation of some property (e.g. a defect energy). Contains a BaseLabel object to allow better implementation of Composite patterns

	"""

	@property
	def label(self):
		""" BaseLabel object """
		raise NotImplementedError("")



class StandardLabelledWorkflowComposite(BaseLabelledWorkflow):

	namespaceAttrs = misc.StandardComponentDescriptor("namespaceAttrs")
	label = misc.StandardComponentDescriptor("label")
	workFolder = misc.StandardComponentDescriptor("workFolder")
	preRunShellComms = misc.StandardComponentDescriptor("preRunShellComms")
	def __init__(self, objs):
		self.objs = objs




