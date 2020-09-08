
from . import core as coreHelp


""" Code used to let user follow the progress of optimisiation to some extent """

def makeObjFunctPrintObjValEveryNIters(printEveryNIters, objFunct):
	""" Causes printing of iter number and objective function value every n-times the objective function is evaluated
	
	Args:
		printEveryNIters: (int) How many iters to pass between each print. E.g. 1 means print every time the objective function is evaluated, 2 every other time etc.
		objFunct: (ObjFunctCalculatorStandard object INSTANCE) 
			 
	"""
	printObserver = PrintIterNumberVsObjVal(printEveryNIters)
	objFunct.addObjValObserver(printObserver)

class PrintIterNumberVsObjVal(coreHelp.ObjFunctObserver):
	""" Class used for printing optimisation progress to screen at times based on how often the objective functio nhas been called

	"""
	def __init__(self, printEveryNIters):
		""" Initializer
		
		Args:
			printEveryNIters: (int) How many iters to pass between each print. E.g. 1 means print every time the objective function is evaluated, 2 every other time etc.
				 
		"""
		self.iterNumb = 0
		self.printEveryNIters = printEveryNIters

	def updateObjVal(self, objVal):
		self.iterNumb += 1
		if (self.iterNumb%self.printEveryNIters)==0:
			self._printVals(self.iterNumb, objVal)

	def _printVals(self, iterNumb, objVal):
		print("Iter Number: {} objVal: {}".format(iterNumb, objVal))
