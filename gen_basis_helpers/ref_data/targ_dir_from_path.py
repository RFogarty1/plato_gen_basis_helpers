
import os

class GetSpecialDirFromSubDir():
	""" Callable class; call redirects to self.getDir; see the docstring for that function for the callable interface. The original use case for this is to get a working directory for jobs based on the path of a jupyter notebook which created the jobs

	"""

	def __init__(self, baseDir, targetSubDir):
		""" Initializer
		
		Args:
			baseDir: (str, path) Path to a base directory; input path to callable interface MUST be downstream from this path
			targetSubDir: (str, path). os.path.join(baseDir, targetSubDir) should get the target directory
				 
		"""
		self.baseDir = os.path.abspath(baseDir)
		self.targetSubDir = targetSubDir


	def getDir(self, inpDir, extension=None):
		""" Gets target directory from input directory based on instance parameters 
		
		Args:
			inpDir: (str, path) Directory containing the jupyter notebook
			extension: (str,path-like, Optional) Added to inpDir path. So could get same effect from just appending to inpDir if required
				 
		Raises:
			AssertionError: If inpDir is not a sub-path of self.baseDir 
		"""
		if os.path.commonprefix([inpDir, self.baseDir]) != self.baseDir:
			raise AssertionError("inpDir {} must be a sub-directory of {}".format(inpDir,self.baseDir))
		extraBit = os.path.relpath(inpDir, self.baseDir)

		if extension is None:
			outPath = os.path.join(self.baseDir, self.targetSubDir, extraBit)
		else:
			outPath = os.path.join(self.baseDir, self.targetSubDir, extraBit, extension)

		return outPath

	def __call__(self, *args, **kwargs):
		return self.getDir(*args, **kwargs)

