
import os

class ChDir(object):
    """
    Step into a directory temporarily. Propagates exceptions upwards (since __exit__ always returns None)
    """
    def __init__(self, startDir, workDir):
        self.old_dir = startDir
        self.new_dir = workDir
 
    def __enter__(self):
        os.chdir(self.new_dir)
 
    def __exit__(self, *args):
        os.chdir(self.old_dir)

