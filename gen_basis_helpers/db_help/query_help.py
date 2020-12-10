
import os
from ..analyse_md import traj_core as trajHelp

def getSingleRecordFromCollectionFromSimpleQuery(collection, queryDict):
	""" Convenience function for returning a single record from a Mongodb collection when given a query dict
	
	Args:
		collection: (Collection object, pymongo package)
		queryDict: (dict) Passed directly as first argument to db.collection.find() function

	Returns
		record: (Dict)
	
	Raises:
		AssertionError: If query dict leads to either one or none returned
	"""
	relDocs = [x for x in collection.find(queryDict)]
	nDocs = len(relDocs)
	assert nDocs==1, "nDocs should equal 1 but is {} for queryDict {}".format(nDocs, queryDict)
	return relDocs[0]



def getTrajectoryFromRecord(record, folderKey="md_traj_folder",fileKey="md_traj_file"):
	""" Gets an MD trajectory from a database record in standard format
	
	Args:
		record: (dict) Record from mongodb
		folderKey: (str) The key in record which contains a folder
		fileKey: (str) The key in record which contains the name of the trajectory file
 
	Returns
		 outTraj: (TrajectoryInMemory object) This object contains all trajectory steps, stored in memory
 
	"""
	folder = record[folderKey]
	fileName = record[fileKey]
	fullPath = os.path.join(folder,fileName)
	outObj = trajHelp.readTrajObjFromFileToTrajectoryInMemory(fullPath)
	return outObj


def getFinalTrajStepFromRecord(record, folderKey="md_traj_folder",fileKey="md_traj_file"):
	""" Gets the final recorded step of an MD trajectory from a database record
	
	Args:
		record: (dict) Record from mongodb
		folderKey: (str) The key in record which contains a folder
		fileKey: (str) The key in record which contains the name of the trajectory file
			 
	Returns
		 outTraj: (TrajStepBase object) This object contains only the final trajectory step for a file
 
	Raises:
		 Errors
	"""
	folder = record[folderKey]
	fileName = record[fileKey]
	fullPath = os.path.join(folder,fileName)
	outObj = trajHelp.readLastTrajStepFromFile(fullPath)
	return outObj
