

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




