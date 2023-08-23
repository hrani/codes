import sys,os
import json
import jsonschema
import argparse
import hillTau
from simError import SimError
import moose

def file_choices(choices,fname,iotype):
	ext = (os.path.splitext(fname)[1][1:]).lower()
	if iotype == "imagetype":
		if fname not in choices:
			parser.error("Requires output filetype {}".format(choices))
	# elif iotype == "outputfile":
	# 	if ext not in choices:
	# 		parser.error("Requires output filetype {}".format(choices))
	else:
		if ext not in ["json",'xml','g']:
			parser.error("Requires HillTau file in JSON format ")
			
	return fname
def dict_raise_on_duplicates(ordered_pairs):
	"""Reject duplicate keys."""
	d = {}
	for k, v in ordered_pairs:
		if k in d:
		   raise ValueError("duplicate key: %r" % (k,))
		else:
		   d[k] = v
	return d
def species_missing_model():
	''' Specie defined in map file but doesn't exist in model '''

def dupe_checking_hook(pairs):
	result = dict()
	for key,val in pairs:
		try:
			if key in result and set(val) != set(result[key]):
					raise SimError(" \'{}\' key has being assigned to 2 different values {} and {} with in the map file which is not allowed".format(key,result[key],val))
			if key in result and val == result[key]:
				raise SimError(" {} {} has definded multiple times in the map file ".format(key,val))
			else:
				if key.lower() not in ["filetype","version"]:  
					result[key] = val
		except SimError as duplicatekeyvalue:
			print("{}".format(duplicatekeyvalue))
			
	return result
def checkmapfile(mapfile):
  
	try:
		'''	    
	    	1. Here map file is looked 
				- if duplicate key,value pair exist its force to remove (to reduce loop time)
				- In map file for a given key,value pair where for given key 2 different value exist then need to be corrected as this is not allowed
		'''
		if not os.path.isfile(mapfile):
			raise SimError( "Map file name {} not found".format( args.map ) )
		else:
			fileName, file_extension = os.path.splitext(mapfile)
			if file_extension == '.json':
				with open(mapfile) as json_file:
					decoder = json.JSONDecoder(object_pairs_hook=dupe_checking_hook)
					maplist = decoder.decode(json_file.read())
					return maplist
			else:
				print("map file should be in json format")
				
	except SimError as mapmsg:
		print("{}".format(mapmsg))

def populate_modelobject(modelfile):
	try:
		'''
		 	2. Here based on the Model file type modelitems are extracted
		'''
		if not os.path.isfile(modelfile):
			raise SimError( "Model file name {} not found".format( args.model ) )
		fileName, file_extension = os.path.splitext(modelfile)
		if file_extension == '.g':
			moose.loadModel(modelfile,'/model')
			return moose.element('/model')
		elif file_extension == '.xml':
			moose.readSBML(modelfile,'/model')
			return moose.element('/model')
		elif file_extension == ".json":
			jsonDict = hillTau.loadHillTau(modelfile)
			htmodel = hillTau.parseModel(jsonDict)
			modelitems=htmodel.grpInfo
			for k,v in htmodel.molInfo.items():
				modelitems.append(k)
			for k,v in htmodel.reacInfo.items():
				modelitems.append(k)
			for k,v in htmodel.eqnInfo.items():
				modelitems.append(k)
			for k,v in htmodel.namedConsts.items():
				modelitems.append(k)
			return modelitems
	except SimError as msg:
		print( "{}".format(msg ))
		return
def main():

	parser = argparse.ArgumentParser( description = 'This program loads genesis/xml/json models into system and map file are validated against the model.\n')
	parser.add_argument('model',type=lambda s:file_choices((".g","xml","json"),s,"input"),help='Required: filename of model, in Genesis,XML,JSON format.')
	parser.add_argument( '-map', '--map', type = lambda s:file_choices(("json"),s,"input"), required="True",help='Optional: mapping file from json names to sim-specific strings. JSON format.')
	args = parser.parse_args()
	
	'''
		 Rule in map file 
			Allowed : we can have multiple keys for a given model object which is valid 
		   			  (b'cos different experiment might have different name which belongs to same model object' )
			Not Allowed: we can't assign for the given key differnt model object
	'''
	mapmodelObjectlist = checkmapfile(args.map)
	modelObjectlist = populate_modelobject(args.model)
	fileName, file_extension = os.path.splitext(args.model)
	itemsnotfound = {}
	multipleitems = {}
	if file_extension == '.json':
		if mapmodelObjectlist and modelObjectlist: 
			for k,v in mapmodelObjectlist.items():
				for i in v:
					if i not in modelObjectlist:
						itemsnotfound[i]={k:v}

	elif file_extension in['.g','.xml']:
		if mapmodelObjectlist and modelObjectlist: 
			for k,v in mapmodelObjectlist.items():
				for vv in v:
					i = []
					if vv.find('/') == -1:
						i = moose.wildcardFind(moose.element(modelObjectlist).path+"/##[FIELD(name)=="+vv+"]")
					else:
						i = moose.wildcardFind(moose.element(modelObjectlist).path+"/##/"+vv)
					
					if len(i) == 0:
						itemsnotfound[vv]={k:v}
					elif(len(i) > 1):
						multipleitems[vv]=i
	
	if itemsnotfound:
		print("These items defined in map file does not found in the model",)
		for k,v in itemsnotfound.items():
			print(k, "in ",v)
	if multipleitems:
		print("Ambiguous name exist in model ",)
		for k1,v1 in multipleitems.items():
			print(k1, "in ",v1)
# Run the 'main' if this script is executed standalone.
if __name__ == '__main__':
	main()
