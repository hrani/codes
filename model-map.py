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
	else:
		if ext not in ["json",'xml','g']:
			parser.error("Requires HillTau file in JSON format ")
			
	return fname


class DuplicateTrackingHook:
	def __init__(self):
		self.multipletimes_keyval = []
		self.multiplevalues_forkey = {}
	def __call__(self, pairs):
		obj = dict(pairs)  # Construct the object
		result = dict()

		for key,value in pairs:
			if key in result and set(value) != set(result[key]):
				self.multiplevalues_forkey[key]=(result[key][0],value[0])
			if key in result and value == result[key]:
				self.multipletimes_keyval.append(" {} {} ".format(key,value))
			else:
				if key.lower() not in ['filetype','version']:
					result[key] = value
		return result


def checkmapfile(mapfile):
	dupli_checking_hook = DuplicateTrackingHook()

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
					decoder = json.JSONDecoder(object_pairs_hook=dupli_checking_hook)
					maplist = decoder.decode(json_file.read())

					return maplist,dupli_checking_hook.multipletimes_keyval,dupli_checking_hook.multiplevalues_forkey
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
			
			modelitemsdict = {}
			modelitemsdict.update(htmodel.molInfo)
			modelitemsdict.update(htmodel.reacInfo)
			modelitemsdict.update(htmodel.eqnInfo)
			modelitemsdict.update(htmodel.namedConsts) 
			
			modelitems = htmodel.grpInfo
			for k,v in modelitemsdict.items():
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
	mapmodelObjectlist,multipletimes_keyval,multiplevalues_forkey = checkmapfile(args.map)
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
	
	return multipletimes_keyval,multiplevalues_forkey,itemsnotfound,multipleitems
# Run the 'main' if this script is executed standalone.
if __name__ == '__main__':
	multipletimes_keyval,multiplevalues_forkey,itemsnotfound,multipleitems=main()
	if multipletimes_keyval:
		print("These key:value pair has definded multiple times in the map file:\n",multipletimes_keyval)
	if multiplevalues_forkey:
		print("\nThese key has being assigned to 2 different values with in the map file which is not allowed:\n",multiplevalues_forkey)
	if itemsnotfound:
		print("\nThese items defined in value field in the map file doesn't exist in the model:\n",itemsnotfound)
	if multipleitems:
		print("\nAmbiguous name exist in model:\n",multipleitems)