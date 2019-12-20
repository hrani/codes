import csv
import json
import os
import argparse
from moose import wildcardFind,element,Neutral,loadModel

def main():
    """ This program takes all the values from the tsv to write Json.
    """
    parser = argparse.ArgumentParser( description = 'This is will take values from the tsv files '
    )

    parser.add_argument( '-m', '--model', type = str, help='Genesis file, .g' )
    parser.add_argument( '-d', '--directory', type = str, help='directory path where experiment sheet file')
    parser.add_argument( '-f', '--file', type = str, help='Opitional tsvfile' )
    #parser.add_argument( '-o', '--output', type = str, help='Optional: directory path to save output tsv files',default="/tmp" )
    args = parser.parse_args()
    tsvfile = ""
    global GmodelLookup,Ggroup,Gentities

    GmodelLookup = {}
    Ggroup = {}
    Gentities = []
    #loadModel('/home/harsha/hrani_FindSim/FindSim_dec2019/models/synSynth7.g','/model',"ee")
    #convert("/home/harsha/hrani_FindSim/FindSim_dec2019/Curated/FindSim-Bhalla1999_fig2B.tsv",'/model')
    #convert("/home/harsha/hrani_FindSim/FindSim_dec2019/Curated/FindSim-Bhalla1999_fig2C.tsv",'/model')
    #convert("/home/harsha/hrani_FindSim/FindSim_dec2019/Curated/FindSim-Bhalla1999_fig4C.tsv",'/model')
    #print "moodelLookup: ", GmodelLookup,len(GmodelLookup)
    #print "group: ",Ggroup,len(Ggroup)
    #print "moose entities: ",Gentities,len(Gentities)
    #print "gmodel",GmodelLookup
    
    # if args.output:
    # 	output = args.output
    # else:
    # 	output = "/tmp"
    try:
        with open(args.model, 'r') as json_file:
        	loadModel(args.model,'/model',"ee")
        	print(os.path.basename(args.model)," model is loaded into moose under /model",)          
    except  IOError:
        print "genesis  file doesn't exist, the program exiting"
        exit()
        	
    if not args.file:
    	if 	args.directory:
    		for f in os.listdir(args.directory):
    			tsvfile = f
    			if '.tsv' in tsvfile:
    				#convert(jsondata,args.directory+tsvfile,output)
    				convert(args.directory+tsvfile,'/model')
    	else:
    		print " directory does not exist"
    		exit()
    else:
    	if args.directory:
    		try:
        		with open(args.directory+args.file, 'r') as fn:
			    	tsvfile = args.file 
			    	#convert(jsondata, args.directory+tsvfile,output)
			    	convert(args.directory+tsvfile,'/model')    
    		except  IOError:
        		print "tsv file doesn't exist "+args.directory+args.file+", the program exiting"
        		exit()
    	else:
    		try:
    			with open(args.file, 'r') as fn:
			    	tsvfile = args.file 
			    	#convert(jsondata,tsvfile,output)
			    	convert(tsvfile,'/model')   
    		except  IOError:
        		print "tsv file doesn't exist in "+args.file+", the program exiting"
        		exit()
	print "moodelLookup: ", GmodelLookup,len(GmodelLookup)
    print "group: ",Ggroup,len(Ggroup)
    print "moose entities: ",Gentities,len(Gentities)
    print "gmodel",GmodelLookup    
def convert(tsv_file,model):
	
#def convert(json_file,tsv_file,output="/tmp"):	
	fname = tsv_file

	writetofile = ""
	head = ""
	with open(fname, 'rt') as fs:
		reader = csv.reader(fs, delimiter='\t')
		for row in reader:
			if 'modelLookup' == row[0]:
				if row[1]:
					mlookup = row[1].split( ',' )
					
				if mlookup:
					modelLookup =  { i.split(':')[0]:i.split(':')[1].strip() for i in mlookup }
					for k,v in modelLookup.iteritems():
						if GmodelLookup.has_key(k):
							if GmodelLookup[k] != v:
								print("value in globalmodellookup is different from this tsvfile",os.path.basename(fname),k,v, "value in gobal",GmodelLookup[k])
						else:
							GmodelLookup[k] = v
					break;

	with open(fname, 'rt') as f:
		reader = csv.reader(f,delimiter="\t")
		errorlist = ""
		para = []
		datawidth = 2
		lineno = 0
		notfound =""
		for line in f:
			row = line.split('\t')
			lineno +=1
			if row[0] == "Stimuli" or row[0] =="Readouts" or row[0] == "Model mapping":
				head = row[0]
			elif  row[0] == 'entities' or row[0] == 'ratioReferenceEntities':
				key = ""
				if row[1]:
					if GmodelLookup.has_key(row[1]):
						if modelLookup.has_key(row[1]):
							if GmodelLookup[row[1]] != modelLookup[row[1]]:
								print ("value in globalmodellookup is not same in modellookup",
								row[1], "->",GmodelLookup[row[1]] ,modelLookup[row[1]])
					else:
							print head, "-->",row[0]," --- No entry found in  modellookupfile for in -->",row[1], "writing as its"
									
				elif row[1] != "":
					print "check ",head, row[0],row[1]
				
			
			elif 'modelSubset' == row[0] or 'itemstodelete' == row[0]:	
				if row[1]:
					modss = row[1].split(',')
					group,others = slicelist(modss,model)
					for k,v in group.iteritems():
						if Ggroup.has_key(k):
							if Ggroup[k] != v:
								print("value in globalmodellookup is different from this tsvfile",os.path.basename(fname),k,v, "value in gobal",Ggroup[k])
						else:
							Ggroup[k] = v
					for l in others:
						if l not in Gentities:
							Gentities.append(l)

				else:
					print row[0], " is empty",
					
			elif 'parameterChange' == row[0]:
				readParameter(f,para,datawidth,lineno)
				for p in para:
					if p not in Gentities:
						Gentities.append(p)

def readParameter(fd, para, width,lineno):
    for ll in fd:
    	cols = ll.split("\t")
        if len( cols ) == 0 or cols[0] == '' or cols[0].isspace():
            break
        if cols[0].strip().lower() == "object":
            continue
        para.append(cols[0])
def slicelist(itemlist,model):
	notfound =""
	modelSubsetlist = ""
	i =1
	group = {}
	others = []
	for ms in itemlist:
		if ms != "" or ms != None:
			rootpath = model
			name = ms
			try1 = wildcardFind( rootpath+'/' + name )
	        try2 = wildcardFind( rootpath+'/##/' + name )
	        try2 = [ i for i in try2 if not '/model[0]/plots[0]' in i.path ]  
	        
	        if len( try1 ) + len( try2 ) > 1:
	            raise SimError( "findObj: ambiguous name: '{}'".format(name) )
	        if len( try1 ) + len( try2 ) == 0:
	            if noRaise:
	                print element('/')
	            else:
	                raise SimError( "findObj: No object found named: '{}'".format( name) )
	        if len( try1 ) == 1:
	            if(try1[0].className== 'Neutral'):
	            	cc = element(try1[0]).name+"_g"
	            	group[cc] = name
	            else:
	            	others.append(name)
	            
	        else:
	            if(try2[0].className== 'Neutral'):
	            	cc = element(try2[0]).name+"_g"
	            	group[cc] = name
	            else:
	            	others.append(name)
	return group,others

# Run the 'main' if this script is executed standalone.
if __name__ == '__main__':
	main()
