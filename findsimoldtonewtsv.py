#Here the findsim old tsv file are converted to new format which is not be moose specific path.
# 1.stimuli->entities,
# 2.readout->entities are taken from Json file,
# 3.In modelmapping, 
#   a.modelSource,filename,citationid,citation,author,modellookup,scoringFormula is removed 
#   b.itemstodelete,modelSubset,paramterChange are taken from json file

import csv
import json
import os
import argparse


def main():
    """ This program convertes findsim old tsv file format to new format with JSON.
    """
    parser = argparse.ArgumentParser( description = 'This is will convert old findsim experiment sheet to new format on the basis of JSon ')
    #parser.add_argument( 'script', type = str, help='Required: Json file and atleast one filename of experiment spec, in tsv format or directory in which experiment file resides.')
    parser.add_argument( '-j', '--json', type = str, help='Json filename, .json' )
    parser.add_argument( '-p', '--path', type = str, help='path of a directory or tsv file')

    #parser.add_argument( '-d', '--directory', type = str, help='directory path where experiment sheet file')
    #parser.add_argument( '-f', '--file', type = str, help='Opitional tsvfile' )
    parser.add_argument( '-o', '--output', type = str, help='Optional: directory path to save output tsv files',default="/tmp" )
    args = parser.parse_args()
    jsonfile = ""
    tsvfile = ""
    global jsondata
    if args.output:
    	output = args.output
    else:
    	output = "/tmp"
    if args.json:
	    try:
	        with open(args.json, 'r') as json_file:
	        	jsondata = json.load(json_file)
	        	print(" json file exist ")          
	    except  IOError:
	        print "json file doesn't exist, the program exiting"
	        exit()
    else:
    	print "json file doesn't exist, the program exiting"
    	exit()

    if not args.path:
    	print ("path to directory or tsv file not provided, the program exiting")
    	exit()
    else:
    	if os.path.isdir(args.path):  
	    		for f in os.listdir(args.path):
	    			tsvfile = f
	    			if '.tsv' in tsvfile:
	    				convert(jsondata,args.path+tsvfile,output)
    	
    	elif os.path.isfile(args.path):
    		try:
    			with open(args.path, 'r') as fn:
    				tsvfile = args.path
    				convert(jsondata,tsvfile,output)
    		except  IOError:
    			print "tsv file doesn't exist in "+args.file+", the program exiting"
    			exit()	 
    	else:
    		print " directory /tsv file for running is not provided, the program exiting"
    		exit()
    
 
def convert(json_file,tsv_file,output="/tmp"):	
	fname = tsv_file
	modelLookup = {}
	mlookup = []
	writetofile = ""
	head = ""
	print("reading ", os.path.basename(fname))
	with open(fname, 'rt') as fs:
		reader = csv.reader(fs, delimiter='\t')
		for row in reader:
			if row:
				if 'modelLookup' == row[0]:
					if row[1]:
						mlookup = row[1].split( ',' )
					if ':' in mlookup:
						modelLookup =  { i.split(':')[0]:i.split(':')[1].strip() for i in mlookup }
						break;
	with open(fname, 'rt') as f:
		reader = csv.reader(f,delimiter="\t")
		#print reader[0]
		errorlist = ""
		para = []
		datawidth = 2
		lineno = 0
		notfound =""
		for line in f:
			#print row
			row = line.split('\t')
			#print fs
			
			lineno +=1
			#print row.split('\t')
			#head = row[0]
			if row[0] == "Experiment metadata" or row[0] == "Stimuli" or row[0] =="Readouts" or row[0] == "Model mapping":
				head = row[0]
				writetofile += "\n"+row[0]+"\n"
			
			elif  row[0] == 'entities' or row[0] == 'ratioReferenceEntities':
				key = ""
				writetofile +=row[0]+"\t"
				#print( "50 ",row[1])
				if row[1]:
					if modelLookup.has_key(row[1]):
						key = [k for k, v in jsondata.iteritems() if modelLookup[row[1]] in v]
						if key:
							if key[0] == row[1]:
								writetofile += key[0]+"\n"
							else:
								print head +" --> ",row[0]
								print " key in model lookup and json file not same modellookup -> key \'",row[1],"\' and json \'",key[0], "\'for the value \'",modelLookup[row[1]],"\' and Json \'",jsondata[key[0]][0],
								print "\' writing as its"
								writetofile += row[1]+"\n"

						else:
							print head, "-->",row[0]," --- No entry found in Json file for in -->",row[1], "writing as its"
							writetofile += row[1]+"\n"		
					else:
						print head, "-->", row[1],"key not found in modelLookup"
						writetofile += "\n"
				else:
					writetofile +="\t\n"

			elif 'modelSubset' == row[0] or 'itemstodelete' == row[0]:	
				writetofile +=row[0]+"\t"
				if row[1]:
					modss = row[1].split(',')
					jsonlist,errorlist = slicelist(modss)
					writetofile += jsonlist+"\n"
					if errorlist != "":
						print "keys not found in -->",head, "-->",row[0], errorlist
				else:
					print row[0], " is empty",
					writetofile += "\t\n"

			elif 'parameterChange' == row[0]:
				readParameter(f,para,datawidth,lineno)
				modelparChange = ""
				writetofile += "\t".join(row)+"Object\tparameter\tValue\n"

				for rr in para:
					if rr[0] != "" or rr[0] != None:
						key = [k for k, v in jsondata.iteritems() if rr[0] in v]
						if key:
							modelparChange +=key[0]+'\t'+str(rr[1])+'\t'+str(rr[2])+"\n"
						else:
							notfound +="\n"+rr[0]
							print ("Parameter change key not found ", rr[0])
				writetofile += modelparChange
			else:
				if head == "Model mapping": 
					if (row[0] not in ['modelSource','fileName','citationId','citation','authors','scoringFormula','solver','modelLookup']):
						writetofile = writetofile+"\t".join(row)
				else:
					writetofile = writetofile+"\t".join(row)
	fileName = os.path.basename(tsv_file)
	newfile = output+"/New_"+fileName
	print("\n")
	print(newfile)
	with open(newfile, "w") as wf:
		wf.write(writetofile)
	print("\n#############\n")
def readParameter(fd, para, width,lineno):
    for ll in fd:
    	cols = ll.split("\t")
        if len( cols ) == 0 or cols[0] == '' or cols[0].isspace():
            break
        if cols[0].strip().lower() == "object":
            continue
        row = []
        lcols = 0
        for c in cols:
            c = c.replace('\n',"")
            if c != '':
                if lcols > 1 :
                    row.append(float(c))
                else:
                    row.append(c)
                if len( row ) > width:
                    break;
                lcols = lcols+1
        para.append( row )
    

def slicelist(itemlist):
	notfound =""
	modelSubsetlist = ""
	i =1
	for ms in itemlist:
		ms = ms.strip(' \n\t\r')
		if ms != "" or ms != None:
			key = [k for k, v in jsondata.iteritems() if ms in v]
			if key:
				if i < len(itemlist):
					modelSubsetlist +=key[0]+','
				elif(i == len(itemlist)):
					modelSubsetlist +=key[0]
			else:
				notfound +="\n"+ms
			i +=1
	return modelSubsetlist, notfound

# Run the 'main' if this script is executed standalone.
if __name__ == '__main__':
	main()
