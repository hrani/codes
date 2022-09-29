# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth
# Floor, Boston, MA 02110-1301, USA.
# 

'''
*******************************************************************
 * File:            htgraph.py
 * Description:
 * Author:          G.V. Harsha Rani, Upinder S. Bhalla
 * E-mail:          hrani@ncbs.res.in, bhalla@ncbs.res.in
 ********************************************************************/

/**********************************************************************
** This program converts HILLTAU models defined in JSON format to 
** reaction diagrams. It draws on the 'dot' program for graphical layout
** of networks.
**           copyright (C) 2021 Harsha Rani, Upinder S. Bhalla. and NCBS
**********************************************************************/
'''
'''
2021
Mar 31 
sub->prd is connected with double arrow,
ligand with single arrow, 
inhibit with tee and 
modifier with diamond
legends and constant are added

Apr 7: group is added

Apr 8: eqns added with pluse and sigma 

May 4: HillTau API is called for reading json file

May 10:
added matplotlib for getting colors
validation of input file type is done
output image can be saved as png or svg

May 15: set function is remove to get Unique items
May 24: Margin for the cluster is increased from 8 to 22. 

May 31: line flags for eliminating the legend and for changing colors to bw.

June 1: added features for adjusting fontsize and height on command line

June 3: Group colors and node colors added

June 15: more option which are Optional 
	-ranksep'   : set rank separation (vertical spacing) in output.
	-group'     : Display group pass multiple group name with comma seperator
	-fontsize'  : set font size for node labels.
	-no_legend' : Turns off generabnc    tion of legend'
	-bw'		: Goes into black-and-white plotting mode

Jun 19: Order of molecules
	#mol Kmod inhibit First-element second-element third-element
	 2	  0	     0 	     Input       Activator         --
	 2    0      1       Input       Inhibitor         --
	 3    1      0       Input       Modifier 		Activator
	 3    1      1       Input       Modifier       Inhibitor

Jun 30: with option -sg or --specific group, one can display specific group from the big model
python htgraph.py model.json -sg "group1_g","group2_g"
- If group name doesn't exist then it just ignores that specific group and display rest 
- If no group, specified in the command line exist then entire model is display like wise if no group is specified then
also entire model is displayed. 

'''

import sys,os
#sys.path.insert(1, 'PythonCode/')
from subprocess import call
import matplotlib
from collections import OrderedDict
import subprocess
from hillTau import *


use_bw = False

matplotcolors = []
for name,hexno in matplotlib.colors.cnames.items():
	matplotcolors.append(name)


def countX(lst, x):
	return lst.count(x)

def unique(list1):
	output = []
	for x in list1:
		if x not in output:
			output.append(x)
	return output
	#return list(set(list1))

def checkdigit(startstringdigit,grp,sp):
	if sp.startswith(tuple('0123456789')):
		if grp in startstringdigit:
			#pass#startstringdigit[grp].append((sp:"s"+sp))
			startstringdigit[grp][sp] = "s"+sp
		else:
			startstringdigit[grp] ={sp:"s"+sp}

def checkdigitEqu(startstringdigit,grp,sp):
	if grp in startstringdigit:
		grpitems = startstringdigit[grp]
		for k,v in grpitems.items():
			if k == sp:
				sp = v
	return(sp)

def getColor(gIndex,fwd_rev="forward"):
	if use_bw:
		return( "black", gIndex )

	if gIndex < len(matplotcolors):
		grpcolor = matplotcolors[gIndex]
		if grpcolor in ["white","wheat","aqua","whitesmoke","mintcream","oldlace","black","snow","aliceblue","azure","cornsilk","beige","bisque","blanchedalmond","antiquewhite","lightyellow","lightsteelblue","ghostwhite","floralwhite","ivory","honeydew"]:#mediumpurple","mediumvioletred","mediumseagreen"]:
			if fwd_rev == "reverse":
				gIndex = gIndex -1
			else:
				gIndex = gIndex +1

			return getColor(gIndex,fwd_rev)
		else:
			if fwd_rev == "reverse":
				gIndex = gIndex -1
			else:
				gIndex = gIndex +1
			return(grpcolor,gIndex)
	else:
		return getColor(0)

def jsontoDb(modelfilename,modelpath, glob_constant,outputfile):#modelpath, outputfile, ranksep = 0, hasLegend = True, fontsize = 18, showGroups = True,specific_group = []):
	group_no = 0;
	#groupmap = dict()
	groupmap = OrderedDict()
	global startstringdigit
	startstringdigit= OrderedDict()
	global node_color
	node_color = {}
	lig_exist = False
	kmod_exist = False
	inhibit_exist = False
	edge_arrowsize = 1.5
	edge_weight = 1
	s = ""
	st = os.path.splitext(outputfile)
	outputfilename = st[0]
	dirname, basename = os.path.split(outputfile)
	accessname  = (os.path.splitext(basename))[0]
	#command = "dot -T"+ outputfiletype + " "+ outputfilename+".dot -o "+outputfile
	#call([command], shell=True)
	#print(" st ",dirname,accessname)
	pngoutput = dirname+"/png/"+accessname+".png"
	#print(pngoutput)
	command = "python3 ~/AutSim/HillTau/htgraph.py "+ modelfilename +" -o " +pngoutput

	#print(command)
	filename =call([command], shell=True)
	#print("sql ->htgraph filename", filename)
	#print("here ",dirname,basename,accessname)
	if len( st ) > 1: 
		outputfiletype = st[1][1:]
	else:
		outputfiletype = "sql"
	
	f_doqcs = open(outputfilename+".sql", "w")
	#print(accessname)
	species = "Rodent"
	tissue = "Brain"
	cellcompartment = "Synapse"
	methodology = "Reduced quantitative fit to experiments"
	source = "<a href=https://doi.org/10.1371/journal.pcbi.1009621>Bhalla US. HillTau: A fast, compact abstraction for model reduction in biochemical signaling networks. PLoS Comput Biol 2021 Nov 29;17(11):e1009621. </a>"
	model_implementation = "HillTau implementation"
	notes= "This model was used to generate fig 2 supp A to F in <a href=https://doi.org/10.1371/journal.pcbi.1009621>Bhalla US. HillTau: A fast, compact abstraction for model reduction in biochemical signaling networks. PLoS Comput Biol 2021 Nov 29;17(11):e1009621.</a>"
	modeltype = "HT"
	model_validation = " "
	f_doqcs.write("USE doqcs\n");
	f_doqcs.write("INSERT INTO accession(\n");
	f_doqcs.write( "\tname,accesstype,entrydate,transcriber,developer,\n");
	f_doqcs.write( "\tspecies, tissue, cellcompartment, \n");
	f_doqcs.write( "\tmethodology, source, model_implementation, \n");
	f_doqcs.write( "\tmodel_validation, notes,modeltype)\n");

	f_doqcs.write( "VALUES(" );
	f_doqcs.write( "\n\t\""+accessname+"\", \"Network\", Now(), \"Upinder S. Bhalla, NCBS\", \"Upinder S. Bhalla, NCBS\",");
	f_doqcs.write( "\n\t \""+species+"\",\""+tissue+"\",\""+cellcompartment+"\",");
	f_doqcs.write( "\n\t\""+methodology+"\",\""+source+"\",\""+model_implementation+"\", ");
	f_doqcs.write( "\n\t\""+model_validation+"\",\""+notes+"\",\""+modeltype+"\");");

	f_doqcs.write( "\nSELECT @accessno := last_insert_id();\n");

	uploadpng = "/var/lib/mysql-files/"+accessname+".png"
		
	f_doqcs.write( "UPDATE accession SET figure = LOAD_FILE("+"'"+uploadpng+"'"+") WHERE accessno = @accessno;\n");

	#f_doqcs.write( "\nINSERT into modelfiles(accessno, model, format, is_native)\n");
	#f_doqcs.write( "VALUES(@accessno, LOAD_FILE("+"'"+"/var/lib/mysql-files/"+accessname+".xml"+"'"+"), 'SBML', '0' );\n");

	f_doqcs.write( "\nINSERT into modelfiles(accessno, model, format, is_native)\n");
	f_doqcs.write( "VALUES(@accessno, LOAD_FILE("+"'"+"/var/lib/mysql-files/"+accessname+".json"+"'"+"), 'JSON', '0' );\n");
	#f_doqcs.write("\nINSERT into geometry(accessno,name,size)\n");
	#f_doqcs.write("VALUES(@accessno,'"+"geometry"+"',"+"1e-18"+");\n");
	#f_doqcs.write( "\nSELECT @comptid0 := last_insert_id();\n");
	#print(" glob_constant ",glob_constant)
	if glob_constant:
		for f,v in glob_constant.items():
			#print("#",f,v)
			f_doqcs.write( "\nINSERT into constant(accessno, name,value) ");
			f_doqcs.write("VALUES(@accessno,'"+str(f)+"',"+str(v)+");\n");
		
	displayGroups = []
	#print("203 grpInfo ",modelpath.grpInfo)
	# if specific_group == None:
	# 	displayGroups = modelpath.grpInfo
	# else:
	# 	if any(i in specific_group for i in modelpath.grpInfo):
	# 		displayGroups = specific_group
	# 	else:
	# 		displayGroups = modelpath.grpInfo
	displayGroups = modelpath.grpInfo
	#print("## displayGroups ",displayGroups)
	specieslist = writeSpecies(modelpath)
	#print("specieslist 2  ",specieslist)
	funclist = writeFunc(modelpath)
	#print("\n funclist ",funclist)
	writeReac(modelpath,groupmap,displayGroups)
	#print("---- ",groupmap)
	#print("------------------Start ----------------")
				
	for grp in displayGroups:
		#print (grp)
		f_doqcs.write("\nINSERT INTO pathway");
		f_doqcs.write("(accessno, name, notes)");
		f_doqcs.write( "	VALUES(");
		f_doqcs.write( "@accessno,'" + grp + "',\"\");");
		f_doqcs.write("\nSELECT @pathwayno := last_insert_id();\n");

		pngoutput = dirname+"/png/"+accessname+"_"+grp+".png"
		command = "python3 ~/AutSim/HillTau/htgraph.py "+ modelfilename +" -sg "+grp+" -o " +pngoutput
		print(command)
		call([command], shell=True)
		uploadpng = "/var/lib/mysql-files/"+accessname+"_"+grp+".png"
		f_doqcs.write("UPDATE pathway SET figure = LOAD_FILE("+"'"+uploadpng+"'"+") WHERE accessno = @accessno AND pathwayno = @pathwayno;\n");
		
		if grp in specieslist:
			for name,convertmillitomicro in specieslist[grp]:
				#f_doqcs.write("\nINSERT INTO molecule(name, accessno, pathwayno, is_buffered, initial_conc, issumtotal_available, D, molwt, Gid, notes)");
				f_doqcs.write("\nINSERT INTO molecule(name, accessno, pathwayno, is_buffered, initial_conc, issumtotal_available, D, molwt, notes)");
				
				f_doqcs.write("  VALUES('"+name+"', @accessno, @pathwayno,'");
				f_doqcs.write('0'+"',"+str(convertmillitomicro)+",'");
				found_name = False
				sumtotal = '0'
				if grp in funclist:
					for fl in funclist[grp]:
						if fl["name"] == name:
							#if name in [name for name,expession in funclist[grp]]:
							found_name = True
							sumtotal = '1'
				#f_doqcs.write(sumtotal+"','0',0, @comptid0,\" \");");
				f_doqcs.write(sumtotal+"','0',0,\" \");");
		if grp in groupmap:
			for l in groupmap[grp]:
				reactionkey = "\nINSERT INTO ht_reaction(accessno,pathwayno"
				reactionvalue = ""
				reactionvalue = "@accessno,@pathwayno"
				for k,v in l.items():
					if k.lower() == "prd":
						reactionkey = reactionkey+", prd"
						reactionvalue=reactionvalue+",'"+str(v)+"'"
					if k.lower() == "sub":
						reactionkey = reactionkey+", input"
						reactionvalue=reactionvalue+",'"+str(v)+"'"
					if k.lower() == "modifier":
						reactionkey = reactionkey+", Modifier"
						reactionvalue=reactionvalue+",'"+str(v)+"'"
					if k.lower() == "activator":
						reactionkey = reactionkey+", ligant,stoich_ligant"
						reactionvalue=reactionvalue+",'"+str(v['s'])+"','"+str(v['count'])+"'"
					if k.lower() == "inhibit":
						reactionkey = reactionkey+", ligant,stoich_ligant,inhibit"
						reactionvalue=reactionvalue+",'"+str(v['s'])+"','"+str(v['count'])+"','1'"
					if k.lower() == "tau":
						reactionkey = reactionkey+",tau"
						reactionvalue = reactionvalue+",'"+str(v)+"'"
					if k.lower() == "tau2":
						reactionkey = reactionkey+",tau2"
						reactionvalue = reactionvalue+",'"+str(v)+"'"
					if k.lower() == "ka":
						reactionkey = reactionkey+",KA"
						reactionvalue = reactionvalue+",'"+str(v)+"'"
					if k.lower() == "baseline":
						reactionkey = reactionkey+",baseline"
						reactionvalue = reactionvalue+",'"+str(v)+"'"
					if k.lower() == "gain":
						reactionkey = reactionkey+",gain"
						reactionvalue = reactionvalue+",'"+str(v)+"'"
					if k.lower() == "amod":
						reactionkey = reactionkey+",Amod"
						reactionvalue = reactionvalue+",'"+str(v)+"'"
					if k.lower() == "kmod":
						reactionkey = reactionkey+",Kmod"
						reactionvalue = reactionvalue+",'"+str(v)+"'"
					if k.lower() == "nmod":
						reactionkey = reactionkey+",Nmod"
						reactionvalue = reactionvalue+",'"+str(v)+"'"
				#print("\n \n reactionkey", reactionkey,"\nreactionvalue",reactionvalue)
				f_doqcs.write(reactionkey+")")
				f_doqcs.write(" VALUES("+reactionvalue+");\n")

	for grp in displayGroups:
		#Sep6 Check what to do with expression table 
		#	Thinking eqnstr will be in 'Expression' field
		#   'Constant' field will be set true or false, if constants exist then display from "Constant" Table "value" field will be filled 
		#	with global constants as its in the notes section in the browser
		#   "sub" field will be set true or false, based on this reedit the eqnstr to link the molecules with in the model in the browser and
		#    add in "moleculelist" all individual substrate
		#    
		if grp in funclist:
			exist_const = '0'
			#print ("funclist grp ",funclist[grp])
			f_doqcs.write("\nSELECT @pathwayno :="+"(SELECT pathwayno FROM `pathway` WHERE name='"+grp+"' and accessno = @accessno);");
			f_doqcs.write("\n SELECT @mre_id :="+"(SELECT id FROM `molecule` WHERE name='"+fl["name"]+"' and accessno = @accessno and pathwayno=@pathwayno);")
			for fl in funclist[grp]:
				#f_doqcs.write( "	VALUES(");
				print(dir(fl))
				print("fl[eqnstr] ",fl["eqnstr"],"C;",fl['consts'],"gl",glob_constant, "sub ",fl["mol"]) 
				if glob_constant != 0:
					exist_const = '1'
				if len(fl["mol"]) != 0:
					exist_sub = '1'
					#print("expression molecule",fl["mol"])
					for x in list(fl["mol"]):
						#print(x)
						f_doqcs.write("\n SELECT @m_id :="+"(SELECT id FROM `molecule` WHERE name='"+x+"' and accessno = @accessno);")
						f_doqcs.write("\nINSERT INTO moleculelist");
						f_doqcs.write("(mer_id, m_id,type) VALUES(@mre_id,@m_id,1);");	
				f_doqcs.write("\nINSERT INTO expression");
				f_doqcs.write("(accessno, pathwayno,m_id,expression,exist_const,exist_sub)");
				f_doqcs.write("VALUES(@accessno,@pathwayno,@mre_id"+",'" + str(fl["eqnstr"]) +"','"+exist_const+"','"+exist_sub+"');");

	
	'''
	# specielist,node_color = writeSpecies(modelpath,groupmap)
	funclist = writeFunc2(modelpath,groupmap,f_doqcs,edge_arrowsize,edge_weight, displayGroups, fontsize = fontsize - 2)
	edgelist,node_color,lig_exist,kmod_exist,inhibit_exist = writeReac(modelpath,groupmap,f_graph,edge_arrowsize,edge_weight,displayGroups,fontsize = fontsize - 2)
	# nIndex = len(matplotcolors)-1
	'''
	f_doqcs.write("\n")
	f_doqcs.close()
	
	# command = "dot -T"+ outputfiletype + " "+ outputfilename+".dot -o "+outputfile
	# call([command], shell=True)

def writeSpecies(modelpath):
	groupmap = {}
	for molname, mol in ( modelpath.molInfo.items() ):
		#print ("species ",molname,mol.grp)
		checkdigit(startstringdigit,mol.grp,molname)
		molname = checkdigitEqu(startstringdigit,mol.grp,molname)
		if mol.grp in groupmap:
			groupmap[mol.grp].append((molname,mol.concInit))
		else:
			groupmap[mol.grp] = [(molname,mol.concInit)]
	return groupmap

def writeFunc(modelpath):
	groupmap = {}
	for e,t in modelpath.eqnInfo.items():
		dict_str = ""
		if hasattr(t, 'consts'):
			dict_str = {"name":e,"eqnstr":t.eqnStr,"consts":t.consts,"mol":t.subs}
		else:
			dict_str = {"name":e,"eqnstr":t.eqnStr,"mol":t.subs}
		if t.grp in groupmap:
			groupmap[t.grp].append(dict_str)
			
		else:
			groupmap[t.grp] = [dict_str]
	return groupmap

def isfloat(n):
    try:
        float(n)
        return True
    except ValueError:
        return False
		
def writeReac(modelpath,groupmap,displayGroups):
	r = 0
	#print("-------Reaction/Enzyme ---------")
	for reacname, reac in ( modelpath.reacInfo.items() ):
		checkdigit(startstringdigit,reac.grp,reacname)
		prd = checkdigitEqu(startstringdigit,reac.grp,reacname)
		reacname = "r"+str(r)
		r = r+1
		if reac.grp in displayGroups:
			# if reac.grp in groupmap:
			# 	groupmap[reac.grp].append(["name":reacname,"prd":prd])
			# else:
			# 	groupmap[reac.grp] = [reacname]
			react_dict = {}
			react_dict = {"prd":prd}
			sublist = reac.subs
			sublistU = unique(reac.subs)
			#print("sublist ",sublist)		
			if (hasattr(reac,"tau")):
				react_dict["tau"]=reac.tau
			if hasattr(reac,"tau1"):
				react_dict["tau2"]=reac.tau2
			if hasattr(reac,"KA"):
				react_dict["KA"]=reac.KA
			if hasattr(reac,"baseline"):
				react_dict["baseline"]=reac.baseline
			if hasattr(reac,"gain"):
				react_dict["gain"]=reac.gain
			if hasattr(reac,"Amod"):
				react_dict["Amod"]=reac.Amod
			if hasattr(reac,"Kmod"):
				react_dict["Kmod"]=reac.Kmod
			if hasattr(reac,"Nmod"):
				react_dict["Nmod"]=reac.Nmod
			for sub in sublistU:
				newsub = sub
				if sub in startstringdigit:
					newsub = startstringdigit[sub]
				checkdigit(startstringdigit,reac.grp,sub)

				if (reac.inhibit == 1.0 and sublistU.index(sub) == len(sublistU)-1 ) :
					#c = countX(sublist,sub)
					sub = checkdigitEqu(startstringdigit,reac.grp,sub)
					react_dict["inhibit"] = {"s":sub,"count":reac.HillCoeff}
				elif len(sublistU) == 3 and sublist.index(sub) == 1:
					''' kmod Modulator odiamond '''
					sub = checkdigitEqu(startstringdigit,reac.grp,sub)
					react_dict["Modifier"] = sub
				else:
					if sublist.index(sub) >= 1:
						#c = countX(sublist,sub)
						lig_exist = True
						sub = checkdigitEqu(startstringdigit,reac.grp,sub)
						react_dict["activator"] = {"s":sub,"count":reac.HillCoeff}
					else:
						if  sublist.index(sub) == 0:
							''' input '''
							sub = checkdigitEqu(startstringdigit,reac.grp,sub)
							react_dict["sub"] =sub
			#print(" reac.group ",reac.grp,react_dict)
			if reac.grp in groupmap:

				groupmap[reac.grp].append(react_dict)
			else:
				groupmap[reac.grp] = [react_dict]

	return
		

def file_choices(choices,fname,iotype):
	ext = (os.path.splitext(fname)[1][1:]).lower()
	if iotype == "outputfile":
		if ext not in choices:
			parser.error("Requires output filetype {}".format(choices))
	else:
		if ext != "json":
			parser.error("Requires HillTau file in JSON format ")
			
	return fname

if __name__ == "__main__":

	parser = argparse.ArgumentParser( description = 'This program generates a db list for a HillTau model. It converts the specified HillTau file in JSON format, to the db format for the database\n')
	parser.add_argument('model',type=lambda s:file_choices(("json"),s,"input"), help='Required: filename of model, in JSON format.')
	#parser.add_argument( 'model', type = str, help='Required: filename of model, in JSON format.')
	parser.add_argument( '-o', '--output', type=lambda out:file_choices(("sql"),out,"outputfile"), help='Optional: writes out the sql model into named file. default takes json filename')
	parser.add_argument( '-r', '--ranksep', type=float, default = 0, help='Optional: set rank separation (vertical spacing) in output.')
	parser.add_argument( '-fs', '--fontsize', type=float, default = 18, help='Optional: set font size for node labels.')
	parser.add_argument( '-nl', '--no_legend', action='store_true', help='Optional: Turns off generation of legend')
	parser.add_argument( '-ng', '--no_groups', action='store_true', help='Optional: Removes grouping. All molecules and reactions sit together.')
	parser.add_argument( '-bw', '--bw', action='store_true', help='Optional: Goes into black-and-white plotting mode')
	parser.add_argument('-sg', '--specific_group', help='Optional: Specfiy group names for display,delimited groupname seprated by comma.',type=lambda s:s.split(","))
	args = parser.parse_args()
	use_bw = args.bw

	if args.output == None:
		dirpath = os.path.dirname(args.model)
		basename = os.path.basename(args.model)
	
	exttype = "sql"
	dirpath = os.path.dirname(args.output)
	basename1 = os.path.basename(args.output)
	if basename1 != "":
		st = os.path.splitext(basename1)
		if len(st[0]) == 1:	
			if not st[0][-1].isalpha():
				outputfilename = st[0][0:len(st[0])-1]
				if len(outputfilename) ==0:
					outputfilename = basestr[0]
			else:
				outputfilename = st[0]
		else:
			outputfilename = st[0]
	else:
		outputfilename = os.path.splitext(os.path.basename(args.model))[0]
		exttype = "sql"
	if dirpath != '/':
		outputfile = dirpath+"/"+outputfilename+"."+exttype
	else:
		outputfile = stt+"."+exttype
	jsonDict = loadHillTau( args.model )
	#print("@@@",jsonDict.get( "QuantityUnits" ))
	qu = jsonDict.get( "QuantityUnits" )
	qs = 1.0
	if not qu:
		qu = jsonDict.get( "quantityUnits" )
	if qu:
		qs = lookupQuantityScale[qu]
	glob_constant = jsonDict.get("Constants")
	modelpath = parseModel( jsonDict )
	jsontoDb(args.model,modelpath, glob_constant,outputfile )
