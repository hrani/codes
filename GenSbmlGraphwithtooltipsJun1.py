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
 * File:            GenSbmlGraph.py
 * Description:
 * Author:          G.V. Harsha Rani
 * E-mail:          hrani@ncbs.res.in
 ********************************************************************/

/**********************************************************************
** This program converts Genesis/SBML models reaction diagrams. 
** It draws on the 'dot' program for graphical layout of networks.
**           copyright (C) 2021 Harsha Rani
**********************************************************************/
'''

import sys,os
from subprocess import call
import matplotlib
from collections import OrderedDict
import argparse 


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

def populate(startstringdigit,grp,sp,nsp):
	if grp in startstringdigit:
		startstringdigit[grp][sp] = nsp 
	else:
		startstringdigit[grp] ={sp:nsp}

def checkSpecialChar(startstringdigit,grp,sp):
	found = False
	spOrg = sp
	if sp.find('.') != -1:
		sp = sp.replace(".","_dot_")
		found = True
	if sp.find("(") != -1:
		sp = sp.replace("(","_bo_")
		found = True
	if sp.find(")") != -1:
		sp = sp.replace(")","_bc_")
		found = True
	if sp.find("_") != -1:
		sp = sp.replace(")","_un_")
		found = True
	if sp.startswith(tuple('0123456789')):
		sp = "s"+sp
		found = True
	# if spOrg != sp:
	# 	print("81 ",grp,spOrg,sp)
	if found:
		populate(startstringdigit,grp,spOrg,sp)
	
	
def checkdigitEqu(startstringdigit,grp,sp1):
	sp = sp1.name
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

def jsontoPng(modelpath, outputfile, ranksep = 0, hasLegend = True, fontsize = 18, showGroups = True,specific_group = []):
	group_no = 0;
	groupmap = OrderedDict()
	global startstringdigit
	startstringdigit= OrderedDict()
	global node_color
	node_color = {}
	edge_arrowsize = 1.5
	edge_weight = 1
	s = ""
	
	st = os.path.splitext(outputfile)
	outputfilename = st[0]
	if len( st ) > 1: 
		outputfiletype = st[1][1:]
	else:
		outputfiletype = "png"
	
	f_graph = open(outputfilename+".dot", "w")
	f_graph.write("digraph mygraph {\n\trank=TB;\n")
	if ranksep > 0.0:
		f_graph.write("\tranksep={};\n".format( ranksep ))
	#f_graph.write("ratio = 1.0\n")
	#f_graph.write("ratio = \"fill\"\n")
	#f_graph.write("size = \"4,4!\"\n")
	#f_graph.write("node [shape=box, penwidth=2, height=0.01, width=0.01 ];")
	f_graph.write("node [shape=box, penwidth=2,fontsize={}];".format( fontsize ) )
	
	
	displayGroups = []
	print("specific_group ",specific_group)
	if specific_group == None:
		displayGroups = moose.wildcardFind(modelpath.path+"/##[TYPE=Neutral],/##[ISA=ChemCompt]")
		
	else:
		foundgroup = False
		fgrp = ""
		for i in specific_group:
			print("i ",i,moose.wildcardFind(modelpath.path+'/##[FIELD(name)=='+i+')'))
			fgrp =  moose.wildcardFind(modelpath.path+'/##[FIELD(name)=='+i+')')
			if fgrp:
				print(fgrp)
				foundgroup = True
				for fgrp1  in fgrp:
					displayGroups.append(fgrp1)
				print(displayGroups)
		if not foundgroup:
			displayGroups = moose.wildcardFind('/modelpath/##[ISA=Neutral]')
		
	print("132 ",displayGroups)
	writeSpecies(modelpath,groupmap,displayGroups)
	chanlist = writechan(modelpath,groupmap,f_graph,edge_arrowsize,edge_weight, displayGroups, fontsize = fontsize - 2)
	
	funclist = writeFunc(modelpath,groupmap,f_graph,edge_arrowsize,edge_weight, displayGroups, fontsize = fontsize - 2)
	reaclist,node_color= writeReac(modelpath,groupmap,f_graph,edge_arrowsize,edge_weight,displayGroups,fontsize = fontsize - 2)
	enzlist,node_color = writeEnz(modelpath,groupmap,f_graph,edge_arrowsize,edge_weight,displayGroups,fontsize = fontsize - 2)
	nIndex = len(matplotcolors)-1
	
	if showGroups:
		color,nIndex = getColor(nIndex,"reverse")
		compt_no = 0
		group_no = 0
		for Cmpt in moose.wildcardFind(modelpath.path+'/##[ISA=ChemCompt]'):
			print("cmpt ",Cmpt)
			s = s+"\nsubgraph cluster_"+str(compt_no)+"\n{"+"\n"+"label=\""+Cmpt.name+"\";\npenwidth=4; margin=10.0\ncolor=\""+color+"\";\nfontsize="+str(fontsize + 2)+";\n"
			sps = ""
			if Cmpt in groupmap:
				items = groupmap[Cmpt]
				
				for sp in items:
						if items.index(sp) != 0:
							if type(sp) is tuple:
								#print("111--- ",sp[1].name,moose.element(sp[1]).concInit, moose.element(sp[1]).n )

								sps = sps+'\n'+sp[0]+' [label=\"'+moose.element(sp[1]).name+'\"'+',tooltip = \"concInit = '+str(float("{:.6f}".format(moose.element(sp[1]).concInit)))+'\nn = '+str(float("{:.6f}".format(moose.element(sp[1]).n)))+'"]'
							else:
								sps = sps+'\n'+sp
						else:
							if type(sp) is tuple:
								sps = sps+sp[0]+' [label=\"'+moose.element(sp[1]).name+'\"'+',tooltip = \"concInit = '+str(float("{:.6f}".format(moose.element(sp[1]).concInit)))+'\nn = '+str(float("{:.6f}".format(moose.element(sp[1]).n)))+'"]'
							else:
								sps = sps+'\n'+sp	
				s = s+sps
				compt_no += 1;
			for grp in moose.wildcardFind(Cmpt.path+'/##[TYPE=Neutral]'):
				if grp in groupmap:
					s = s + "\nsubgraph cluster_"+str(group_no)+"_"+str(compt_no)+"i\n{"
					s = s+"\nsubgraph cluster_"+str(group_no)+"\n{"+"\n"+"label=\""+grp.name+"\";\npenwidth=4; margin=10.0\ncolor=\""+color+"\";\nfontsize="+str(fontsize + 2)+";\n"
					sps = ""
					#print("groupmap ",groupmap,grp)
					items = groupmap[grp]
					items = list(unique(items))
					for sp in items:
						if items.index(sp) != 0:
							if type(sp) is tuple:
								sps = sps+'\n'+sp[0]+' [label=\"'+moose.element(sp[1]).name+'\"]'
							else:
								sps = sps+'\n'+sp
						else:
							if type(sp) is tuple:
								sps = sps+sp[0]+' [label=\"'+moose.element(sp[1]).name+'\"]'
							else:
								sps = sps+'\n'+sp	
					s = s+sps+"\n} style=invisible\n}"
					group_no += 1;
			s = s +"\n}"
			
	f_graph.write(s)
	f_graph.write(chanlist)
	f_graph.write(funclist)
	f_graph.write(reaclist)
	f_graph.write(enzlist)
	
	nodeIndex = 0
	for k,uu in groupmap.items():
		#print("s ",k)
		if k in displayGroups:
			
			for vl in uu:
				if type(vl) is tuple:
					l = vl[0]
					if l in node_color:
						v = node_color[l]
						v,nodeIndex = getColor(nodeIndex)
						ii = modelpath.path+'/##[FIELD(name)='+l+']'
						# for ll in moose.wildcardFind(ii):
						# 	if(len(ll.neighbors['nOut']) == 0  and  len(ll.neighbors['reacDest']) == 0):
						# 		#print("pools not connected to any obj ",ll)
						# 		pass
						# 	else:
						# 		f_graph.write("\n"+l+"[color=\""+v+"\"]")
						f_graph.write("\n"+l+"[color=\""+v+"\"]")
	for p,q in startstringdigit.items():
		if p in displayGroups:
			for m,n in q.items():
				pass #print("n",n,m)
				#f_graph.write("\n"+n+"[label=\""+m+"\"]")	
	
	if hasLegend:
		
		f_graph.write("\nnode [shape=plaintext]\nsubgraph cluster_01 {\n\tlabel = \"Legend\";\n\t{ rank=sink;\n\tkey [label=<<table border=\"0\" cellpadding=\"2\" cellspacing=\"0\" cellborder=\"0\">\n\t\
			<tr><td align=\"right\" port=\"i1\">Input/Output</td></tr>\
			<tr><td align=\"right\" port=\"i4\">Enzyme parent </td></tr>\
			\n")
		
		f_graph.write("\t</table>>]\n\tkey2 [label=<<table border=\"0\" cellpadding=\"2\" cellspacing=\"0\" cellborder=\"0\">\
			\n\t<tr><td port=\"i1\">&nbsp;</td></tr>\
			<tr><td port=\"i4\">&nbsp;</td></tr>\n")
		
		f_graph.write("\t</table>>]\n\tkey:i1:e -> key2:i1:w [arrowhead=normal  color=\"black:black\" style=bold]\
								   \n\tkey:i4:e -> key2:i4:w [arrowhead=none color=\"black:black\" style=bold]\
								   \n")

		
		f_graph.write("\t}\n\t}")

	f_graph.write("\n}")
	f_graph.close()
	
	command = "dot -T"+ outputfiletype + " "+ outputfilename+".dot -o "+outputfile
	call([command], shell=True)

def mooseIsInstance(melement, classNames):
	return moose.element(melement).className in classNames

def findGroup_compt(melement):
	while not (mooseIsInstance(melement, ["Neutral","CubeMesh", "CyclMesh"])):
		melement = melement.parent
	return melement

def writeSpecies(modelpath, groupmap,displayGroups):
	# getting all the species
	mIndex = 0
	numMol = 0 
	for mol in ( moose.wildcardFind(modelpath.path+'/##[ISA=PoolBase]') ):
		molname = mol.name
		molsize = "mol"+str(numMol)
		numMol+=1
		
		if (mol.parent).className != 'Enz':
			if displayGroups != "None":
				print("mol ",mol,mol.parent,"Display ",displayGroups)
				if mol.parent in displayGroups:
					molgrp = findGroup_compt(mol)
					checkSpecialChar(startstringdigit,molgrp,molname)
					molname = checkdigitEqu(startstringdigit,molgrp,mol)
					#print(" 288 ",mol)
					if molname not in node_color:
						spe_color,mIndex = getColor(mIndex)
						#node_color[molname] = spe_color
						node_color[molsize] = spe_color
					if molgrp in groupmap:
						groupmap[molgrp].append((molsize,mol))
					else:
						groupmap[molgrp] = [(molsize,mol)]
	
def writeFunc(modelpath,groupmap,f_graph,edge_arrowsize,edge_weight,displayGroups, fontsize = 16):
	equation_pluse = 0
	equation_sigma = 0
	edgelist = ""
	funcno = 0
	for t in moose.wildcardFind(modelpath.path+'/##[ISA=Function]'):
		tgrp = findGroup_compt(t)
		tname_no = t.name+str(funcno)
		funcno = funcno+1
		# checkSpecialChar(startstringdigit,tgrp,tname_no)
		# t.name = checkdigitEqu(startstringdigit,tgrp,tname_no)
		allpluse = True
		mathSym = []
		tt = moose.element(t.path+'/x')
		#print("tt ",moose.showfields(moose.element(t.path)))
		plusesize = "sigma"+str(equation_sigma)
		equation_sigma+=1
		print(tgrp)
		if tgrp in displayGroups:
			edgelist = edgelist+"\n"+plusesize+"[label=<&Sigma; "+str(equation_sigma)+ ">,shape=circle]"
		groupmap[tgrp].append(plusesize)
		for tsubs in unique(tt.neighbors['input']):
			#print(tsubs,tsubs.parent,groupmap[tsubs.parent])
			for oo in groupmap[tsubs.parent]:
				#print(" k, l ",oo)
				if type(oo) is tuple:
					if oo[1] == moose.element(tsubs):
						tsubs1 = oo[0]
						#print (tsubs1)
			input_color = node_color[tsubs1]
			#c = countX(tsubs,tt.neighbors['input'])
			c = 1
			tsubsname = checkdigitEqu(startstringdigit,tgrp,tsubs)

			if tgrp in displayGroups:
				#print(tsubsname,"#####",tsubs1)
				edgelist = edgelist+"\n"+tsubs1+"->"+plusesize+"[arrowhead=vee weight = "+str(edge_weight)+" color=\""+input_color+ "\" arrowsize = "+str(edge_arrowsize)+""
				if c > 1:
					edgelist = edgelist+ " label=\" "+str(c)+"\" fontsize={}".format( fontsize )
				edgelist = edgelist+"]"
		if tgrp in displayGroups:
			outputpool = t.neighbors['valueOut']
			#print("### ",outputpool,outputpool[0].parent)
			for oo in groupmap[outputpool[0].parent]:
				#print(" k, l ",oo)
				if type(oo) is tuple:
					if oo[1] == moose.element(outputpool[0]):
						outputpoolname1 = oo[0]
						#print (tsubs1)
			outputpoolname = checkdigitEqu(startstringdigit,tgrp,outputpool[0])
			edgelist = edgelist+"\n"+plusesize+"->"+outputpoolname1+"[arrowhead=vee weight = "+str(edge_weight)+ " arrowsize = "+str(edge_arrowsize)+"]"
	return edgelist
def writechan(modelpath,groupmap,f_graph,edge_arrowsize,edge_weight,displayGroups,fontsize = 16):
	chanlist = ""
	chan_color_list =[]
	sIndex = 0
	pIndex = 0 
	numchan = 0
	
	for chan in moose.wildcardFind( modelpath.path+'/##[ISA=ConcChan]'):
		changrp = findGroup_compt(chan)
		channame = chan.name 
		parChan = chan.parent
		chansize = "chan"+str(numchan)
		numchan+=1
		groupmap[changrp].append(chansize)
		for oo in groupmap[parChan.parent]:
				#print(" k, l ",oo)
				if type(oo) is tuple:
					if oo[1] == moose.element(parChan):
						parchan1 = oo[0]
				chanlist = chanlist+"\n"+chansize+"[label=<> shape=restrictionsite]"
		chanlist = chanlist+"\n"+parchan1+"->"+chansize+"[arrowhead=none]"
		inputChan = chan.neighbors['inPoolOut']
		for oo in groupmap[moose.element(inputChan[0]).parent]:
				if type(oo) is tuple:
					if oo[1] == moose.element(inputChan[0]):
						inpchan1 = oo[0]
		
		chanlist = chanlist +"\n"+inpchan1+"->"+chansize+"[arrowhead=normal weight = "+str(edge_weight)+ " arrowsize = "+str(edge_arrowsize)+ " style=bold]"
		outputChan = chan.neighbors['outPoolOut']
		for oo in groupmap[moose.element(outputChan[0]).parent]:
				if type(oo) is tuple:
					if oo[1] == moose.element(outputChan[0]):
						outchan1 = oo[0]
		
		chanlist = chanlist +"\n"+chansize+"->"+outchan1+"[arrowhead=normal weight = "+str(edge_weight)+ " arrowsize = "+str(edge_arrowsize)+ " style=bold]"
		
	return(chanlist)
	
def writeEnz(modelpath,groupmap,f_graph,edge_arrowsize,edge_weight,displayGroups,fontsize = 16):
	edgelist1 = ""
	enzyme_color_list =[]
	sIndex = 0
	pIndex = 0 
	numenz = 0
	for enz in moose.wildcardFind( modelpath.path+'/##[ISA=EnzBase]'):
		enzgrp = findGroup_compt(enz)
		enzname = enz.name
		sublist = enz.neighbors['sub']
		sublistU = unique(enz.neighbors['sub'])
		prdlist = enz.neighbors['prd']
		prdlistU = unique(enz.neighbors['prd'])
		enzsize = "Enz"+str(numenz)
		numenz+=1
		#print("enx ",moose.showfields(enz))
		edgelist1 = edgelist1+"\n"+enzsize+"[label=<""> shape=oval width=0.5 height=0.2, tooltip = \"km = "+str(float("{:.5f}".format(enz.Km)))+"\nkcat = "+str(float("{:.5f}".format(enz.kcat)))+"\"]"
		#edgelist1 = edgelist1+"\n"+enzsize+"[label=<^E^> shape=oval]"
		#edgelist1 = edgelist1+"\n"+enzsize+"[label=<^E^"+enzsize+"-"+enzname+"> shape=oval]"
		
		groupmap[enzgrp].append(enzsize)
		enzpar = checkdigitEqu(startstringdigit,enzgrp,enz.parent)
		#print("### ",enz.parent)
		for oo in groupmap[(enz.parent).parent]:
				#print(" k, l ",oo)
				if type(oo) is tuple:
					if oo[1] == moose.element(enz.parent):
						enzpar1 = oo[0]
		edgelist1 = edgelist1+"\n"+enzpar1+"->"+enzsize+"[arrowhead=none]"
				
		for sub in sublistU:
			c = countX(sublist,sub)
			for oo in groupmap[sub.parent]:
				#print(" k, l ",oo)
				if type(oo) is tuple:
					if oo[1] == moose.element(sub):
						newsub1 = oo[0]
			subname = sub.name
			subgrp = findGroup_compt(sub)	
			newsub = checkdigitEqu(startstringdigit,subgrp,sub)
			if newsub1 in node_color:
				enzyme_color = node_color[newsub1]
			else:
				enzyme_color,sIndex = getColor(sIndex)
				node_color[newsub1] = enzyme_color
			if c >1:
				edgelist1 = edgelist1 +"\n"+newsub1+"->"+enzsize+"[arrowhead=normal weight = "+str(edge_weight)+ " arrowsize = "+str(edge_arrowsize)+" color=\""+enzyme_color+":"+enzyme_color+"\""+ " label=\" "+str(c)+"\" style=bold]"
			else:
				edgelist1 = edgelist1 +"\n"+newsub1+"->"+enzsize+"[arrowhead=normal weight = "+str(edge_weight)+ " arrowsize = "+str(edge_arrowsize)+" color=\""+enzyme_color+":"+enzyme_color+"\" style=bold]"
		for prd in prdlistU:
			c = countX(prdlist,prd)
			for oo in groupmap[prd.parent]:
				#print(" k, l ",oo)
				if type(oo) is tuple:
					if oo[1] == moose.element(prd):
						newprd1 = oo[0]
			prdname = prd.name
			prdgrp = findGroup_compt(prd)
			newprd = checkdigitEqu(startstringdigit,prdgrp,prd)
			if newprd1 in node_color:
				enzyme_color = node_color[newprd1]
			else:
				enzyme_color,sIndex = getColor(sIndex)
				node_color[newprd1] = enzyme_color
			if c >1:
				edgelist1 = edgelist1 +"\n"+enzsize+"->"+newprd1+"[arrowhead=normal weight = "+str(edge_weight)+ " arrowsize = "+str(edge_arrowsize)+" color=\""+enzyme_color+":"+enzyme_color+"\""+ " label=\" "+str(c)+"\" style=bold]"
			else:
				edgelist1 = edgelist1 +"\n"+enzsize+"->"+newprd1+"[arrowhead=normal weight = "+str(edge_weight)+ " arrowsize = "+str(edge_arrowsize)+" color=\""+enzyme_color+":"+enzyme_color+"\" style=bold]"
	return(edgelist1,node_color)

def writeReac(modelpath,groupmap,f_graph,edge_arrowsize,edge_weight,displayGroups,fontsize = 16):
	edgelist1 = ""
	reaction_color_list =[]
	sIndex = 0
	pIndex = 0 
	numReac = 0
	
	for reac in moose.wildcardFind( modelpath.path+'/##[ISA=Reac]'):
		reacgrp = findGroup_compt(reac)
		reacname = reac.name
		sublist = reac.neighbors['sub']
		sublistU = unique(reac.neighbors['sub'])
		prdlist = reac.neighbors['prd']
		prdlistU = unique(reac.neighbors['prd'])
		reacsize = "Reac"+str(numReac)
		numReac+=1
		edgelist1 = edgelist1+"\n"+reacsize+"[label=<" "> shape=square,width=0.2, tooltip = \"kf = "+str(float("{:.2f}".format(reac.Kf)))+"\nkb = "+str(float("{:.2f}".format(reac.Kb)))+"\"]"
		groupmap[reacgrp].append(reacsize)
		for sub in sublistU:
			c = countX(sublist,sub)
			#print(sub, sub.parent)
			#print ("\n",sub,groupmap[sub.parent])
			#iii = [k for k,l in groupmap[sub.parent].items() if l == moose.element(sub)]
			for oo in groupmap[sub.parent]:
				#print(" k, l ",oo)
				if type(oo) is tuple:
					if oo[1] == moose.element(sub):
						newsub1 = oo[0]
			#print("|n ",newsub)
			subname = sub.name
			subgrp = findGroup_compt(sub)
			newsub = checkdigitEqu(startstringdigit,subgrp,sub)
			#print("sub ",reacsize,reacgrp, subname, "## ",newsub)
			if newsub1 in node_color:
				reaction_color = node_color[newsub1]
			else:
				reaction_color,sIndex = getColor(sIndex)
				node_color[newsub1] = reaction_color

			if c >1:
				#print("416 ",newsub,newsub1,reacsize)

				edgelist1 = edgelist1 +"\n"+newsub1+"->"+reacsize+"[arrowhead=normal weight = "+str(edge_weight)+ " arrowsize = "+str(edge_arrowsize)+" color=\""+reaction_color+":"+reaction_color+"\""+ " label=\" "+str(c)+"\" style=bold]"
			else:
				edgelist1 = edgelist1 +"\n"+newsub1+"->"+reacsize+"[arrowhead=normal weight = "+str(edge_weight)+ " arrowsize = "+str(edge_arrowsize)+" color=\""+reaction_color+":"+reaction_color+"\" style=bold]"
		for prd in prdlistU:
			c = countX(prdlist,prd)
			for oo in groupmap[prd.parent]:
				#print(" k, l ",oo)
				if type(oo) is tuple:
					if oo[1] == moose.element(prd):
						newprd1 = oo[0]
			prdname = prd.name
			prdgrp = findGroup_compt(prd)
			newprd = checkdigitEqu(startstringdigit,prdgrp,prd)
			if newprd1 in node_color:
				reaction_color = node_color[newprd1]
			else:
				reaction_color,sIndex = getColor(sIndex)
				node_color[newprd1] = reaction_color
			if c >1:
				edgelist1 = edgelist1 +"\n"+reacsize+"->"+newprd1+"[arrowhead=normal weight = "+str(edge_weight)+ " arrowsize = "+str(edge_arrowsize)+" color=\""+reaction_color+":"+reaction_color+"\""+ " label=\" "+str(c)+"\" style=bold]"
			else:
				edgelist1 = edgelist1 +"\n"+reacsize+"->"+newprd1+"[arrowhead=normal weight = "+str(edge_weight)+ " arrowsize = "+str(edge_arrowsize)+" color=\""+reaction_color+":"+reaction_color+"\" style=bold]"
	return(edgelist1,node_color)
		

def file_choices(choices,fname,iotype):
	ext = (os.path.splitext(fname)[1][1:]).lower()
	if iotype == "outputfile":
		if ext not in choices:
			parser.error("Requires output filetype {}".format(choices))
	else:
		if ext != "xml" and ext != "g":
			parser.error("Requires Genesis or SBML file format ")
			
	return fname

if __name__ == "__main__":

	parser = argparse.ArgumentParser( description = 'This program generates a reaction diagram for a SBML/Genesis model. It converts the specified Genesis/SBML, to the dot format. The dot file is further converted to an image in png/svg format\n')
	parser.add_argument('model',type=lambda s:file_choices((".xml",".g"),s,"input"), help='Required: filename of model, in Genesis format.')
	#parser.add_argument( 'model', type = str, help='Required: filename of model, in JSON format.')
	parser.add_argument( '-o', '--output', type=lambda out:file_choices(("png","svg"),out,"outputfile"), help='Optional: writes out the png model into named file. default takes json filename')
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
		if dirpath:
			outputfile = os.path.dirname(args.model)+"/"+os.path.splitext(os.path.basename(args.model))[0]+".png"	
		else:
			outputfile = os.path.splitext(args.model)[0]+".png"
	else:
		outputfile = args.output
	import moose
	#jsonDict = loadModel( args.model )
	#modelpath = parseModel( jsonDict )
	#modelpath = moose.loadModel(args.model,'/model')
	ext = (os.path.splitext(args.model)[1][1:]).lower()
	if ext == "xml":
		modelpath,errormsg = moose.readSBML(args.model,'/model' ,"ee")
	elif ext == "g":
		modelpath = moose.loadModel(args.model,'/model',"ee")
	else:
		print ("Input file should be genesis or SBML")
		exit()
	#print("#######",modelpath)
	jsontoPng(modelpath, outputfile, ranksep = args.ranksep, hasLegend = not args.no_legend, fontsize = args.fontsize, showGroups = not args.no_groups,specific_group = args.specific_group )
