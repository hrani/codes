#This file reads kkit/SBML file into moose and writes out to sql for dumping the model into doqcs database

from datetime import date
import os,sys
import moose

filepath = sys.argv[1]
print (filepath)
print(os.path.basename(filepath))
print(os.path.basename(filepath).split('.')[0]," -1 	",os.path.basename(filepath).split('.')[-1])
filename = os.path.basename(filepath).split('.')[0]
fileext = os.path.basename(filepath).split('.')[-1]

accessname = filename
if fileext == 'g':
	moose.loadModel(filepath,"/"+filename)
elif fileext == '.xml':
	moose.readSBML(filepath,"/"+filename)
else:
	print ("Provide GENESIS or SBML files")
f_doqcs = open(filename+".sql", "w")
f_doqcs.write("USE doqcs\n");
f_doqcs.write("INSERT INTO accession(\n");
f_doqcs.write( "\tname,accesstype,entrydate,transcriber,developer,\n");
f_doqcs.write( "\tspecies, tissue, cellcompartment, \n");
f_doqcs.write( "\tmethodology, source, model_implementation, \n");
f_doqcs.write( "\tmodel_validation, notes)\n");

f_doqcs.write( "VALUES(" );
f_doqcs.write( "\n\t\""+accessname+"\", \"Pathway\", Now(), \"Upinder S. Bhalla, NCBS\", \"Upinder S. Bhalla, NCBS\",");
f_doqcs.write( "\n\t\"\", \"\", \"\",");
f_doqcs.write( "\n\t\"\", \"\", \"\", ");
f_doqcs.write( "\n\t\"\",\"\");");

f_doqcs.write( "\nSELECT @accessno := last_insert_id();\n");

f_doqcs.write( "UPDATE accession SET figure = LOAD_FILE("+"'"+os.getcwd()+'/'+accessname+'.png'+"'"+") WHERE accessno = @accessno;\n");

f_doqcs.write( "\nINSERT into modelfiles(accessno, model, format, is_native)\n");
f_doqcs.write( "VALUES(@accessno, LOAD_FILE("+"'"+os.getcwd()+"/"+accessname+".g"+"'"+"), 'GENESIS', '1' );\n");

f_doqcs.write( "\nINSERT into modelfiles(accessno, model, format, is_native)\n");
f_doqcs.write( "VALUES(@accessno, LOAD_FILE("+"'"+os.getcwd()+"/"+accessname+".xml"+"'"+"), 'SBML', '0' );\n");

f_doqcs.write( "\nINSERT into modelfiles(accessno, model, format, is_native)\n");
f_doqcs.write( "VALUES(@accessno, LOAD_FILE("+"'"+os.getcwd()+"/"+accessname+".json"+"'"+"), 'JSON', '0' );\n");

f_doqcs.write("INSERT INTO pathway\n");
f_doqcs.write("	(accessno, name, notes)\n");
f_doqcs.write( "	VALUES(\n");
	
# if len(moose.wildcardFind("/"+filename+"/##[TYPE=Neutral]"))>0:
# 	print (moose.wildcardFind('/'+filename+"/##[TYPE=Neutral]"))
# 	f_doqcs.write( "@accessno, 'Shared_Object_"+accessname +","")");
# else:
f_doqcs.write( "@accessno,'" + accessname + "',\"\");");

f_doqcs.write("\nSELECT @pathwayno := last_insert_id();\n");
f_doqcs.write("UPDATE pathway SET figure = LOAD_FILE("+"'"+os.getcwd()+"/shared.png"+"'"+") WHERE accessno = @accessno AND pathwayno = @pathwayno;\n");
comptlist = []
for compt in moose.wildcardFind('/'+filename+"/##[ISA=ChemCompt]"):
	f_doqcs.write("INSERT INTO geometry(\n");
	f_doqcs.write("accessno, name, size, dimensions, shape, outside, organelle, notes)\n");
	f_doqcs.write(" VALUES(\n");
	f_doqcs.write("	@accessno,\""+compt.name+"\","+str(compt.volume)+",\""+str(compt.numDimensions)+"\",\"""sphere\",\"\",\""+compt.name+"\",\"\"); ");
	f_doqcs.write("\nSELECT @comptid0 := last_insert_id();");
	comptlist.append(compt)

for cmplist in comptlist:
	print(" c ",cmplist,moose.element(cmplist).path)
	sumtotallist = moose.wildcardFind(cmplist.path+"/##[ISA=Function]")
	for molecule in moose.wildcardFind(cmplist.path+"/##[ISA=PoolBase]"):
		print("molecule ",molecule.name, molecule.parent)
		if not molecule.parent.isA('Enz'):
			convertmillitomicro = molecule.concInit*1000
			f_doqcs.write("\nINSERT INTO molecule(name, accessno, pathwayno, is_buffered, initial_conc, issumtotal_available, D, molwt, Gid, notes)");
			f_doqcs.write("  VALUES('"+molecule.name+"', @accessno, @pathwayno,'");
			if molecule.isBuffered == False:
				f_doqcs.write('0'+"',"+str(convertmillitomicro)+",'");
			else:
				f_doqcs.write('1'+"',"+str(convertmillitomicro)+",'");
			if molecule.name in [x.parent.name for x in sumtotallist]:
				f_doqcs.write("1")
			else:
				f_doqcs.write("0")
			f_doqcs.write("', '0', 0, @comptid0,\" \");");
	
	for reaction in moose.wildcardFind(cmplist.path+"/##[ISA=Reac]"):
		print("reaction ",reaction.name,reaction.Kf,reaction.Kb)
		f_doqcs.write("\nINSERT INTO reaction(name, accessno, pathwayno, kf,kb)")
		sub = moose.element(reaction).neighbors["sub"]
		prd = moose.element(reaction).neighbors["prd"]
		no_sub = len(sub)
		no_prd = len(prd)
		Kf = reaction.Kf*pow(0.001,no_sub-1)
		f_doqcs.write("\nVALUES('"+reaction.name+"', @accessno, @pathwayno,"+ str(Kf)+","+str(reaction.Kb)+");")
		f_doqcs.write("\nSELECT @rid := last_insert_id();\n");

		for rsub in sub:
			print (reaction.name, rsub)
			f_doqcs.write("\nSELECT @RSid := id FROM molecule WHERE\n");
			f_doqcs.write("	name = '" +rsub.name + "' AND\n");
			f_doqcs.write("	accessno = @accessno AND pathwayno = @pathwayno;")
			f_doqcs.write("\nINSERT INTO moleculelist (mer_id, m_id, type)\n");
			f_doqcs.write("	VALUES(@rid, @RSid, '4');\n");
		for rprd in prd:
			print (reaction.name, rprd)
			f_doqcs.write("\nSELECT @RSid := id FROM molecule WHERE\n");
			f_doqcs.write("	name = '" +rprd.name + "' AND\n");
			f_doqcs.write("	accessno = @accessno AND pathwayno = @pathwayno;")
			f_doqcs.write("\nINSERT INTO moleculelist (mer_id, m_id, type)\n");
			f_doqcs.write("	VALUES(@rid, @RSid, '5');\n");

	for enzyme in moose.wildcardFind(cmplist.path+'/##[ISA=EnzBase]'):
		print("Enz ",enzyme)
		print( moose.showfields(enzyme))
		if enzyme.className == "MMenz":
			Km = enzyme.numKm
		else:
			k1 = enzyme.k1
			k2 = enzyme.k2
			k3 = enzyme.k3
			Km = ((k2 + k3)/(k1 * enzvol))*1000;

		f_doqcs.write("SELECT @pmid := id FROM molecule WHERE\n");
		f_doqcs.write("	name = '" + (enzyme.parent).name + "' AND\n");
		f_doqcs.write("	accessno = @accessno AND pathwayno = @pathwayno;");
		f_doqcs.write("\nINSERT INTO enzyme(");
		f_doqcs.write("name, accessno, pathwayno, parent_mid, km, vmax, ratio, is_available, notes)\n");
		f_doqcs.write("  VALUES	('" + enzyme.name+"',@accessno,@pathwayno,@pmid,");
		if enzyme.className == "enz":
			f_doqcs.write(str(Km) +", " + str(k3) + ", " + str(k2/k3) + ",");
		else:
			f_doqcs.write(str(Km) +", " + str(enzyme.kcat) + ", " + str('4') + ",");
		if enzyme.className == "MMenz":
			f_doqcs.write( "1,'"+ "');\n");
		else:
			f_doqcs.write( "0,'"+ "');\n");
		#f_doqcs.write(" " +");\n")
		f_doqcs.write("SELECT @eid := last_insert_id();\n");
		for esub in moose.element(enzyme).neighbors["sub"]:
			print (enzyme.name, rsub)
			f_doqcs.write("\nSELECT @ESid := id FROM molecule WHERE\n");
			f_doqcs.write("	name = '" +esub.name + "' AND\n");
			f_doqcs.write("	accessno = @accessno AND pathwayno = @pathwayno;")
			f_doqcs.write("\nINSERT INTO moleculelist (mer_id, m_id, type)\n");
			f_doqcs.write("	VALUES(@eid, @ESid, '2');\n");
		for eprd in moose.element(enzyme).neighbors["prd"]:
			print (enzyme.name, rprd)
			f_doqcs.write("\nSELECT @ESid := id FROM molecule WHERE\n");
			f_doqcs.write("	name = '" +eprd.name + "' AND\n");
			f_doqcs.write("	accessno = @accessno AND pathwayno = @pathwayno;")
			f_doqcs.write("\nINSERT INTO moleculelist (mer_id, m_id, type)\n");
			f_doqcs.write("	VALUES(@eid, @ESid, '3');\n");
