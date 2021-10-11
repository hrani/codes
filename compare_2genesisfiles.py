import moose
#from setsolver import *
from moose.SBML import *
import numpy

moose.loadModel('/tmp/fsDev_MC_b4c.g','/MC')
moose.loadModel('/tmp/fsDev_ZG_b4c.g','/ZG')

print ("\nCube Mesh")
for c in moose.wildcardFind('/MC/##[ISA=CubeMesh]'):
	p2str = c.path.replace('/MC','/ZG')
	print (c.className,c.name," ",c.volume, " Second ",moose.element(p2str), moose.element(p2str).volume)

pdoesnotexist, pclass, pConcInit,pninit = "", " ", " ", " "

for p1 in moose.wildcardFind('/MC/##[ISA=PoolBase]'):
	p2str = p1.path.replace('/MC','/ZG')
	if moose.exists(p2str) == False:
		pdoesnotexist = pdoesnotexist+"\n"+p2str
	else:
		p2 = moose.element(p2str)
		if p1.className != p2.className:
			pclass = pclass +"\n"+ p1.path+" "+ p1.className + " p2:  "+ p2.path + p2.className
		if not numpy.allclose(p1.concInit, p2.concInit):
			pConcInit = pConcInit+"\n"+ p1.path + " :"+str(p1.concInit) +" p2: "+p2.path + ": "+str(p2.concInit)

if pdoesnotexist:
	print (" # pool doesn't exist ",pdoesnotexist)
if pclass != " ":
	print (" \t #pools have differenct classType ", pclass)
if pConcInit:
	print (" \t \t #pools have differenct coninit", pConcInit )


reacdoesntexist = ""
ratediff = " "
subprddiff = ""
reacdoesntexist = ""
for r1 in moose.wildcardFind('/MC/##[ISA=ReacBase]'):
	r2str = r1.path.replace('/MC','/ZG')
	
	if moose.exists(r2str) == False:
		reacdoesntexist = reacdoesntexist+ "\n"+r2str
	else:
		r2 = moose.element(r2str)
		if(len(r1.neighbors['sub']) != len(r2.neighbors['sub']) or len(r1.neighbors['prd']) != len(r2.neighbors['prd']) ):
			subprddiff = subprddiff + "\n"+  r1, "sub " ,len(r1.neighbors['sub'])," prd", len(r1.neighbors['prd']), " r2 ",r2, "sub " ,len(r2.neighbors['sub'])," prd", len(r2.neighbors['prd'])
		#print " numpy ", numpy.around(r1.Kf, decimals=2), " ",numpy.around(r2.Kf, decimals=2)		
		#if not numpy.allclose(r1.Kf, r2.Kf) or not numpy.allclose(r1.Kb, r2.Kb):
		if (not numpy.allclose(r1.Kf,r2.Kf) or not numpy.allclose(r1.Kb,r2.Kb)):
			ratediff = ratediff + "\n"+  r1.parent.name+'/'+r1.name+ " Kf " +str(r1.Kf)+ " Kb "+ str(r1.Kb)+ " r2 "+r2.parent.name+'/'+r2.name+ " Kf " +str(r2.Kf)+" Kb "+ str(r2.Kb)

#exit()
if reacdoesntexist != " ":
	print ("Reaction doesn't exist ", reacdoesntexist)
if subprddiff != "":
	print ("sub or prd is differebt", subprddiff)
if ratediff != "":

	print ("rate are different ",ratediff)

enzdoesntexist = ""
ratediff = " "
subprddiff = ""
enzdoesntexist = ""
for e1 in moose.wildcardFind('/MC/##[ISA=EnzBase]'):
	e2str = e1.path.replace('/MC','/ZG')
	
	if moose.exists(e2str) == False:
		enzdoesntexist = enzdoesntexist+ "\n"+e2str
	else:
		e2 = moose.element(e2str)
		if(len(e1.neighbors['sub']) != len(e2.neighbors['sub']) or len(e1.neighbors['prd']) != len(e2.neighbors['prd']) ):
			subprddiff = subprddiff + "\n"+  e1, "sub " ,len(e1.neighbors['sub'])," prd", len(e1.neighbors['prd']), " e2 ",e2, "sub " ,len(e2.neighbors['sub'])," prd", len(e2.neighbors['prd'])
		
		#if not numpy.allclose(r1.Kf, r2.Kf) or not numpy.allclose(r1.Kb, r2.Kb):
		if not numpy.allclose(numpy.around(e1.Km, decimals=5),numpy.around(e2.Km, decimals=5)) or not numpy.allclose(numpy.around(e1.kcat, decimals=5),numpy.around(e2.kcat, decimals=5)):
			ratediff = ratediff + "\n"+  e1.parent.name+'/'+e1.name+ " Km " +str(e1.Km)+ " kcat "+ str(e1.kcat)+" k1 "+ str(e1.k1)+" k2 "+ str(e1.k2)+" k3 "+ str(e1.k3)+ " e2 "+e2.parent.name+'/'+e2.name+ " Km " +str(e2.Km)+" kcat "+ str(e2.kcat)+" k1 "+ str(e2.k1)+" k2 "+ str(e2.k2)+" k3 "+ str(e2.k3)


if enzdoesntexist != " ":
	print ("\n \t Enz doesn't exist ", reacdoesntexist)
if subprddiff != "":
	print ("\n \t \t sub or prd is differebt", subprddiff)
if ratediff != "":
	print ("\n \t \t rate are different ",ratediff)
# '''
# print "\nReacBase"
# for r1 in moose.wildcardFind('/MC/##[ISA=ReacBase]'):
# 	print r1.className,r1.name," ",r1.Kf," ",r1.Kb

# print "EnzBase"
# for e1 in moose.wildcardFind('/MC/##[ISA=EnzBase]'):
# 	print e1.className,e1.name," ",e1.kcat," ",e1.Km

# print "\nCube Mesh"
# for c in moose.wildcardFind('/ZG/##[ISA=CubeMesh]'):
# 	print c.className,c.name," ",c.volume
# print "\nPoolBase"
# for p2 in moose.wildcardFind('/ZG/##[ISA=PoolBase]'):
# 	print p2.className,p2.name," ",p2.concInit," ",p2.nInit
# print "\nReacBase"
# for r2 in moose.wildcardFind('/ZG/##[ISA=ReacBase]'):
# 	print r2.className,r2.name," ",r2.Kf," ",r2.Kb
# print "\nEnzBase "
# for e2 in moose.wildcardFind('/ZG/##[ISA=EnzBase]'):
# 	print e2.className,e2.name," ",e2.kcat," ",e2.Km

# '''
# '''
# moose.loadModel('/home/harsha/genesis_files/gfile/anno/Anno_acc50.g','/gen',"gsl")
# moose.reinit()
# moose.start(400)
# for i in range(10,19):
# 	print i,moose.element('/clock').dts[i]
# gp = moose.wildcardFind('/gen/##[ISA=PoolBase]')
# print " path conc n"
# for gps in gp:
# 	print gps.path,gps.concInit,gps.nInit, " ",gps.conc,gps.n

# mooseWriteSBML('/gen','/home/harsha/Trash/acc50_oct19cmd.xml')

# mooseReadSBML('/home/harsha/Trash/acc50_oct17.xml','/xml',"gsl")
# addSolver('/xml',"gsl")
# moose.reinit()
# moose.start(400)
# for j in range(10,19):
# 	print j,moose.element('/clock').dts[j]
# #function
# # xf = moose.wildcardFind('/gen/##[ISA=Function]')
# # print " genesis function ",moose.wildcardFind('/gen/##[ISA=Function]')
# # for xpf in xf:
# # 	print xpf,xpf.expr

# # xfx = moose.wildcardFind('/xml/##[ISA=Function]')
# # print " xml function",moose.wildcardFind('/xml/##[ISA=Function]')
# # for xpfx in xfx:
# # 	print xpfx,xpfx.expr
# '''
# #Reac
# '''
# print " kf kb numkf numkb"
# gr = moose.wildcardFind('/gen/##[ISA=ReacBase]')
# for grs in gr:
# 	print grs.name,grs.Kf," ",grs.Kb," ",grs.numKf," ",grs.numKb
# print " ----- "
# xr = moose.wildcardFind('/gen/##[ISA=ReacBase]')
# for xrs in xr:
# 	print xrs.name,xrs.Kf," ",xrs.Kb," ",xrs.numKf," ",xrs.numKb
# 	'''
# #Enz
# '''
# ge = moose.wildcardFind('/gen/##[ISA=EnzBase]')
# print " className name kcat Km numKm k1 k2 k3"
# for xps in ge:
# 	if xps.className == "ZombieEnz":
# 		print xps.className,xps.name," " ,xps.kcat," " ,xps.Km," " ,xps.numKm," " ,xps.k1," " ,xps.k2," " ,xps.k3
# 	else:
# 		print xps.className,xps.name," " ,xps.kcat," " ,xps.Km," " ,xps.numKm
# print "--------"
# xe = moose.wildcardFind('/xml/##[ISA=EnzBase]')
# for xps in xe:
# 	if xps.className == "ZombieEnz":
# 		print xps.className,xps.name," " ,xps.kcat," " ,xps.Km," " ,xps.numKm," " ,xps.k1," " ,xps.k2," " ,xps.k3
# 	else:
# 		print xps.className,xps.name," " ,xps.kcat," " ,xps.Km," " ,xps.numKm

# '''
# '''
# print " --- xml --"
# xp = moose.wildcardFind('/xml/##[ISA=PoolBase]')
# for xps in xp:
# 	print xps.path,xps.concInit,xps.nInit, " ",xps.conc,xps.n
# '''
