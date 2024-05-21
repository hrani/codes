''' Convert json file fro new schema Jan30,Feb1
    1.Readout->entities and Readout->normalization->entities array replaced to string
    2.Experiment->BarChart ->sampling =start, then settleTime at Readout level is retained
    3.DoseResponse -> sampling=dose, then "dose" and settleTime is retained
    4.TimeSeries -> sampling=presetTime, 'time' field ->replaced 'samptime' (settleTime will be discarded for TS and DP)
       b.TimeSeries -> sampling=[start,end,min,max,each] 'settleTime',"time" will be removed
    5.DoseResponse -> dose then “time” from Readout->normalization will be deleted and settleTime will be kept at Readout level
        B. same way doseResponse->start,max then what is the significant for settleTime? Time will be removed
    "Feb1"
    6.According to /home/harsharani/Desktop/codes/HTpipeline_jan30/newMappingScheme.txt replaced old name to new names in
        -stimuli entity
        -Readout parameter entities
        -Readout enitities
        -Readout normlaization entity
        -Modification parameter change entity
        -Modification subset
    "Feb5 "

    7."entity": {"type": "string" }, "details": {"type": "string" }, "identifier-type": {"enum": ["UnitProt","ChEMBL"] },"identifier": {"type":"string"}
        Now  entity/entities will become a object now
        entity/entities
        { "name": {"type": "string" }, "alias": {"type": "string" },"identifier-type": {"enum": ["UnitProt","ChEMBL"] },"identifier": {"type":"string"},"notes": {"type": "string" }}
    8. Added FileType and Version and Made changes name and alias, if name and alias is defined, then later part of the file will contains only name to
      reduce the size of file 
'''
import os
import json
import argparse
from collections import OrderedDict
import subprocess
import re

def modify_entities(json_content):
    testing = []
    start_index = json_content.find('"entity"')
    while start_index != -1:
        # Find end index of "field" after the start index
        end_index = json_content.find('"field"', start_index)
        if end_index != -1:
            # Find the space index after the end index of "field"
            #print("40 ",start_index)
            # Extract and modify the entity string
            entity_str = str(json_content[start_index:end_index])# + '\n' + ' ' * (space_index - end_index))
            #print("@@",entity_str)
            # Replace the entity in the JSON content
            json_content = json_content[:start_index] + entity_str + json_content[end_index:]
            testing.append([start_index,end_index,entity_str])
            #print("54 ",json_content)
            # Update the start index for the next iteration
            start_index = json_content.find('"entity"', end_index)
            
        else:
            # If "field" is not found after the start index, break the loop
            break
    #print("56 ",testing)
    return testing,json_content

def findNewname(entity_dict,data):
    eoI = {key if key != 'entity' else 'name': value for key, value in entity_dict.items()}
    if isinstance(eoI["name"],list):
        eoI_name = eoI['name'][0]
        alias = eoI['name'][0]
    if isinstance(eoI["name"],str):
        eoI_name = eoI['name']
        alias = eoI['name']

    if eoI_name not in old_new_mapdict:
        notfound.append(eoI_name)
        name = eoI_name

    if eoI_name in old_new_mapdict and eoI_name != old_new_mapdict[eoI_name]:
        eoI_name = ''.join(old_new_mapdict[eoI_name])
        name = eoI_name                            
    else:
        name = eoI_name
    return name,alias
def main():
    ''' 
    This program loads a map file and then goes through a directory
    of experiments to run on the model(s) associated with the map file.
    It picks out those experiments whose inputs and outputs are present
    in the map file.
    '''
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument( "-e", "--exptDir", type = str, 
        help = "Optional: Directory in which to look for all experiments to convert.", default = "" )
    parser.add_argument( "-f", "--exptfile", type = str, 
        help = "Optional: File to convert.", default = "" )
    parser.add_argument( "-o", "--outputDir", type = str, 
        help = "Optional: Directory in output file need to save.", default = "" )
    parser.add_argument('-a', '--alias',type=str,help="Optional: True alias will be written, False then will removed",default="True")
    args = parser.parse_args()
    listOfFiles =[]
    #aliasneeded = args.alias.lower()
    aliasneeded = args.alias.lower() == "true"

    if args.exptDir == "" and args.exptfile == "":
        print("Directory and filename is null")
    elif args.exptDir != "" and args.exptfile != "":
        print(" Pass either Directory or file ")
    else:
        if args.exptDir != "":
            if os.path.isdir(args.exptDir):
                exptDir = args.exptDir
                listOfFiles = [f for f in os.listdir(args.exptDir) if os.path.isfile(os.path.join(args.exptDir, f) )]
                listOfFiles = sorted(listOfFiles)
            else:
                print(args.exptDir ," is not directory")
        elif args.exptfile != "":
            if os.path.isfile(args.exptfile):
                exptDir = os.path.dirname(args.exptfile)
                
                listOfFiles =  [os.path.basename(args.exptfile)]
            else:
                print(args.exptfile," is not filename")
    global old_new_mapdict
    with open('/home/harsharani/AutSim/development_branch_FS_HT/Nisha_htpipeline_Mar112024/findsim_model_map4mNisha/newMappingScheme.txt','r') as file:
        old_new_mapdict = {}
        for line in file:
            key, value = line.strip().split(':')
            key = key.strip().strip('"')
            value = value.strip().strip('\"').strip(',')
            value = value.replace('"', '')
            old_new_mapdict[key.strip()] = value.strip()

    TS_readout_noNormalization = []
    global notfound
    notfound = []
    if listOfFiles:
        for filename in listOfFiles:
            if filename.endswith(".json") and filename not in ["Antion2008_Fig5B_pS6_DHPG_0126.json","Antion2008_Fig5B_pS6_DHPG_u0126.json","Antion2008_Fig4B_pS6_DHPG_Rapamycin.json","Antion2008_Fig5B_pS6_DHPG_UO126.json"]: #["Abe2000_Fig1A_n.json","Aditi_S6K_n.json","Akam1997_n.json","Alessi1998_Fig4B_n.json","Antion2008_Fig1B_pERK2_n.json","Antion2008_Fig3B_DHPG_n.json","Abe2000_Fig1B_n.json","Akam1997_Fig1_inset_n.json","Albasanz2005_Fig1_n.json","Alessi1998_Fig6A_n.json","Antion2008_Fig1B_pS6_n.json","Antion2008_Fig3B_DHPG_Wortmannin_n.json","Aditi_AKT_n.json","Akam1997_Fig3_DHPG_n.json","Alessi1997_Fig3A_n.json","Alessi1998_n.json","Antion2008_Fig1C_n.json","Antion2008_Fig4B_pS6_DHPG_n.json","Aditi_ERK_n.json","Akam1997_Fig3_Glu_n.json","Alessi1997_Fig3B_n.json","Antion2008_Fig1B_mTOR_n.json","Antion2008_Fig3B_DHPG_LY294002_n.json"]:
                try:
                    exptDir
                    file_path = os.path.join( exptDir, filename )
                except:
                    file_path = filename
                if os.path.exists(file_path):

                    with open( file_path ) as fd:
                        headdata ={}
                        pop_list = []
                        namedefined=[]
                        findsim = json.load( fd )
                        print("---------- File ",os.path.splitext(filename)[0]+".json","",findsim["Experiment"]["design"],"--------------------------")
                        if "FileType" not in findsim:
                            headdata["FileType"] = "Findsim"
                        
                        if "Version" in findsim:
                            findsim["Version"]="2.0"
                        else:
                            print("\tFindsim Version missing ")
                            headdata["Version"] = "2.0"
                        
                        if "Stimuli" in findsim:
                            stimuliIndex = 0
                            for soI_dict in findsim["Stimuli"]:
                                new_soI = {}
                                name,alias = findNewname(soI_dict,"Stimuli")
                                for k,v in soI_dict.items():
                                    if k in ["timeUnits",'quantityUnits']:
                                        new_soI.update({k:v})
                                    elif k == 'entity':
                                        if name not in namedefined:
                                            namedefined.append(name)
                                            if aliasneeded:
                                                new_soI.update({"entity": {"name": name, "alias": alias}})
                                            else:
                                                new_soI.update({"entity": {"name": name}})
                                        else:
                                            new_soI.update({"entity": {"name": name}})
                                    else:
                                        new_soI.update({k:v})
                                findsim["Stimuli"][stimuliIndex] = dict(new_soI)
                                stimuliIndex+=1

                        # ---------Readouts Section---------
                        if findsim["Experiment"]["design"] == "BarChart":
                            if "bardata" in findsim["Readouts"]:
                                for bl in findsim["Readouts"]["bardata"]:
                                    for k,v in bl.items():
                                        if k == 'stimulus':
                                            if len(v)>0:
                                                for vv in v:
                                                    if vv in old_new_mapdict and vv != old_new_mapdict[vv]:
                                                        index_to_insert_alias = next((index for index, item in enumerate(v) if item ==vv),None)
                                                        v[index_to_insert_alias] = old_new_mapdict[vv]
                        if findsim["Experiment"]["design"] == "DirectParameter":
                            index_to_paramdata = 0
                            for l in findsim["Readouts"]["paramdata"]:
                                name,alias = findNewname(l,"Readouts")
                                newRdout = {}    
                                for k,v in l.items():
                                    if k == 'entity':
                                        if name not in namedefined:
                                            namedefined.append(name)
                                            if aliasneeded:
                                                newRdout.update({"entity": {"name": name, "alias": alias}})
                                            else:
                                                newRdout.update({"entity": {"name": name}})
                                        else:
                                            newRdout.update({"entity": {"name": name}})
                                    else:
                                        newRdout.update({k:v })
                                findsim["Readouts"]["paramdata"][index_to_paramdata] = dict(newRdout)
                                index_to_paramdata +=1
                        else:
                            findsim["Readouts"] = {key if key != 'entities' else 'entity': value for key, value in findsim["Readouts"].items()}
                            readout_new = {}
                            for roI,roI_value in findsim["Readouts"].items():
                                if roI == "data":
                                    rdata = findsim["Readouts"]
                                    for i in range(len(rdata['data'])):
                                        rdata['data'][i] = '[' + ', '.join(map(str, rdata['data'][i])) + ']'
                                if roI == "entity":
                                    if isinstance(findsim["Readouts"]["entity"],list):
                                        Rd_entity = findsim["Readouts"]["entity"][0]
                                    else:
                                        Rd_entity = findsim["Readouts"]["entity"]
                                    name,alias = findNewname(findsim["Readouts"],"Readout")
                                    
                                    if name not in namedefined:
                                        namedefined.append(name)
                                        if aliasneeded:
                                            readout_new.update({"name": name, "alias": alias})
                                        else:
                                            readout_new.update({"name": name})
                                    else:
                                        readout_new.update({"name": name})
                                    
                                if roI == "normalization":
                                    findsim["Readouts"]["normalization"] = {key if key != 'entities' else 'name': value for key, value in findsim["Readouts"]["normalization"].items()}
                                    alias = findsim["Readouts"]["normalization"]["name"][0]
                                    details_fieldR_nor = ('alias',findsim["Readouts"]["normalization"]["name"][0]+"}")

                                    if findsim["Readouts"]["normalization"]["name"][0] not in old_new_mapdict:
                                        print("\t\t Readouts->normalization->name ",findsim["Readouts"]["normalization"]["name"][0]," not found in the map file program exited")
                                        notfound.append(findsim["Readouts"]["normalization"]["name"][0])
                                        exit()
                                    if findsim["Readouts"]["normalization"]["name"][0] in old_new_mapdict and findsim["Readouts"]["normalization"]["name"] != old_new_mapdict[findsim["Readouts"]["normalization"]["name"][0]]:
                                        findsim["Readouts"]["normalization"]["name"] = ''.join(old_new_mapdict[findsim["Readouts"]["normalization"]["name"][0]])

                                    else:
                                        findsim["Readouts"]["normalization"]["name"] = findsim["Readouts"]["normalization"]["name"][0]
                                    nor_value = findsim["Readouts"]["normalization"]
                                    new_nor = {}
                                    for k,v in nor_value.items():
                                        if k != "name":
                                            new_nor.update({k:v})
                                        else:
                                            if v not in namedefined:
                                                new_nor.update({ "entity": {'name': findsim['Readouts']['normalization']['name'],'alias':alias} })
                                            else: 
                                                new_nor.update({ "entity" : {'name':findsim["Readouts"]["normalization"]['name']} })
                                    findsim["Readouts"]["normalization"] = new_nor
                            findsim["Readouts"]["entity"] = readout_new

                        if findsim["Experiment"]["design"] in ["TimeSeries"]:
                                if "normalization" in findsim["Readouts"]:
                                    if findsim["Readouts"]["normalization"]["sampling"] == "presetTime":
                                        if "time" in findsim["Readouts"]["normalization"]:
                                            findsim["Readouts"]["normalization"]["sampTime"] = findsim["Readouts"]["normalization"].pop("time")
                                        if "settleTime" in findsim["Readouts"]:
                                            #findsim["Readouts"].pop("settleTime")
                                            if "settleTime" not in pop_list:
                                                pop_list.append("settleTime")
                                    if findsim["Readouts"]["normalization"]["sampling"] in ["start","end","min","max","each"]:
                                        if "dose" in findsim["Readouts"]["normalization"]:
                                            pop_list.append("dose")
                                        if "time" in findsim["Readouts"]["normalization"]:
                                            if "time" not in pop_list:
                                                pop_list.append("time")
                                        if "settleTime" in findsim["Readouts"]:
                                            #findsim["Readouts"].pop("settleTime")
                                            if "settleTime" not in pop_list:
                                                pop_list.append("settleTime")    
                                else:
                                    if filename not in TS_readout_noNormalization: 
                                        TS_readout_noNormalization.append(filename)
            
                        elif findsim["Experiment"]["design"] in ["BarChart","DoseResponse"]:
                                if "normalization" in findsim["Readouts"]:
                                    if 'dose' in findsim["Readouts"]["normalization"]:
                                        if "time" in findsim["Readouts"]["normalization"]:
                                            if 'time' not in pop_list:
                                                pop_list.append("time")
                                    if findsim["Readouts"]["normalization"]["sampling"] in ["start","end","min","max","each"]:
                                        if "time" in findsim["Readouts"]["normalization"]:
                                            if "time" not in pop_list:
                                                pop_list.append("time")            

                        for i in pop_list:
                            if i == "settleTime":
                                findsim["Readouts"].pop(i)
                            if i in ["time","dose"]:
                                findsim["Readouts"]["normalization"].pop(i)
                        
                        if "Modifications" in findsim:
                                ssStr = []
                                if "subset" in findsim["Modifications"]:
                                    if len(findsim["Modifications"]["subset"]) == 1 and findsim["Modifications"]["subset"][0] == 'all':
                                        ssStr.append({"name":"all"})
                                    else:
                                        for index,i in enumerate(findsim["Modifications"]["subset"]):
                                            if i != "all":
                                                #print("not all i ",i)
                                                if i not in old_new_mapdict:
                                                    print("\t\t Modification->subset ",i,"not found in the map file exited")
                                                    notfound.append(i)
                                                    exit()
                                                if i in old_new_mapdict and i != old_new_mapdict[i]:
                                                    name = old_new_mapdict[i]
                                                else:
                                                    name = i
                                                #for now subset is filling with name only
                                                if name not in namedefined:
                                                    if name != i:
                                                        ssStr.append({"name":name,"alias":i})
                                                    else:
                                                        ssStr.append({"name":name})
                                                else:
                                                     ssStr.append({"name":name})
                                                
                                    findsim["Modifications"]["subset"] = ssStr
                                
                                i2Dstr = []
                                itemstodelete_notfound=[]
                                if "itemsToDelete" in findsim["Modifications"]:
                                    for index,i2D in enumerate(findsim["Modifications"]["itemsToDelete"]):
                                        if i2D not in old_new_mapdict:
                                            notfound.append(i2D)
                                            itemstodelete_notfound.append(i2D)
                                        if i2D in old_new_mapdict and i2D != old_new_mapdict[i2D]:
                                            name = old_new_mapdict[i2D] 
                                        else:
                                            name = i2D

                                        if name not in namedefined:
                                            if name != i2D:
                                                i2Dstr.append({"name":name,"alias":i2D})
                                            else:
                                                i2Dstr.append({"name":name})
                                        else:
                                            i2Dstr.append({"name":name})
                                    
                                    findsim["Modifications"]["itemsToDelete"] = i2Dstr

                                if "parameterChange" in findsim["Modifications"]:
                                    paraChange = []
                                    for pc in findsim["Modifications"]["parameterChange"]:
                                        pc = {key if key != 'entity' else 'name': value for key, value in pc.items()}
                                        new_pc = {}
                                        
                                        if pc["name"] not in old_new_mapdict:
                                            print("\t\t Modification->ParamterChange->name ",pc["name"]," not found in the map file program exited")
                                            notfound.append(pc["name"])
                                            exit()
                                        
                                        if pc["name"] in old_new_mapdict and pc["name"] != old_new_mapdict[pc["name"]]:
                                            name = old_new_mapdict[pc['name']]
                                        
                                        elif pc["name"] in old_new_mapdict and pc["name"] == old_new_mapdict[pc["name"]]:
                                            name = old_new_mapdict[pc['name']]
                                        
                                        if pc['name'] not in namedefined:
                                            namedefined.append(pc['name'])
                                            new_pc.update({"entity": { "name": name, 
                                                        "alias":pc['name']}
                                                        })
                                        else:
                                            new_pc.update({"entity": { "name": name}})
                                        new_pc.update({
                                            'field': pc['field'],
                                            'value': pc['value'],
                                            'units': pc['units']
                                            })
                                        paraChange.append(new_pc)
                                    findsim["Modifications"]["parameterChange"] = paraChange
                                    

                outputdirectory = args.outputDir
                outputfilename = outputdirectory+os.path.splitext(filename)[0]+"_n.json" 
                print("outputdirectory ",outputfilename)
                #del outputfilename

                with open(outputfilename, 'w') as fout:
                    if headdata != "":
                        merged_json = {**headdata, **findsim}
                    else:
                        merged_json = {**findsim}
                    json_dumps_str = json.dumps(merged_json)#, indent=4)
                    
                    # start_index = json_dumps_str.find('"entity"') + len('"entity"')
                    # end_index = json_dumps_str.find('"field"')
                    # space_index = json_dumps_str.find(' ',end_index)
                    # entity_str = json_dumps_str[start_index:space_index-(space_index-end_index)].replace('\n', '').replace(' ', '')+"\n"+" "*(space_index-end_index+4)
                    
                    #modified_json_content = modify_entities(json_dumps_str)
                    #print("394 ",modified_json_content)
                    #modified_json_content = json_dumps_str[:start_index] + entity_str + json_dumps_str[end_index:]
                    f ,j = modify_entities(json_dumps_str)
                    #print("j ",type(j))
                    j_dict = json.loads(j)
                    json.dump(j_dict, fout, indent=4) 
                    #print(j, file=fout, indent=4)
                    #indented_json_string1 = j.replace('\n                ]',  ']').replace('\n                    ','')
                    #indented_json_string2 = indented_json_string1.replace('\"[',  '[').replace(']\"' ,']')
                    #print(json.dumps(indented_json_string2,indent=4), file=fout)
                    #print(json.dumps(j, indent=4), file=fout)
                    #print(f,type(f))
                    i = 1
                    '''
                    modified_json_content =""
                    for ff in f:
                        entity_str = ff[2].replace('\n', '').replace(' ', '')#+"\n"+" "*(4)
                        if i == 1:
                            modified_json_content += json_dumps_str[:ff[0]-1] + entity_str
                            previos_end = ff[1]
                            i = 0 
                        else:
                            modified_json_content += json_dumps_str[previos_end:ff[0]-1]+entity_str
                            previos_end = ff[1]

                        print("\n\n! ",modified_json_content)
                        #+ json_dumps_str[ff[1]:]
                
                    indented_json_string1 = modified_json_content.replace('\n                ]',  ']').replace('\n                    ','')
                    indented_json_string2 = indented_json_string1.replace('\"[',  '[').replace(']\"' ,']')
                    print(indented_json_string2, file=fout)
                    
                    #file_path = outputdirectory+os.path.splitext(filename)[0]+"_n.json"
                    '''
                # with open(file_path, 'r') as file:
                #     lines = file.readlines()

                # # # Perform replacements
                
                # substrings_to_match = ["\"name\"", "\"alias\""]
                # indentation = 0
                # for i, line in enumerate(lines):
                # #     # if '"{\'' in line:
                # #     #     lines[i] = line.replace('\"{\'',"{\'").replace('}\"','}').replace('\'','\"')
                #     if '},{' in line:
                #         print("364 ",i,line)
                #         #line[i] = line.replace('},\n{','}')
                #         lines[i] = line.replace('},{','},\n             {').replace('\"{',"{").replace('},\"','}').replace('\'','\"')
                    
                # # # Write the modified content back to the file
                # with open(file_path, 'w') as file:
                #     file.writelines(lines)            
        
                # print("File converted and saved at",outputdirectory+os.path.splitext(filename)[0]+"_n.json")
                if notfound:
                    print("Not Found ",str(set(notfound)))
                if itemstodelete_notfound:
                    print("These object in itemstodelete section not found please check ",itemstodelete_notfound)
                #return file_path
if __name__ == "__main__":
    main()
    #file_path = outputdirectory+os.path.splitext(filename)[0]+"_n.json"
   
