''' Convert json file fro new schema Nov20 
    1.Readout->entities and Readout->normalization->entities array replaced to string
    2.Experiment->BarChart ->sampling =start, then settleTime at Readout level is retained
    3.DoseResponse -> sampling=dose, then "dose" and settleTime is retained
    4.TimeSeries -> sampling=presetTime, 'time' field ->replaced 'samptime' (settleTime will be discarded for TS and DP)
       b.TimeSeries -> sampling=[start,end,min,max,each] 'settleTime',"time" will be removed
    5.DoseResponse -> dose then “time” from Readout->normalization will be deleted and settleTime will be kept at Readout level
        B. same way doseResponse->start,max then what is the significant for settleTime? Time will be removed
''' 

import os
import json
import argparse


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
    
    args = parser.parse_args()
    listOfFiles =[]
    
    if args.exptDir == "" and args.exptfile == "":
        print("Directory and filename is null")
    elif args.exptDir != "" and args.exptfile != "":
        print(" Pass either Directory or file ")
    else:
        if args.exptDir != "":
            if os.path.isdir(args.exptDir):
                exptDir = args.exptDir
                listOfFiles = [f for f in os.listdir(args.exptDir) if os.path.isfile(os.path.join(args.exptDir, f) )]
            else:
                print(args.exptDir ," is not directory")
        elif args.exptfile != "":
            if os.path.isfile(args.exptfile):
                exptDir = os.path.dirname(args.exptfile)
                
                listOfFiles =  [os.path.basename(args.exptfile)]
            else:
                print(args.exptfile," is not filename")

    TS_readout_noNormalization = []

    if listOfFiles:
        for filename in listOfFiles:
            if filename.endswith(".json"):
                try:
                    exptDir
                    file_path = os.path.join( exptDir, filename )
                except:
                    file_path = filename
                if os.path.exists(file_path):
                    with open( file_path ) as fd:

                        pop_list = []
                        findsim = json.load( fd )
                        if "Version" in findsim:
                            findsim["Version"]="2.0"
                        else:
                            print("Findsim Version missing ")
                        
                        for roI in findsim["Readouts"]:
                            if roI == "entities":
                                findsim["Readouts"]["entities"] = ''.join(findsim["Readouts"]["entities"])
                            if roI == "normalization":
                                findsim["Readouts"]["normalization"]["entities"] = ''.join(findsim["Readouts"]["normalization"]["entities"])    
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
                            if i == "time":
                                findsim["Readouts"]["normalization"].pop(i)
                outputdirectory = args.outputDir
                with open(outputdirectory+os.path.splitext(filename)[0]+"_n.json", 'w') as fout:
                    json_dumps_str = json.dumps(findsim, indent=4)
                    print(json_dumps_str, file=fout)
                    print("File converted and saved at",outputdirectory+os.path.splitext(filename)[0]+"_n.json")
if __name__ == "__main__":
    main()
