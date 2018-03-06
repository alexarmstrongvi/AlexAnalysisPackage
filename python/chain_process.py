import glob
import sys
import shutil
import os

import ROOT as r
import global_variables as g

directory_with_inputs = g.input_files
directory_with_files =  g.analysis_run_dir+"ntuples/"
output_directory =      g.analysis_run_dir+"ntuples/"

missing_dsid = open(g.missing_dsids_file,'w')

def get_files(dsid_list, directory) :
    """
    Grab all files in a given directory that have a desired DSID.
    Write all missing file DSIDs to output file missing_dsid 

    Args:
        dsid_list ( list(int) ): list of desired DSID  
        directory: directory containing files
    
    Returns:
        list(str): list of all file names that have a desired DSID
    """
    mc_files = glob.glob(directory + "mc/" + "CENTRAL*.root")
    data_files = glob.glob(directory + "data/" + "CENTRAL*.root")
    all_files = mc_files + data_files
    process_files = []
    for dsid in dsid_list :
        iii = 0
        for root_file in all_files :
            if str(dsid) in root_file :
                process_files.append(root_file)
                break # move to next dsid 
            iii += 1
        if iii == len(all_files): 
            missing_dsid.write(str(dsid)+"\n")
    return process_files

def main() :

    for group, dsids in g.groups.iteritems():
        print "Running over " + group
        print "\tFound %d dataset ids"%len(dsids)
        
        #missing_dsid.write("\n"+"="*10+" Sample:  "+group+"  "+"="*10+"\n\n")
        root_files = get_files(dsids, directory_with_files)
        print "\tFound %d root files for the process"%len(root_files)

        if len(root_files) != len(dsids) :
            print "\tERROR Did not find expected number of files!"
            #sys.exit()
        chain = r.TChain("superNt","superNt")
        for root_f in root_files :
            chain.Add(root_f)
       
        # Create TChain file for each group
        #output_file = r.TFile("%s%s.root"%(output_directory,group),"RECREATE")
        #chain.SetName("%s"%group)
        #chain.Write()
        ##chain.CloneTree(-1,"fast");
        #output_file.Write() 
        print "\tProcess tree has %d total entries"%chain.GetEntries()

if __name__ == '__main__':
    main()
