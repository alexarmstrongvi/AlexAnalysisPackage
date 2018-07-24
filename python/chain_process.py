import glob
import sys
import shutil
import os

import ROOT as r
import global_variables as g

directory_with_inputs = g.input_files
directory_with_files =  g.analysis_run_dir+"outputs/"
output_directory =      g.analysis_run_dir+"ntuples/grouped_samples/"

missing_dsids = set()


def get_files(dsid_list, directory) :
    """
    Grab all files in a given directory that have a desired DSID.
    Write all missing file DSIDs to output file missing_dsids 

    Args:
        dsid_list ( list(int) ): list of desired DSID  
        directory: directory containing files
    
    Returns:
        list(str): list of all file names that have a desired DSID
    """
    all_files = glob.glob(directory + "CENTRAL*.root")
    process_files = []
    for dsid in dsid_list :
        for root_file in all_files :
            if str(dsid) in root_file :
                process_files.append(root_file)
                break # move to next dsid 
        else: 
            missing_dsids.add(dsid)
    return process_files

def main() :

    for group, dsids in g.groups.iteritems():
        print "Running over " + group
        print "\tFound %d dataset ids"%len(dsids)
        
        root_files = get_files(dsids, directory_with_files)
        print "\tFound %d root files for the process"%len(root_files)

        if len(root_files) != len(dsids) :
            print "\tERROR Did not find expected number of files!"
            #sys.exit()
        chain = r.TChain("superNt","superNt")
        for root_f in root_files :
            chain.Add(root_f)
       
        # Create TChain file for each group
        
        file_name = "%s%s.root"%(output_directory,group)
        output_file = r.TFile(file_name,"RECREATE")
        chain.SetName("%s"%group)
        chain.Write()
        #chain.CloneTree(-1,"fast");
        #print "TChain written to", file_name 
        #output_file.Write() 
        print "\tProcess tree has %d total entries"%chain.GetEntries()
        output_file.Close()

        
    with  open(g.missing_dsids_file,'w') as ofile:
        for dsid in missing_dsids:
            ofile.write(str(dsid)+"\n")

if __name__ == '__main__':
    main()
