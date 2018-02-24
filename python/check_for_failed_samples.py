#!/bin/bash/env python

import python_tools as Tools
import os
import global_variables as g

directory_with_logs = g.logs_dir
missing_dsid_file = g.missing_dsids_file

def main():
    """
    Finds output log files without a success message (e.g. "Files OK") 
    and adds the DSIDs to the missing DSID text file
    """
    print "Checking for failed samples..." 

    # Open missing dsid text file
    # If exists, store DSIDs and append additional DSIDs 
    missing_dsids = ""
    if os.path.exists(missing_dsid_file):
        append_write = 'a'
        ofile = open(missing_dsid_file,'r')
        missing_dsids = ofile.read()
        ofile.close()
    else:
        append_write = 'w'
    ofile = open(missing_dsid_file,append_write)

    # Get output log files
    log_files = Tools.get_list_of_files(directory_with_logs)
    out_files = [x for x in log_files if x.endswith('.out')]
    
    # Check which output files have success message
    # Add missing/failed DSIDs to output file if not already there
    failed_file_count = 0
    stored_dsid_count = 0
    for out_file in out_files:
        if 'Files OK' in open(directory_with_logs+out_file).read(): continue
        DSID = Tools.strip_string_to_substring(out_file,'[1-9][0-9][0-9][0-9][0-9][0-9]')
        if DSID: 
            failed_file_count+=1
            if DSID not in missing_dsids: 
                ofile.write(DSID+'\n')
                stored_dsid_count+=1
    
    # Wrap up
    ofile.close()
    print "\t%i/%i failed samples stored"%(stored_dsid_count,failed_file_count)
    if append_write == 'a': print "\tDSIDs added to " + missing_dsid_file.split('/')[-1]
    if append_write == 'w': print "\tDSIDs stored at " + missing_dsid_file
if __name__ == "__main__":
    main()
