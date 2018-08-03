#!/bin/bash/env python

import pyTools as Tools
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
    # Add failed DSIDs to output file if not already there (don't add empty files)
    failed_file_count = 0
    stored_dsid_count = 0
    empty_file_count = 0
    print "Checking %d files for 'Files OK'" % len(out_files)
    for out_file in out_files:
        output_info = open(directory_with_logs+out_file).read()
        if 'Files OK' in output_info:
            continue
        if not output_info: 
            empty_file_count+=1
            continue
        DSID = Tools.strip_string_to_substring(out_file,'[1-9][0-9][0-9][0-9][0-9][0-9]')
        if DSID: 
            failed_file_count+=1
            if DSID not in missing_dsids: 
                ofile.write(DSID+'\n')
                stored_dsid_count+=1
    
    # Wrap up
    ofile.close()  
    if empty_file_count:
        print "\t%i files are still running or empty" % empty_file_count
    if failed_file_count:
        print "\t%i/%i failed samples stored"%(stored_dsid_count,failed_file_count)
    elif not empty_file_count:
        print "\tAll files are good :)"
    if stored_dsid_count > 0:
        if append_write == 'a': print "\tDSIDs added to " + missing_dsid_file.split('/')[-1]
        if append_write == 'w': print "\tDSIDs stored at " + missing_dsid_file
if __name__ == "__main__":
    main()
