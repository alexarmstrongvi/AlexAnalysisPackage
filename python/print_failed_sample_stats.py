#!/bin/bash/env python

import python_tools as Tools
import os
import sys
import global_variables as g

directory_with_logs = g.logs_dir
missing_dsid_file = g.missing_dsids_file

def main():
    print "Checking failed samples..." 
    failure_messages = ["No servers are available to read the file",\
                        "no url specified",\
                        "Redirect limit has been reached",\
                        "Operation expired",\
                        "no such file or directory",\
                        "no message of desired type",\
                        "permission denied",\
                        "is disallowed",\
                        "Invalid address",\
                        "Object is in \'zombie\' state",\
                        "-3 in inflate (zlib)",\
                        "Register StreamerInfo for Susy::Muon on non-empty slot",\
                        "Socket error",\
                        #"which: no python in ((null))",\
                        "FATAL Cross-section is negative",\
                        "<EmptyFile>"]
    fail_count = {}
    for fm in failure_messages:
        fail_count[fm] = 0
    unmatched_dsids = []
    missed_dsids = []
    log_files = Tools.get_list_of_files(directory_with_logs)
    err_files = [x for x in log_files if x.endswith('.err')]
    out_files = [x for x in log_files if x.endswith('.out')]
    dsid_list = open(missing_dsid_file,'r').readlines()

    for dsid in dsid_list: 
        dsid = dsid.strip()
        found_dsid = False
        for err_file in err_files:
            if dsid not in err_file: continue
            if found_dsid: print 'Duplicate err file = ', dsid
            else: found_dsid = True
            with open(directory_with_logs+err_file,'r') as tmp_file:
                flag = False
                if len(tmp_file.readlines())==0:
                    fail_count["<EmptyFile>"] += 1
                    continue
                for fm in failure_messages:
                    tmp_file.seek(0)
                    if fm not in tmp_file.read(): continue
                    if flag: print "Duplicate with: " + fm
                    else: flag = True
                    fail_count[fm] += 1
                if not flag: unmatched_dsids.append(dsid)
        if not found_dsid: missed_dsids.append(dsid)


    tmp_list = list(unmatched_dsids)
    print "Missing Cross Section:"
    for dsid in tmp_list:
        for out_file in out_files:
            if dsid not in out_file: continue
            with open(directory_with_logs+out_file,'r') as tmp_file:
                if "FATAL Cross-section is negative" in tmp_file.read():
                    fail_count["FATAL Cross-section is negative"] += 1 
                    unmatched_dsids.remove(dsid)
                    print "\t",out_file
                else:
                    print "? " + out_file
    print "Results:"
    for fm in failure_messages:
        percentage = 100*fail_count[fm]/float(len(dsid_list))
        print "%3i/%3i (%3.1f%%)\twith %s"%(fail_count[fm],len(dsid_list),percentage,fm)
    percentage  = 100*len(unmatched_dsids)/float(len(dsid_list))
    percentage2 = 100*len(missed_dsids)/float(len(dsid_list))
    print "%3i/%3i (%3.1f%%)\t     Unmatched - unrecognized error message"%(
            len(unmatched_dsids),len(dsid_list),percentage)
    print "%3i/%3i (%3.1f%%)\t     Missed - No log file for DSID\n"%(
            len(missed_dsids),len(dsid_list),percentage2)
    print "Unmatched:", unmatched_dsids
    print "Missed:", missed_dsids

if __name__ == "__main__":
    main()

