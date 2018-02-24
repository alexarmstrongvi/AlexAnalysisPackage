#!/bin/env python

#################################################
# This script will make a directory structure
# wherein each directory corresponds to a
# specific background and in it there will be
# the condor *.txt files made by the script
# "make_condor_lists.py" corresponding to that
# background process
#
# REQUIREMENTS BEFORE RUNNING :
#   You must have a directory containing
#   condor filelists of the same format
#   as those provided by the script
#   "make_condor_lists.py":
#    mc15_13TeV.<DSID>.foo.bar.txt
#
# daniel.joseph.antrim@cern.ch
# March 1 2016
##################################################

import os
import sys
import argparse
import subprocess
import glob
import global_variables as g

def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help = "Provide the directory " \
                        "containing the condor lists you wish to regroup",
                        required = True)
    args = parser.parse_args()
    indir = args.input
    outdir = g.input_files

    input_lists = []
    if os.path.isdir(indir) :
        input_lists = get_lists(indir)

    # find/build the groups in the user-provided input directory
    out_samples = {}
    for key, dsids in g.groups.iteritems() :
        out_samples[key] = []
        n_expected = len(dsids)
        for ds in dsids :
            for sample in input_lists :
                if str(ds) in sample :
                    out_samples[key].append(sample)
                    break
            else:
                print "WARNING Did not find DSID %d from sample group "\
                      "%s"%(ds, key)

    print "++ " + 75*"-" + " ++"
    # now that we have the groupings, make the directories 
    # to store the lists in
    if not outdir.endswith("/") :
        outdir = outdir + "/"
    for key in out_samples.keys() :
        dir_for_this_sample = outdir + key
        cmd = "mkdir -p %s"%dir_for_this_sample
        subprocess.call(cmd, shell=True)
        print "Storing sample %s with %d datasets in %s"%(key,
                                                        len(out_samples[key]),
                                                        dir_for_this_sample)
        for sample in out_samples[key] :
            cmd = "cp %s %s"%(sample, dir_for_this_sample)
        subprocess.call(cmd, shell=True)
    print "++ " + 75*"-" + " ++"

#########################################
# open up a directory and get all *.txt
# files
def get_lists(indir = "") :
    if indir == "" :
        print "get_lists ERROR: input directory is \"\""
        sys.exit()
    check_dir = indir
    if not indir.endswith("/") :
        check_dir += "/"
    lists_ = glob.glob(check_dir + "*.txt")
    print "get_lists Found %d lists in directory %s"%(len(lists_),
                                                             check_dir)
    return lists_

#########################################
if __name__=="__main__" :
    main()

