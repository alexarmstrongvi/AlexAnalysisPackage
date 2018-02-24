#!/usr/bin/env python

##############################
#
# Makes a new filelist directory 
# with samples in the current filelist
# directory that are available in on the brick.
# It also sets the correct prefix 
#
#
#
##############################

import subprocess
import glob
import sys
import os

current_filelist_dir = G.analysis_dir+"BACKUP_inputs_LFV/"
new_filelist_dir = G.analysis_dir+"inputs_LFV_new/"
new_prefix = "root://${ATLAS_XROOTD_CACHE}//data/uclhc/uci/user/dantrim/susynt_productions/n0233/"

samples = ["data_all", "MC"]#, "data"]

def get_sample_names() :

    out = []
    tmp = glob.glob(current_filelist_dir + "*")

    for t in tmp :
        t = t.split("/")[-1]
        out.append(t)

    return out

def get_dsids(process_name) :
    file_to_check = "%s%s/"%(current_filelist_dir, process_name)

    files = glob.glob(file_to_check + "*.txt")

    out = []
    for f in files :
        is_data = "data" and "Main" in f

        f = f.split("/")[-1].replace(".txt","")
        out.append(f)
        #if is_data :
        #    f = f.split(".00")[1].split(".")[0]
        #else :
        #    f = f.split("mc15_13TeV.")[1].split(".")[0]
        #out.append(f)

    return out

def get_new_files(dataset) :

    #if "284154" in dataset :
    #    return []

    
    check_name = "/data/uclhc/uci/user/dantrim/susynt_productions/n0233/"
    suffix = "mc/"
    if "data" and "Main" in dataset :
        suffix = "data/"
    check_name += suffix

    check_name += dataset
    if not check_name.endswith("n0233"): check_name = check_name[:-1]
    check_name += "c_nt/"
    samples = glob.glob(check_name + "*.root")
    if len(samples) == 0: 
        samples = glob.glob(check_name.replace("c_nt/","b_nt/") + "*.root")
    if len(samples) == 0: 
        samples = glob.glob(check_name.replace("c_nt/","_nt/") + "*.root")
    
    return samples

def make_new_dir(sample) :

    cmd = "mkdir -p %s%s/"%(new_filelist_dir, sample)
    subprocess.call(cmd, shell=True)


def main() :

    sample_names = get_sample_names()

    tmp_samples = []
    for s in sample_names :
        if s not in samples : continue
        tmp_samples.append(s)
    sample_names = tmp_samples

    #sample_dict = {} # { process name : list of dsids }
    #root://${ATLAS_XROOTD_CACHE}//data/uclhc/uci/user/dantrim/susynt_productions/n0232/mc/

    for s in sample_names :
        dataset_names = get_dsids(s)

        make_new_dir(s)

        counter = 0
        for ds in dataset_names :
            list_of_new_files = get_new_files(ds)
            if len(list_of_new_files)==0 : 
                print "No match for " + ds
                counter += 1
                continue
                

            oname = new_filelist_dir + s + "/" + ds + ".txt" 
            ofile = open(oname, 'w')

            prefix = "root://${ATLAS_XROOTD_CACHE}/"
            for f in list_of_new_files :
                line = prefix + f + '\n'
                ofile.write(line)
            ofile.close

        print "Number of skipped files = ", counter        

        #print 55*"-"
        #print "%s"%s
        #for x in dsids :
        #    print "   > %s"%x



#___________________________________________________
if __name__ == "__main__":
    main()
