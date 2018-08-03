#!/usr/bin/env python
"""
== make_local_condor_lists ==
Make text files for each local sample for use in condor submissions

Examples:
    python make_local_condor_lists.py
Author: 
    Alex Armstrong <alarmstr@cern.ch>
License:
    Copyright: (C) <Mar 7, 2018>; University of California, Irvine 
"""

import sys, os, traceback, argparse
import time
import global_variables as g
import pyTools as tools
import subprocess

################################################################################
def main ():
    """ Main Function """
    # Get list of sample names 
    local_files = g.get_local_files()
    condor_lists_made = False

    # Remove anything that is not a SusyNt container
    local_files = [x[:-3] for x in local_files if x.endswith("_nt")]

    missing_dsids = dict( (key, []) for key in g.groups)
    # Loop over all the desired DSIDs and grab the most up to date sample name
    # available in the input file list. Record all DSIDs that are not found
    for group, dsids in g.groups.iteritems():
        for dsid in dsids:
            matched_files = []
            for sample in local_files:
                if str(dsid) in sample:
                    matched_files.append(sample)
            if matched_files:
                if len(matched_files) == 1:
                    sample = matched_files[0]
                else:
                    # Use latest tag version (e.g. n0235c over n0235b)
                    # This method is format sensative
                    def get_tag_ver(name):
                        search_str = g.tag
                        check_spot = name.find(search_str)+len(search_str)
                        return name[check_spot]
                    sample = sorted(matched_files, key=get_tag_ver)[-1]
                # Create list file for condor
                group_dir = "%s/"%group 
                sample_dir = "%s_nt/"%sample
                if create_condor_list(group_dir, sample_dir):
                    condor_lists_made = True
            else:
                missing_dsids[group].append(dsid)
    if condor_lists_made:            
        print "Local file condor lists stored at", g.local_input_files
    else:
        print "All lists up-to-date"
    
    # Save missing information
    missing_name = 'missing_local_dsids.txt'
    with tools.cd(g.analysis_run_dir):
        print "Saving missing dsids at",g.analysis_run_dir+missing_name
        ofile = open(missing_name,'w')
        for group, dsids in missing_dsids.iteritems():
            ofile.write("===== %s (%d) =====\n"%(group,len(dsids)))
            for dsid in dsids:
                ofile.write(str(dsid)+'\n')
        ofile.close()



################################################################################
# FUNCTIONS
################################################################################

def create_condor_list(group_dir, path_to_sample):
    """
    Create a condor text file for the provided sample directory
    args:
        group_dir (str): group name for directory holding sample
        path_to_sample (str): full path name to sample directory
    return:
        (bool): indicates if any condor lists made or updated 

    """
    def get_file_number(file_name):
        for c in reversed(root_file):
            if c.isdigit():
                return int(c)

    condor_lists_made = False 
    root_files = sorted(os.listdir(path_to_sample))
    if not len(root_files):
        print "No root files found in", path_to_sample
        return condor_lists_made

    with tools.cd(g.local_input_files):
        sample_name = path_to_sample.split('/')[-2]
        if sample_name.endswith('_nt'):
            sample_name = sample_name[:-3]
        ofile_name = group_dir+sample_name+'.txt'
        if os.path.exists(ofile_name):
            nfiles_old = sum(1 for line in open(ofile_name))
            if nfiles_old < len(root_files):
                print "Updating", ofile_name
                shell_cmd = 'rm %s'%ofile_name
                subprocess.call(shell_cmd, shell=True)
                condor_lists_made = True
            elif nfiles_old > len(root_files):
                print "Current file has more links. Not updated", ofile_name
                return condor_lists_made
            else:
                return condor_lists_made

        if not os.path.exists(group_dir):
            os.makedirs(group_dir)
        ofile = open(ofile_name,'w')
        condor_lists_made = True
        for root_file in root_files:
            file_num = get_file_number(root_file)

            condor_link = "%s%s%s"%(g.local_prefix, path_to_sample, root_file)
            ofile.write("%s\n"%condor_link)
        ofile.close
        return condor_lists_made

    

################################################################################
if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = argparse.ArgumentParser(
                description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('--fax', 
                            action='store_true', default=False, 
                            help='make fax list for inputs')
        parser.add_argument('-v', '--verbose', 
                            action='store_true', default=False, 
                            help='verbose output')
        args = parser.parse_args()
        if args.verbose: 
            print '>'*40
            print 'Running {}...'.format(os.path.basename(__file__))
            print time.asctime()
        main()
        if args.verbose: 
            print time.asctime()
            time = (time.time() - start_time)
            print 'TOTAL TIME: %fs'%time,
            print ''
            print '<'*40
        #sys.exit(0)
    except KeyboardInterrupt, e: # Ctrl-C
        print 'Program ended by keyboard interruption'
        raise e
    except SystemExit, e: # sys.exit()
        print 'Program ended by system exit'
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)
