#!/usr/bin/env python
"""
================================================================================
Submit condor jobs for making flat ntuples
Examples
    ./submit_ntuple_maker.py -e makeFlatNtuples --zll_cr

Author:
    Alex Armstrong <alarmstr@cern.ch>
License:
    Copyright: (C) <June 23rd, 2018>; University of California, Irvine
================================================================================
"""



import sys, os, traceback, argparse
import time
import glob
import subprocess
import re
import global_variables as g

################################################################################
# Configuration options
doBrick = True
doLocal = False # Usually the same as the brick
doSDSC = False # We do not have the necessary permissions, jobs will hang
doUC = True 

FAX_FILES = True
################################################################################
# Format sensative functions
# Naming should match with what the executable expects
def get_region_name(args):
    # Names should match the options for the executable
    if args.do_baseline :
        return "baseline_sel"
    elif args.do_baseline_den :
            return "baseline_den"
    elif args.do_zjets_num_fakes :
            return "fake_num"
    elif args.do_zjets_den_fakes :
            return "fake_den"
    elif args.do_zll_cr :
            return "zll_cr"

def add_executable_ops(path_to_dataset, args):
        ops = "-i %s " % path_to_dataset
        ops += "--%s " % get_region_name(args)
        if args.apply_ff: ops += "--apply_ff "
        if args.suffix:
            ops += '-s %s ' % args.suffix
        return ops
################################################################################
def main ():
    """ Main Function """

    global args

    # Make needed condor scripts
    make_condor_script(brick = doBrick, local = doLocal, sdsc = doSDSC, uc = doUC)
    make_condor_executable()

    ############################################################################
    # Get dataset file names for submission

    # Get DSID set with samples to submit
    if args.dsid_list:
        with open(args.dsid_list) as f:
            dsids = {d.strip() for d in f if d.strip().isdigit()}
    else:
        dsids = {str(d) for d in g.get_all_dsids()}

    if args.apply_ff:
        data_dsids = set()
        import pdb
        for group, dsid_ls in g.get_data_groups().items():
            tmp_dsids = {str(s) for s in dsid_ls}
            data_dsids = data_dsids | set(tmp_dsids).intersection(dsids)
        dsids = data_dsids

    # Get list of sample files
    samples_dir = os.path.normpath(args.file_dir)
    all_sample_files = set(glob.glob("%s/*.txt"%samples_dir))

    # Select samples from list requested in DSID list
    # Check if duplicate DSIDs in the sample list
    sample_files = set()
    for d in dsids:
        flag = False
        for s in all_sample_files:
            if d in s and s in sample_files:
                print "WARNING :: DSID %s already found in %s" % (d, samples_dir)
            if d in s:
                sample_files.add(s)
                flag = True
        if not flag:
            print "WARNING :: DSID %s not found in %s" % (d, samples_dir)


    ############################################################################
    # Submit jobs for each sample
    for ii, dataset in enumerate(sample_files, 1):
        if args.apply_ff and 'data' not in dataset.split("/")[-1]: continue
        # TODO
        #
        # Get sample DSID
        # If sample in split list
        #   Open file and grab root file name
        #   Define relpath_to_dataset 

        d_name = dataset.split("/")[-1].replace(".txt", "")
        r_name = get_region_name(args)
        output_name = r_name + "_" + d_name
        if args.apply_ff:
            output_name = 'fakes_' + output_name
        if args.suffix: output_name += "_" + args.suffix
        log_file_path = os.path.normpath(os.path.join(args.log_dir, output_name))

        # Convert paths to work on remote site
        relpath_to_dataset = dataset.strip('/')

        #

        ########################################################################
        # Create run command

        # Positional arguments for condor executable
        # Order is important. See "make_condor_executable" for expected order
        # They get imported as an environment variable
        cmd_args = '%s '%args.out_dir
        cmd_args += '%s '%args.log_dir
        cmd_args += '%s '%args.executable
        cmd_args += '%s '%args.area_dir
        cmd_args += "%s " % add_executable_ops(relpath_to_dataset, args)
        out_file = log_file_path + ".out"
        log_file = log_file_path + ".log"
        err_file = log_file_path + ".err"
        append_queue_submission_to_condor_submit(cmd_args, err_file, out_file, log_file)

        # Execute run command
        print "[%d/%d] Adding to queue: %s" % (ii, len(sample_files), d_name)
        if args.dry_run:
            print
            print cmd_args
            print
    
    if not args.dry_run:
        subprocess.call('condor_submit submitFile.condor', shell=True)
        subprocess.call("condor_q $USER", shell=True)

################################################################################
# FUNCTIONS
def check_inputs(args):
    """ Check the input arguments are as expected """

    dirs_to_check = [args.out_dir, args.log_dir, args.file_dir]
    for d in dirs_to_check:
        if os.path.isdir(d): continue
        print "ERROR :: Could not find directory:", d
        return False

    if not os.path.isfile(args.area_tar):
        print "ERROR :: Cannot find tar file:", args.area_tar
        return False

    if args.dsid_list and not os.path.isfile(args.dsid_list):
        print "ERROR :: Cannot find dsid list:", args.dsid_list
        return False

    fake_region = args.do_baseline_den or args.do_zjets_den_fakes
    if args.apply_ff and not fake_region:
        return False

    all_regions = []
    all_regions.append(args.do_baseline)
    all_regions.append(args.do_baseline_den)
    all_regions.append(args.do_zjets_num_fakes)
    all_regions.append(args.do_zjets_den_fakes)
    all_regions.append(args.do_zll_cr)

    if all_regions.count(True) == 0:
        print "ERROR :: No region selected. Must pick one"
        return False
    elif all_regions.count(True) > 1:
        print "ERROR :: Multiple regions selected. Can only pick one"
        return False

    den_region = args.do_baseline_den or args.do_zjets_den_fakes
    if args.apply_ff and not den_region:
        print "ERROR :: Cannot apply fake factor on non-denominator region"
        return False

    pwd = os.environ['PWD']
    out_dir = os.path.abspath(args.out_dir)
    if pwd != out_dir:
        print "You must call this script from the output directory ",
        print "where the ntuples will be stored!"
        print " >>> Expected submission directory :", out_dir
        print " >>> You are running from PWD =", pwd
        return False

    return True

def check_environment():
    """ Check if the shell environment is setup as expected """
    assert os.environ['USER'], "USER variable not set"

    python_ver = sys.version_info[0] + 0.1*sys.version_info[1]
    assert python_ver >= 2.7, ("Running old version of python\n", sys.version)

def make_condor_script(brick = False, local = False, sdsc = False, uc = False) :
    f = open('submitFile.condor', 'w')
    f.write('universe = vanilla\n')
    f.write('+local=%s\n'%brick)
    f.write('+site_local=%s\n'%local)
    f.write('+sdsc=%s\n'%sdsc)
    f.write('+uc=%s\n'%uc)
    f.write('transfer_input_files = %s\n' % args.area_tar)
    f.write('should_transfer_files = YES\n')
    f.write('when_to_transfer_output = ON_EXIT\n')
    f.write('use_x509userproxy = True\n')
    f.write('notification = Never\n')
    f.close()

def append_queue_submission_to_condor_submit(arguments, err_file, out_file, log_file):
    f = open("submitFile.condor", 'a')
    f.write('\n')
    f.write('executable = RunCondorSF.sh\n')
    f.write('arguments = %s\n' % arguments)
    f.write('error = %s\n' % err_file)
    f.write('output = %s\n' % out_file)
    f.write('log = %s\n' % log_file)
    f.write('queue\n')
    f.close()

def make_condor_executable() :
    f = open('RunCondorSF.sh', 'w')
    f.write('#!/bin/bash\n\n\n')
    f.write('echo " ------- RunCondorSF -------- "\n')
    f.write('output_dir=${1}\n')
    f.write('log_dir=${2}\n')
    f.write('sflow_exec=${3}\n')
    f.write('stored_dir=${4}\n')
    f.write('sflow_options=${@:5}\n\n')
    f.write('echo "    output directory   : ${output_dir}"\n')
    f.write('echo "    log directory      : ${log_dir}"\n')
    f.write('echo "    sflow executable   : ${sflow_exec}"\n')
    f.write('echo "    tarred dir         : ${stored_dir}"\n')
    f.write('echo "    sflow options      : ${sflow_options}"\n\n')
    f.write('while (( "$#" )); do\n')
    f.write('    shift\n')
    f.write('done\n\n')
    f.write('work_dir=${PWD}\n')
    f.write('echo "current directory structure:"\n')
    f.write('ls -ltrh\n\n')
    f.write('echo "untarring area.tgz"\n')
    f.write('tar -xzf area.tgz\n\n')
    f.write('echo "done untarring"\n')
    f.write('export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n')
    f.write('source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh\n')
    f.write('echo "1 current directory structure:"\n')
    f.write('ls -ltrh\n\n')
    f.write('echo "Calling : cd ${stored_dir}"\n')
    f.write('cd ${stored_dir}\n')
    f.write('echo "Directory structure:"\n')
    f.write('ls -ltrh\n')
    f.write('lsetup fax\n')
    f.write('source susynt-read/bash/setup_root.sh\n')
    f.write('echo "Calling : source RootCoreBin/local_setup.sh"\n')
    f.write('source RootCoreBin/local_setup.sh\n')
    f.write('echo "Calling : cd AlexAnalysisPackage/bash"\n')
    f.write('cd AlexAnalysisPackage/bash\n')
    f.write('source setRestFrames.sh\n')
    f.write('echo "Calling : cd ${work_dir}"\n')
    f.write('cd ${work_dir}\n')
    f.write('echo "Calling : ${sflow_exec} ${sflow_options}"\n')
    f.write('${sflow_exec} ${sflow_options}\n')
    f.write('echo "final directory structure:"\n')
    f.write('ls -ltrh\n')
    #f.write('echo "Printing environment"\n\n')
    #f.write('( set -o posix ; set ) | sort')
    f.close()

def is_data_sample(did):
    return re.search("data1[5-8]", dataset.split("/")[-1])

def get_args():
    _executable = "makeFlatNtuples"
    _area_dir   = g.analysis_dir[1:] if g.analysis_dir.startswith('/') else g.analysis_dir
    _area_tar   = os.path.normpath(os.path.join(g.analysis_path,'area.tgz'))
    _out_dir    = g.output_dir
    _log_dir    = g.logs_dir
    if FAX_FILES:
        _file_dir   = g.fax_input_files
    else:
        _file_dir   = g.local_input_files

    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-e', '--executable',
                        default=_executable,
                        help='name of executable')
    parser.add_argument('-a', '--area_dir',
                        default=_area_dir,
                        help='path to executable area in tarred file')
    parser.add_argument('-t', '--area_tar',
                        default=_area_tar,
                        help='file of tarred area for running executable')
    parser.add_argument('-o', '--out_dir',
                        default=_out_dir,
                        help='directory for returned executable outputs')
    parser.add_argument('-l', '--log_dir',
                        default=_log_dir,
                        help='directory for returned executable logs (.out, .log, .err)')
    parser.add_argument('-f', '--file_dir',
                        default=_file_dir,
                        help= 'directory of sample files with sample file links')
    parser.add_argument('--dsid_list',
                        help='list of DSID files to look for in sample directory')
    parser.add_argument('-s','--suffix',
                        default='',
                        help='suffix to append to output and log names')
    parser.add_argument('--do_baseline', action='store_true')
    parser.add_argument('--do_baseline_den', action='store_true')
    parser.add_argument('--do_zjets_num_fakes', action='store_true')
    parser.add_argument('--do_zjets_den_fakes', action='store_true')
    parser.add_argument('--do_zll_cr', action='store_true')
    parser.add_argument('--apply_ff', action='store_true')
    parser.add_argument('--dry_run', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    args.area_tar = os.path.abspath(args.area_tar)

    if not check_inputs(args):
        sys.exit()

    return args


################################################################################
# Run main when not imported
if __name__ == '__main__':
    try:
        start_time = time.time()
        args = get_args()
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

