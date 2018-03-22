#!/usr/bin/env python
"""
================================================================================
Make a yield table
Examples
    python yield_table_maker.py -c config_file.conf

Author:
    Alex Armstrong <alarmstr@cern.ch>
TODO: Licence:
    Copyright: (C) <Mar 20th; University of California, Irvine
================================================================================
"""

import sys, os, traceback, argparse
import time
import ROOT

################################################################################
def main ():
    """ Main Function """
    global args
    print_user_arguments()

    # get the config file
    conf_file = get_plotConfig(g_plotConfig)
    print "Found the configuration file: %s"%conf_file

    # variables to be filled
    g_data = None
    g_backgrounds = []
    g_cuts = []
    execfile(conf_file)

    make_yield_table(g_data, g_backgrounds, g_cuts)


################################################################################
# FUNCTIONS
def print_user_arguments():
    """ print arguments given as inputs """

    print " ++ ------------------------- ++ "
    print "      User Inputs                "
    print "                                 "
    print " config file      :  %s          "%args.conf_file_path
    print " output directory :  %s          "%args.output_path
    print " verbose          :  %s          "%args.verbose
    print "                                 "
    print " ++ ------------------------- ++ \n"

    if not os.path.exists(args.conf_file_path):
        print "ERROR :: Unable to find file:", args.conf_file_path
        sys.exit()

    if os.path.exists(args.output_path) and not args.overwrite:
        print "ERROR :: Output file already exists:", args.output_path
        print "         Either delete the file, pick a different name,"
        print "         or use the --overwrite option."
        sys.exit()

def make_yield_table(data, backgrounds, cuts):


    if args.steps:
        cumulative_cut = ROOT.TCut("")
        for cut in cuts:
            for bkgd in backgrounds:
                # The TTree in each background keeps the same event list
                # from before so only needs to be updated with new cut
                list_name = "tmp_list_" b.treename
                draw_list = ">> " + list_name
                b.tree.Draw(draw_list, cut)
                elist = r.gROOT.FindObject(list_name)
                b.tree.SetEventList(elist)
                ylow = args.variable_range[0]
                yhigh = args.variable_range[1]
                hist_name = 'hist'
                hist = TH1D(hist_name,hist_name,1,ylow,yhigh)


            if data :
                data_list_name = "tmp_list_" + data.treename
                draw_list = ">> " + data_list_name
                data.tree.Draw(draw_list, cut)
                data_list = r.gROOT.FindObject(data_list_name)
                data.tree.SetEventList(data_list)
################################################################################
# Run main when not imported
if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = argparse.ArgumentParser(
                description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('conf_file_path',
               help='configuration file')
        parser.add_argument('-v', '--variable'
               default='isMC',
               help='Plotting variable, assumed boolean [default: isMC]')
        parser.add_argument('-r', '--variable_range'
               default=[-1,2], type=list
               help='y range for variable [default: [-1,2] ]')
        parser.add_argument('-o', '--output_path',
               default='yield_table.txt',
               help='path name for output file')
        parser.add_argument('--steps',
               action='store_true', default=False,
               help='print yield table for at each cut')
        parser.add_argument('--latexify',
               action='store_true', default=False,
               help='print yield table in latex format')
        parser.add_argument('--overwrite',
               action='store_true', default=False,
               help='overwite output file if it exists')
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
