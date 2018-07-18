#!/bash/env python

import sys, os, traceback, argparse
import time
from glob import glob
import ROOT as r

################################################################################
def main ():
    """ Main Function """

    global args

    # Get all files in directory

    # Loop over each file in the directoty
    nevents = {}
    for f in  glob("./*root"):
        print f
        rf = r.TFile(f, 'r')
        t = rf.Get('superNt')
        nevents[f] = t.GetEntries()
        rf.Close()

    # Print map into output file
    with open(args.output, 'w') as ofile:
        for name, yld in nevents.iteritems():
            ofile.write("%s : %d\n" % (name, yld))

    print "DONE"


################################################################################
# FUNCTIONS
def get_args():
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--output', help='output')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='verbose')
    args = parser.parse_args()
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
