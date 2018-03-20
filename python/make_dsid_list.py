#!/bin/bash/env python
"""
Program make_dsid_list.py
Author: Alex Armstrong <alarmstr@cern.ch>
Copyright: (C) Jan 6th, 2018; University of California, Irvine
"""
import sys
import os
from argparse import ArgumentParser
from collections import defaultdict
import global_variables as g
import python_tools as tools

def get_desired_samples_map(args):
    """ Get map of DSID lists from global list"""
    if args.all:
        return g.groups
    elif args.data:
        return g.get_data_groups()
    elif args.bkgd:
        return g.get_bkgd_groups()
    elif args.mc:
        return g.get_mc_groups()
    return g.groups

def get_local_samples_map():
    """ Get map of samples available locally.
        Map keys are the sample groups.
        Directory storing local files defined in global module.
    """
    samples_map = defaultdict(list)
    for key, group in g.LOCAL_DSID_SUBDIRS.items():
        #reject_dir = ((args.data and group != g.LOCAL_DSID_DIR['data'])
        #        or ((args.mc or args.bkgd) and group != g.LOCAL_DSID_DIR['mc']))
        #if reject_dir:
        #    continue
        path = '%s/%s'%(g.LOCAL_DSID_DIR,group)
        for sample in os.listdir(path):
            dsid = tools.get_dsid_from_sample(sample)
            group_name = tools.get_sample_group(sample)
            if dsid:
                samples_map[group_name].append(sample)
    return samples_map

def get_fax_samples_map(store_empty=False):
    """ Get map of samples available locally.
        Map keys are the sample groups.
        Directory storing local files defined in global module.
    """
    samples_map = defaultdict(list)
    path1 = g.input_files
    for groups in tools.get_list_of_subdirectories(path1):
        path2 = '%s%s'%(path1,groups)
        for sample in tools.get_list_of_files(path2):
            path_to_file = '%s%s'%(path2,sample)
            file_size = os.path.getsize(path_to_file)
            if store_empty != (file_size <= 0):
                continue
            if sample.endswith('.txt'):
                sample = sample[:-4]
            group_name = tools.get_sample_group(sample)
            samples_map[group_name].append(sample)
    return samples_map

def get_todo_samples_map(args,to_dwnld=False,to_make=False):
    # Get map from dsid to sample name for main categories
    category_map = {'desired' : get_desired_samples_map(args),
                   'local' : get_local_samples_map(),
                   'fax' : get_fax_samples_map()}
    dsid_sample_map = {}
    for key, samples_map in category_map.iteritems():
        dsid_sample_map [key] = {}
        for group, samples in samples_map.iteritems():
            for sample in samples:
                if key == 'desired':
                    dsid = str(sample)
                else:
                    dsid = tools.get_dsid_from_sample(sample)
                if not dsid:
                    continue
                dsid_sample_map[key][dsid] = sample

    # Get relevent sets of DSIDs
    samples_map = defaultdict(list)
    fax_dsids = dsid_sample_map['fax']
    local_dsids = dsid_sample_map['local']
    for dsid in dsid_sample_map['desired']:
        for_dwnld = dsid in fax_dsids and dsid not in local_dsids
        for_making = dsid not in fax_dsids and dsid not in local_dsids
        if (to_dwnld and not for_dwnld) or (to_make and not for_making):
            continue

        # Get full sample name from dsid if possible
        if to_dwnld:
            sample_name = dsid_sample_map['fax'][dsid]
            sample_name += '_nt'
            if 'SusyNt' in sample_name:
                sample_name = 'group.phys-susy:' + sample_name
            else:
                print 'Unknown scope name for', sample_name
                sys.exit()
            group = tools.get_sample_group(sample_name)
        elif to_make:
            sample_name = dsid
            group = 'Find full DAOD sample name'
        samples_map[group].append(sample_name)
    return samples_map

def main():
    """ Main function """
    print "\nMaking DSID list..."

    parser = ArgumentParser()
    parser.add_argument('-o', '--output',
                        help='output path and file name')
    parser.add_argument('--local',
                        action='store_true',
                        help='output only locally available samples')
    parser.add_argument('--empty',
                        action='store_true',
                        help='output only empty rucio samples')
    parser.add_argument('--fax',
                        action='store_true',
                        help='output only fax available sample DSIDs')
    parser.add_argument('--dwnld',
                        action='store_true',
                        help='output only sample DSIDs to be downloaded')
    parser.add_argument('--daod',
                        action='store_true',
                        help='output only sample DSIDs to be produced')
    parser.add_argument('--all',
                        action='store_true',
                        help='output all DSIDs in global_variables')
    parser.add_argument('--data',
                        action='store_true',
                        help='output only data DSIDs')
    parser.add_argument('--bkgd',
                        action='store_true',
                        help='output only background MC DSIDs')
    parser.add_argument('--mc',
                        action='store_true',
                        help='output only mc DSIDs')
    args = parser.parse_args()
    output = args.output

    # Checks
    if sum([args.data, args.bkgd, args.mc, args.dwnld]) > 1:
        print 'ERROR :: More than one DSID type specified.'\
              'Specify at most one of data, bkgd, mc, and dwnld.'
        sys.exit()

    # Determine output name
    if not output:
        if args.local:
            prefix = 'local'
        elif args.dwnld:
            prefix = 'download_these'
        elif args.daod:
            prefix = 'daod_production'
        elif args.fax:
            prefix = 'fax'
        elif args.empty:
            prefix = 'empty'
        elif args.all:
            prefix = 'all'
        elif args.data:
            prefix = 'data'
        elif args.bkgd:
            prefix = 'bkgd'
        elif args.mc:
            prefix = 'mc'
        else:
            prefix = 'analysis'
        output = '%s%s_dsids.txt'%(g.dsid_dir, prefix)
    with open(output,'w') as ofile:
        if args.local:
            groups = get_local_samples_map()
        elif args.dwnld:
            groups = get_todo_samples_map(args,to_dwnld=True)
        elif args.daod:
            groups = get_todo_samples_map(args,to_make=True)
        elif args.fax:
            groups = get_fax_samples_map()
        elif args.empty:
            groups = get_fax_samples_map(store_empty=True)
        else:
            groups = get_desired_samples_map(args)

        # Write to output
        for group in sorted(groups, key=lambda k:len(groups[k])):
            lst = groups[group]
            ofile.write('===== %s (%d) ===== \n'%(group, len(lst)))
            for dsid in lst:
                ofile.write('%s\n'%dsid)
            ofile.write('\n')
    print 'Output saved at %s\n'%output

if __name__ == '__main__':
    # Do not execute main() when script is imported as a module
    main()
