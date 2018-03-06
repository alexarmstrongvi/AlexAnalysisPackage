#1/bin/bash/env python

import os
import glob
import re

def main():
    print "Testing Tools"
    directory_with_subdir = 'LOCAL_inputs_LFV'
    directory_with_files = '/data/uclhc/uci/user/armstro1/SusyNt/analysis_n0232_run/scripts'
    list_of_subdirectories = get_list_of_subdirectories(directory_with_subdir)
    list_of_files = get_list_of_files(directory_with_files)
    print list_of_subdirectories
    print list_of_files
    list_of_files = strip_strings_to_substrings(list_of_files,'[aeiou]')
    print list_of_files

# output list of directories in directory
def get_list_of_subdirectories(directory,search_string='*/'):
    if not directory.endswith('/'): directory+='/'
    if not search_string.endswith('/'): search_string += '/'
    subdirectories = [(x.split('/')[-2]+'/') for x in glob.glob(directory+search_string)]
    if len(subdirectories)==0: print "WARNING: %s has no subdirectories"%directory
    return subdirectories

# output list of files in directory
def get_list_of_files(directory, search_string='*'):
    if not directory.endswith('/'): directory+='/'
    files = [x.split('/')[-1] for x in glob.glob(directory+search_string)]
    if len(files)==0: print "WARNING: %s has no files"%directory
    return files

# strip input string list into substring list
def strip_strings_to_substrings(string_list,substring):
    substring_list = []
    for s in string_list:
        substring_list.append(strip_string_to_substring(s,substring))
    if all(x == '' for x in substring_list): print "WARNING: %s not found in any entry"%substring
    return substring_list

def strip_string_to_substring(string,substring):
    match = re.search(r'%s'%(substring),string)
    if match:
        return match.group()
    else:
        return ''
def pprint_dict(input_dict, tabs='', max_display=20):
    if type(input_dict) is None:
        return
    elif type(input_dict) is not dict:
        print('{}{}'.format(tabs,input_dict))
        return

    for key,val in input_dict.iteritems():
        val_t = type(val)
        print '{}{})'.format(tabs,key)
        tabs += '\t'
        if val_t is dict:
            pprint_dict(val,tabs)
        elif val_t is list:
            counter = 0
            prior_value = None
            for x in val:
                # Print dicts no matter what
                if type(x) is dict:
                    pprint_dict(prior_value,tabs)
                    pprint_dict(x,tabs)
                    counter = max_display-1
                elif type(x) is list:
                    pprint_dict(prior_value,tabs)
                    print 'Not able to handle lists within lists'
                    print '{}{})'.format(tabs,x[:max_display])
                    counter = max_display-1
                    return
                # Skip non-dicts if max reached
                if counter <= max_display:
                    pprint_dict(x,tabs)

                if counter == max_display:
                    print '{}. . .'.format(tabs)
                counter += 1
                prior_value = x
            else:
                pprint_dict(x,tabs)
                
        else: 
            pprint_dict(val,tabs)
        tabs = tabs[1:] # remove one tab
    return

def to_ordinal(value):
    """
    Converts zero or a *postive* integer (or their string 
    representations) to an ordinal value.

    """
    try:
        value = int(value)
    except ValueError:
        print "ERROR (to_ordinal) :: Input value must be an integer"
        return value

    if value % 100//10 != 1:
        if value % 10 == 1:
            ordval = u"%d%s" % (value, "st")
        elif value % 10 == 2:
            ordval = u"%d%s" % (value, "nd")
        elif value % 10 == 3:
            ordval = u"%d%s" % (value, "rd")
        else:
            ordval = u"%d%s" % (value, "th")
    else:
        ordval = u"%d%s" % (value, "th")

    return ordval

def get_unique_elements(list1, list2):
    unique_list1 = []
    for x in list1:
        if x not in list2:
            unique_list1.append(x)
    unique_list2 = []
    for x in list2:
        if x not in list1:
            unique_list2.append(x)
    return unique_list1, unique_list2

def get_shared_elements(list1, list2):
    shared_list = []
    for x in list1:
        if x in list2:
            shared_list.append(x)
    return shared_list

def get_duplicates(lst):
    """ Find and return a list of duplicates for an input list"""
    seen = set()
    dups = []
    for x in lst:
        if x in seen:
            dups.append(x)
        else:
            seen.add(x)
    return dups

def all_elements_equal(input_list):
    """ Checks if all elements in an iterable are the same"""
    return len(set(input_list)) <= 1

def good_matrix_shape(matrix):
    """ Check if matrix has the same number of elements in each row"""
    column_len = [len(x) for x in matrix]
    return all_elements_equal(column_len)

def get_matrix_shape(matrix):
    """ Return matrix dimensions """
    if not good_matrix_shape(matrix):
        print "WARNING (get_matrix_shape) :: "\
              "Input list of lists is not a matrix. "\
              "Returning [0,0]"
        return 0, 0
    row_size = len(matrix)
    column_size = len(matrix[0])
    return row_size, column_size

# Analysis specific
def get_dsid_from_sample(sample):
    """ Extract DSID from sample name"""
    sample = sample.strip()
    search_str = '[0-9]{6}(?=\.)'
    dsid = strip_string_to_substring(sample,search_str)
    if not dsid:
        dsid = strip_string_to_substring(sample,'[0-9]{6}')
    return dsid

def get_dsid_sample_map(sample_list):
    """
    Get map of DSIDs to full sample name
    given a list of full sample names
    """
    dsid_sample_map = {}
    for sample in sample_list:
        dsid = get_dsid_from_sample(sample)
        if not dsid:
            continue
        if dsid in dsid_sample_map:
            print 'WARNING :: DSID %s is duplicated'%dsid
        dsid_sample_map[dsid] = sample
    return dsid_sample_map

def trim_sample_name(sample):
    dsid = get_dsid_from_sample(sample)
    if not dsid:
        print 'WARNING (trim_sample_name) :: '\
              'Input sample name has no DSID'
        return ''
    pos = sample.find(dsid)
    sample = sample[pos:]
    # sample names are just DSIDs
    if len(sample) == 6:
        return sample
    # unexpected format with spaces
    if sample.find(' ') != -1:
        return sample
    
    sample = sample.split('.')[:-1]
    sample = '.'.join(sample)
    return sample

def get_sample_group(sample):
    searches = {}
    # Match pattern order
    # Elements in the innermost lists are OR'd
    # The list elements in main list are AND'd
    # Data
    searches['data15'] = [['data15']]
    searches['data16'] = [['data16']]
    # Higgs -> tautau
    searches['Htt'] = [ ['H125'],['tautau','tt'] ]
    # Higgs -> WW
    searches['HWW'] = [ ['H125'],['WW'] ]
    # Z -> tautau + jets
    # Z -> tautau + jets
    searches['Ztt'] = [['Ztautau','Ztt']]
    # Z -> ee + jets
    searches['Zee'] = [['Zee']]
    # Z -> mumu + jets
    searches['Zmumu'] = [['Zmumu','Zmm']]
    # Z -> vv + jets
    searches['Znunu'] = [['Znunu']]
    # W + jets
    searches['Wjets'] = [['Wplus','Wminus','W(e|mu|tau)nu']]
    # ttbar
    searches['ttbar'] = [['ttbar']]
    # single top
    searches['Single top'] = [['top','antitop','atop']]
    # Diboson
    searches['Diboson'] = [['!H125'],['WW','ZZ','WZ','ZW','[lv]{4}','[WZ][pqlv]*[WZ]']]
    # Signal: LFV Higgs
    searches['LFV Higgs'] = [['H125'],['taue','etau','mutau','taumu']]

    matched_category = []
    for category, search in searches.iteritems():
        match_results = []
        for must_match in search:
            for possible_match in must_match:
                if possible_match[0] == '!':
                    found = not re.search(r'%s'%(possible_match[1:]),sample)
                else:
                    found = re.search(r'%s'%(possible_match),sample)
                # if found, skip to next must_match
                if found:
                    match_results.append(True)
                    break
            else:
                # if not found, skip to next search
                match_results.append(False)
                break
            # if match with all must_match, store category
        if all(match_results):
            matched_category.append(category)
    if len(matched_category) == 0:
        return 'Unknown'
    elif len(matched_category) == 1:
        return matched_category[0]
    else:
        comb_str = 'Double Matched: '
        for i, category in enumerate(matched_category):
            if i:
                comb_str += ', '
            comb_str += category
        return comb_str


if __name__ == '__main__':
    main()
