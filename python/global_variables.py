#!/bin/bash/env python
import ROOT

# User-defined analysis specific settings
tag = "n0235"
analysis_path = '/data/uclhc/uci/user/armstro1/SusyNt/SusyNt_n0235_LFV_analysis/' 
LOCAL_DSID_DIR = '/data/uclhc/uci/user/dantrim/susynt_productions/n0235/'
LOCAL_DSID_SUBDIRS = {'data':'data','mc':'mc'}

analysis_dir = analysis_path+'analysis/'
analysis_run_dir = analysis_path+'analysis_run/'

fax_input_files     = analysis_dir+'inputs/'
local_input_files   = analysis_dir+'inputs_local/'
output_dir          = analysis_run_dir+'outputs/'
logs_dir            = analysis_run_dir+'logs/'
dsid_dir            = analysis_run_dir+'dsid_filelist/'
missing_dsids_file  = output_dir+'missing.txt'

event_list_dir      = analysis_run_dir+"lists/"
plots_dir           = analysis_run_dir+"plots/"

mc_ntuples = analysis_run_dir+"ntuples/mc/"
data_ntuples = analysis_run_dir+"ntuples/data/"

input_files = local_input_files


# Samples recommended by Duc Bao Ta (xiedebao@googlemail.com)
# subset of full sample list used in Htt analysis (ATL-COM-PHYS-2017-446)
# additions where indicated
groups = {}
groups['higgs_lfv'] = [
    344084, 344085, 344086, 344087, 344088, 344089, 344090, 344091,
    345077, 345078, 
    345124, 345125, 
    345213, 345214, 345215, 345216,
    345218, 345219
]
groups['htt'] = [
    341122, 
    341155,
    342170, 342171, 342172,
    345073,
    345120,
    #343323,
    #343334,
    #343942,
    #343953,
]
groups['hww'] = [
    341079, 341080,
]
groups['ttbar'] = [
    410009, 
]
groups['singletop'] = [
    410011, 410012, 
    410015, 410016, 
    410025, 410026,
]
groups['ww'] = [
    361600,
]
groups['wz'] = [
    361601,
    361607,
]
groups['zz'] = [
    361603, 361604, 
    361610,
]
groups['zll'] = [
    364100, 364101, 364102, 364103, 364104, 364105, 364106, 364107, 364108,
    364109, 364110, 364111, 364112, 364113,
    364114, 364115, 364116, 364117, 364118, 364119, 364120, 364121, 364122,
    364123, 364124, 364125, 364126, 364127,
    364198, 364199, 364200, 364201, 364202, 364203,
    364204, 364205, 364206, 364207, 364208, 364209,
    #304018, 304019, 304020, replaced with 344441-3
    344441, 344442, 344443,
    345099, 345100, 345101, 345102,
]
groups['ztt'] = [
    344772, 344776, 344780,
    #304021,
    #364128, 364129, 364130, 364131, 364132, 364133, 364134, 364135, 364136,
    364137, 364138, 364139, 364140, 364141, 364210, 364211, 364212, 364213,
    364214, 364215,
] 
groups['wjets'] = [
    # Temporary fake sample
    361100, 361101, 361102, 361103, 361104, 361105,
]
groups['data15'] = [
    276262, 276329, 276336, 276416, 276511, 276689, 276778, 276790, 276952, 
    276954, 278880, 278912, 278968, 279169, 279259, 279279, 279284, 279345, 
    279515, 279598, 279685, 279813, 279867, 279928, 279932, 279984, 280231, 
    280273, 280319, 280368, 280423, 280464, 280500, 280520, 280614, 280673, 
    280753, 280853, 280862, 280950, 280977, 281070, 281074, 281075, 281317, 
    281385, 281411, 282625, 282631, 282712, 282784, 282992, 283074, 283155, 
    283270, 283429, 283608, 283780, 284006, 284154, 284213, 284285, 284420, 
    284427, 284484,
]
groups['data16'] = [
    297730, 298595, 298609, 298633, 298687, 298690, 298771, 298773, 298862, 
    298967, 299055, 299144, 299147, 299184, 299243, 299584, 300279, 300345, 
    300415, 300418, 300487, 300540, 300571, 300600, 300655, 300687, 300784, 
    300800, 300863, 300908, 301912, 301918, 301932, 301973, 302053, 302137, 
    302265, 302269, 302300, 302347, 302380, 302391, 302393, 302737, 302831, 
    302872, 302919, 302925, 302956, 303007, 303079, 303201, 303208, 303264, 
    303266, 303291, 303304, 303338, 303421, 303499, 303560, 303638, 303832, 
    303846, 303892, 303943, 304006, 304008, 304128, 304178, 304198, 304211, 
    304243, 304308, 304337, 304409, 304431, 304494, 305380, 305543, 305571, 
    305618, 305671, 305674, 305723, 305727, 305735, 305777, 305811, 305920, 
    306269, 306278, 306310, 306384, 306419, 306442, 306448, 306451, 307126, 
    307195, 307259, 307306, 307354, 307358, 307394, 307454, 307514, 307539, 
    307569, 307601, 307619, 307656, 307710, 307716, 307732, 307861, 307935, 
    308047, 308084, 309375, 309390, 309440, 309516, 309640, 309674, 309759, 
    310015, 310247, 310249, 310341, 310370, 310405, 310468, 310473, 310634, 
    310691, 310738, 310809, 310863, 310872, 310969, 311071, 311170, 311244, 
    311287, 311321, 311365, 311402, 311473, 311481,
]

def get_mc_groups():
    mc_groups = {}
    for key, lst in groups.iteritems():
        if 'data' in key.lower():
            continue
        mc_groups[key] = lst
    return mc_groups
def get_bkgd_groups():
    bkgd_groups = {}
    for key, lst in groups.iteritems():
        if 'data' in key.lower() or 'signal' in key.lower():
            continue
        bkgd_groups[key] = lst
    return bkgd_groups
def get_data_groups():
    data_groups = {}
    for key, lst in groups.iteritems():
        if 'data' not in key.lower():
            continue
        data_groups[key] = lst
    return data_groups

# Set ATLAS style
def setATLASStyle(path="/data/uclhc/uci/user/armstro1/ATLASStyle/"):
    ROOT.gROOT.SetMacroPath(path)
    ROOT.gROOT.LoadMacro("AtlasStyle.C")
    ROOT.gROOT.LoadMacro("AtlasLabels.C")
    ROOT.gROOT.LoadMacro("AtlasUtils.C")
    ROOT.SetAtlasStyle()
