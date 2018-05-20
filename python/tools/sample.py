"""
================================================================================
Class for handleing samples for plotting

Author:
    Alex Armstrong <alarmstr@cern.ch>
    ... with lots of ideas borrowed from Danny Antrim <dantrim@cern.ch>

License:
    Copyright: (C) <May 20th, 2018>; University of California, Irvine
================================================================================
"""
# General python
import sys, os
import glob

# Root data analysis framework
import ROOT as r
r.TColor.__init__._creates = False
r.TEventList.__init__._creates = False

################################################################################
# Main Base Class
################################################################################
class Sample :
    def __init__(self, name = "", displayname = ""):
        self.dbg = False #TODO: Remove debug
        self.name = name
        self.displayname = displayname
        self.treename = ""
        self.color = r.kRed
        self.tree = None
        self.isMC = False

    def isMC(self):
        return self.isMC

    # Setup TChain/TTree
    def set_chain_from_root_file(self, file_name, flat_ntuple_dir):
        full_name = flat_ntuple_dir+file_name
        self.set_chain_from_list([full_name])

    def set_chain_from_dsid_list(self, dsid_list, flat_ntuple_dir, search_strs='', exclude_strs=''):
        '''
        Build chain of background ntuples

        Args:
            dsid_list (list(int)): list of DSIDs to include in TChain
            flat_ntuple_dir (str): Name of dir containing all .root flat ntuples

        Returns:
            TChain: TChain of flat ntuples from input directory
        '''
        # Get list of flat ntuples file names from sample directory
        flat_ntuples = glob.glob(flat_ntuple_dir + "*.root")
        assert len(flat_ntuples), "No root files found at %s"%flat_ntuple_dir

        # Select out flat ntuples found in DSID list
        chosen_ntuples = []
        #TODO: Make inputs lists instead of converting here

        search_str_list = [s.strip() for s in search_strs.split(",")]
        exclude_str_list = [s.strip() for s in exclude_strs.split(",")]
        for dsid in dsid_list:
            for fname in flat_ntuples:
                if any(s not in fname for s in search_str_list if s): continue
                if any(s in fname for s in exclude_str_list if s) : continue
                if str(dsid) in fname :
                    chosen_ntuples.append(fname);
                    break
            else:
                print "WARNING :: Unable to find file for DSID =", dsid
        if not chosen_ntuples:
            print "WARNING :: No samples found for", self.name
        
        self.set_chain_from_list(chosen_ntuples)

    def set_chain_from_file_list(self, file_list, flat_ntuple_dir):
        full_file_list = [flat_ntuple_dir+f for f in file_list]
        self.set_chain_from_list(full_file_list)

    def set_chain_from_list(self, files):
        chain = r.TChain("superNt") #TODO: Remove hard-coded value
        for n_files, fname in enumerate(files):
            chain.Add(fname)
        print "%10s : ADDED %d FILES"%(self.name, n_files+1)
        self.tree = chain

    # Comparison
    def __eq__(self, other) :
        return (self.displayname == other.displayname
            and self.name == other.name
            and self.treename == other.treename)

     # Comparison
    def Print(self) :
        print 'Sample "%s" (tree %s)'%(self.displayname, self.treename)

    def __gt__(self, other, tcut) :
        '''
        Comparison operator to order background samples by
        their yields in a given region defined by tcut
        '''
        cut = r.TCut(tcut)
        sel = r.TCut("1")
        return ( self.tree.Draw("isMC", cut * sel, "goff") > other.tree.Draw("isMC", cut * sel, "goff") )

    # Sanity checks
    def check_for_duplicates(self):
        if not self.tree:
            print "ERROR (sample.CheckForDuplicates) :: tree not yet defined"
            return
        print "Checking for duplicate events in", self.name
        events = [x.event_number for x in self.tree]
        n_events = len(events)
        n_events_no_dup = len(set(events))
        if n_events_no_dup != n_events:
            dup_evts = n_events - n_events_no_dup
            print "There are %d/%d duplicate events"%(dup_evts, n_events)

################################################################################
# Data class
################################################################################
class Data(Sample):
    def __init__(self) :
        Sample.__init__(self, "Data", "Data")
        self.color = r.kBlack
        self.isMC = False

    def Print(self) :
        print 'Data (tree %s)'%(self.treename)

################################################################################
# MC classes
################################################################################
class MCsample(Sample):
    def __init__(self, name = "", displayname = ""):
        Sample.__init__(self, name, displayname)
        self.dsid = ""
        self.line_style = 1
        self.fill_style = 0
        self.scale_factor = 1.0
        self.isSignal = False
        self.isMC = True

    def isSignal(self) :
        return self.is_signal
    def Print(self) :
        print 'MC Sample "%s" (tree %s)'%(self.displayname, self.treename)

class Background(MCsample) :
    def __init__(self, name = "", displayname = "") :
        Sample.__init__(self, name, displayname)
        self.is_signal = False

    def Print(self) :
        print 'Background "%s" (tree %s)'%(self.displayname, self.treename)

class Signal(MCsample) :
    def __init__(self, name = "", displayname ="") :
        Sample.__init__(self, name, displayname)
        self.is_signal = True

    def Print(self) :
        print 'Signal "%s" (tree %s)'%(self.displayname, self.treename)


