import ROOT as r
import glob
import sys
sys.path.append('../..')

import supersusy.utils.background as background

class Signal() :
    def __init__(self, name = "", displayname ="") :
        self.dbg = False
        self.file = ""
        self.grid = ""
        #self.dsid = str(self.dsid_from_file(file))
        self.dsid = 0
        self.tree = None

        self.color = r.kBlue
        self.displayname = ""
        self.scale_factor = 1.0
        
        # mass info
        self.mA = 0.0
        self.mB = 0.0
        self.mC = 0.0
    
    def __eq__(self, other) :
        '''
        Equality method to compare two signal points to test
        if they are the same
        '''
        return (self.mA==other.mA) and (self.mB==other.mB) and (self.mC==other.mC) 
   
    def set_debug(self) :
        self.dbg = True
    def set_file(self, file) :
        self.file = file    
    def set_grid(self, grid) :
        self.grid = grid

    def set_dsid_from_file(self, file) :
        dsid = ""
        if "CENTRAL" in file :
            dsid = file[file.find('CENTRAL_')+8 : self.file.find('.root')]
        else :
            print "Signal.set_dsid_from_file ERROR    File format for signal sample is unexpected. Exitting."
            sys.exit()
        self.dsid = dsid

    def set_tree(self) :
        tree_name = ""
        if "CENTRAL" in self.file :
            tree_name = 'id_' + str(self.dsid)
        else :
            print "Signal.set_tree ERROR    File format for signal sample is unexpected. Exitting."
            sys.exit()
        infile = r.TFile.Open(self.file)
        chain = r.TChain(tree_name)
        chain.Add(self.file)
        self.tree = chain

    def set_mass_info(self) :
        if self.grid=="tN1" or self.grid=="bWN" :
            txtfile="../susyinfo/tN1/masses_tN1.txt"
            lines = open(txtfile).readlines()
            fields = lines[0].split()
            dsIdx = fields.index('DSID')
            mtIdx = fields.index('mt[GeV]')
            mn1Idx = fields.index('mN1[GeV]')
            for line in lines[1:] :
                line = line.strip()
                if not line : continue
                fields = line.split()
                if fields[0]!=str(self.dsid) : continue
                self.mA = float(fields[mtIdx])
                self.mB = float(fields[mn1Idx])

                if self.dbg : print "Signal %s found at (m_st, m_c1)=(%s, %s)"%(self.dsid, self.mA, self.mB)
        else :
            print "Signal.get_mass_info ERROR    Requested grid (%s) not supported. Exitting."%(self.grid)
            sys.exit()

    def set_displayname_from_masses(self) :
        if (self.mA==0.0) and (self.mB==0.0) and (self.mC==0.0) :
            print "Signal.set_displayname_from_masses ERROR    You must set the signal point masses before calling this function! Exitting."
            sys.exit()

        if self.grid=="tN1" or self.grid=="bWN" :
            self.displayname = "(%.0f, %.0f)"%(float(self.mA), float(self.mB)) 
        else :
            print "Signal.set_displayname_from_masses ERROR    Requested grid (%s) not supported. Exitting."%(self.grid)
            sys.exit()
                
        
