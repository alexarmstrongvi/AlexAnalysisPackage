import ROOT as r
import sys
sys.path.append('../..')

r.TH1F.__init__._creates = False
r.TTree.__init__._creates = False

class Systematic :
    def __init__(self, syst_name="",up_name="",down_name="") :

        # these are the names that will
        # be used to define the treenames
        # i.e.
        #    <BKG>_<name><up_name>
        #  --> Zjets_MSUP
        # if weight systematic these will
        # be the names appended to "syst_"
        # in the CENTRAL tree
        # i.e.
        #    syst_<name><up_name>
        self.name = syst_name
        self.up_name = ""
        self.down_name = ""

        # toggles to set whether to treat as a weight
        # or object/kin sys
        self.weight_sys = False
        self.kin_sys = False
        self.one_side = False

        # placeholders for various
        # trees
        self.tree = None
        self.tree_up = None
        self.tree_down = None

        self.up_histo = None
        self.down_histo = None

        self.up_yield = 0.0
        self.down_yield = 0.0

    def setWeightSys(self) :
        '''
        Toggle whether this systematic is
        a weight systematic
        '''
        self.weight_sys = True
        self.name = "syst_" + self.name
        self.up_name = self.name + "UP"
        self.down_name = self.name + "DOWN"
        #if "syst_" not in self.name :
        #    self.name = "syst_" + self.name
        #    self.up_name = self.name + self.up_name
        #    self.down_name = self.name + self.down_name

    def isWeightSys(self) :
        return self.weight_sys

    def setKinSys(self) :
        '''
        Toggle whether this systematic is
        a kinematic/object systematic
        '''
        self.kin_sys = True
        if "SoftTrk_Scale" not in self.name :
            self.up_name = self.name + "_UP"
            self.down_name = self.name + "_DN"
        else :
            self.up_name = self.name + "Up"
            self.down_name = self.name + "Down"

    def isKinSys(self) :
        return self.kin_sys

    def setOneSided(self) :
        '''
        Set this systematic as onesided
        '''
        self.one_side = True

    def isOneSided(self) :
        return self.one_side

    def setUpYield(self, up_yield) :
        self.up_yield = up_yield

    def upYield(self) :
        return self.up_yield

    def setDownYield(self, down_yield) :
        self.down_yield = down_yield

    def downYield(self) :
        return self.down_yield

    def check(self) :
        if self.weight_sys and self.kin_sys :
            print "Systematic::check ERROR    you have configured a systematic (%s) as both a weight systematic AND kinematic systematic. Exitting."%self.name
            sys.exit()

    def Print(self) :
        desc = ""
        if self.weight_sys : desc = "weights"
        if self.kin_sys : desc = "kinematics"
        print "Systematic -- %s -- %s (%s, %s)"%(desc, self.name, self.up_name, self.down_name)

 
