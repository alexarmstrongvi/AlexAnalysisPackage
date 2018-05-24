"""
================================================================================
Plot classes for ROOT plotting

Author:
    Alex Armstrong <alarmstr@cern.ch>
    ... with almost the entire code borrowed exactly from Danny Antrim
        <dantrim@cern.ch>

License:
    Copyright: (C) <May 20th, 2018>; University of California, Irvine
================================================================================
"""
# General python
import sys
import glob
from math import floor
from enum import Enum
from collections import OrderedDict

# Root data analysis framework
import ROOT as r
r.TCanvas.__init__._creates = False

################################################################################
# Enumerating plot types
################################################################################
class Types(Enum):
    default = 0
    stack = 1
    ratio = 2
    double_ratio = 3
    comparison = 4
    undefined = 5

################################################################################
# Plot Classes
################################################################################
class Plot1D :
    def __init__(self) :
        # Descriptors
        self.region = ""
        self.name = ""
        self.variable = ""

        # Objects
        self.pads = None

        # Properties
        self.x_label = "x-Label"
        self.y_label = "Events"
        self.x_bin_width = 1.0
        self.x_range_min = 0.0
        self.x_range_max = 50.0
        self.y_range_min = 0.0
        self.y_range_max = 50.0
        self.nbins = 20

        # Flags
        self.is2D = False
        self.doLogY = False
        self.doNorm = False
        self.add_overflow = True
        self.add_underflow = False
        self.leg_is_left = False
        self.leg_is_bottom_right = False
        self.leg_is_bottom_left = False
        self.ptype = Types.undefined

    def initialize(self, region = "", variable = "", name = "") :
        self.region = region
        self.variable = variable
        self.name = name

    def labels(self, x="", y="Entries") :
        self.x_label = x
        self.y_label = y

    def xax(self, width=1.0, min = 0.0, max = 500) :
        self.x_bin_width = width
        self.x_range_min = min
        self.x_range_max = max
        # call this to reset nbins
        self.nbins = self.get_n_bins()

    def yax(self, min = 0.0, max = 500) :
        self.y_range_min = min
        self.y_range_max = max

    def setDefaultCanvas(self, name) :
        self.pads = Canvas(name)
        self.ptype = Types.default

    def setStackCanvas(self, name):
        self.pads = StackCanvas(name)
        self.ptype = Types.stack

    def setRatioCanvas(self, name) :
        self.pads = RatioCanvas(name)
        self.ptype = Types.ratio

    def setDoubleRatioCanvas(self, name) :
        self.pads = DoubleRatioCanvas(name)
        self.ptype = Types.double_ratio

    def get_n_bins(self) :
        '''
        From the user-provided bin width and (min,max) get
        the number of bins
        '''
        max = self.x_range_max
        min = self.x_range_min
        width = self.x_bin_width
        nbins = floor( (max - min) / (width) + 0.5 )
        return nbins

    def Print(self) :
        print "Plot1D    plot: %s  (region: %s  var: %s)"%(self.name, self.region, self.variable)

class Plot2D :
    def __init__(self) :
        # Descriptors
        self.region = ""
        self.name = ""

        # Objects
        self.canvas = None

        # Properties
        self.style = "colz"
        self.xVariable = ""
        self.yVariable = ""
        self.x_label = "x-Label"
        self.y_label = "y-Label"
        self.x_bin_width = 1.0
        self.x_range_min = 0.0
        self.x_range_max = 50.0
        self.y_range_min = 0.0
        self.y_range_max = 50.0
        self.y_bin_width = 1.0

        # Flags
        self.do_profile = False
        self.do_profileRMS = False

    def initialize(self, region="", xvar="", yvar="", name="") :
        '''
        Initalize the selection ('region'), x- and y-variables
        to be plotted, as well as the name of the plot
        '''
        self.region = region
        self.xVariable = xvar
        self.yVariable = yvar
        self.name = name

    def setDefaultCanvas(self, name) :
        c = r.TCanvas("c_"+name, "c_"+name, 800, 600)
        self.canvas = c

    def labels(self, x="",y="") :
        '''
        Set the x- and y-axis titles
        '''
        self.x_label = x
        self.y_label = y

    def xax(self, width=1.0, min=0.0, max=50.0) :
        '''
        Set the x-axis attributes
        '''
        self.x_bin_width = width
        self.x_range_min = min
        self.x_range_max = max
        self.n_binsX = self.get_n_bins(width, min, max)

    def yax(self, width=1.0, min=0.0, max=50.0) :
        '''
        Set the y-axis attributes
        '''
        self.y_bin_width = width
        self.y_range_min = min
        self.y_range_max = max
        self.n_binsY = self.get_n_bins(width, min, max)

    def doProfile(self) :
        '''
        Set whether to do a profile plot
        '''
        self.do_profile = True

    def doProfileRMS(self) :
        '''
        Set whether to do a profile plot
        with RMS on the y-axis
        '''
        self.do_profileRMS = True

    def get_n_bins(self, width, min, max) :
        '''
        From the user-provided bin width and (min,max) get
        the number of bins for the specified axis
        '''
        nbins = floor( (max - min) / (width) + 0.5 )
        return nbins

    def set_style(self, style="") :
        '''
        Override the default style of "colz"
        '''
        self.style = style

    def Print(self) :
        print "Plot2D    plot: %s  (region: %s  xVar: %s  yVar: %s)"%(self.name, self.region, self.xVariable, self.yVariable)


################################################################################
# TPad handler classes
################################################################################
#TODO: rename Canvas -> Pad based names
class Canvas :
    def __init__(self,name):
        self.name = "c_" + name
        self.canvas = r.TCanvas(self.name, self.name, 800, 600)
        self.set_pad_dimensions()

    def set_pad_dimensions(self):
        pass

class StackCanvas(Canvas):
    def __init__(self, name):
        Canvas.__init__(self, name)

    def set_pad_dimensions(self):
        can = self.canvas
        can.cd()

        # Color
        can.SetFrameFillColor(0)
        can.SetFillColor(0)

        # Margins
        can.SetRightMargin(0.05)
        can.SetLeftMargin(0.14)
        can.SetBottomMargin(1.3*can.GetBottomMargin())

        can.Update()
        self.canvas = can

class RatioCanvas :
    def __init__(self, name) :
        self.name = "c_" + name
        self.canvas = r.TCanvas(self.name, self.name, 768, 768)
        self.upper_pad = r.TPad("upper", "upper", 0.0, 0.0, 1.0, 1.0)
        self.lower_pad = r.TPad("lower", "lower", 0.0, 0.0, 1.0, 1.0)
        self.set_pad_dimensions()

    def set_pad_dimensions(self) :
        can = self.canvas
        up  = self.upper_pad
        dn  = self.lower_pad

        can.cd()
        up_height = 0.75
        dn_height = 0.30
        up.SetPad(0.0, 1.0-up_height, 1.0, 1.0)
        dn.SetPad(0.0, 0.0, 1.0, dn_height)

        up.SetTickx(0)
        dn.SetGrid(0)
        dn.SetTicky(0)

        up.SetFrameFillColor(0)
        up.SetFillColor(0)

        # set margins
        up.SetRightMargin(0.05)
        up.SetLeftMargin(0.14)
        up.SetTopMargin(0.7 * up.GetTopMargin())
        up.SetBottomMargin(0.09)

        dn.SetRightMargin(up.GetRightMargin())
        dn.SetLeftMargin(up.GetLeftMargin())
        dn.SetBottomMargin(0.4)

        up.Draw()
        dn.Draw()
        can.Update()

        self.canvas = can
        self.upper_pad = up
        self.lower_pad = dn

class DoubleRatioCanvas :
    def __init__(self, name) :
        self.name = "c_" + name
        self.canvas = r.TCanvas(self.name, self.name, 300, 350)
        self.upper_pad = r.TPad("upper", "upper", 0.0, 0.0, 1.0, 1.0)
        self.middle_pad = r.TPad("middle", "middle", 0.0, 0.0, 1.0, 1.0)
        self.lower_pad = r.TPad("lower", "lower", 0.0, 0.0, 1.0, 1.0)
        self.set_pad_dimensions()

    def set_pad_dimensions(self) :
        can = self.canvas
        up  = self.upper_pad
        mid = self.middle_pad
        dn = self.lower_pad

        can.cd()
        up_height = 0.90
        mid_height_low = 0.25
        mid_height_high = 0.40
        dn_height = 0.25

        up.SetPad(0.0, mid_height_high, 1.0, 1.0)
        mid.SetPad(0.0, mid_height_low, 1.0, mid_height_high)
        dn.SetPad(0.0, 0.0, 1.0, mid_height_low)

        up.SetTickx(0)
        mid.SetGrid(0)
        mid.SetTicky(0)
        dn.SetGrid(0)
        dn.SetTicky(0)

        up.SetFrameFillColor(0)
        up.SetFillColor(0)

        # set right margins
        right_margin = 0.05
        up .SetRightMargin(right_margin)
        mid.SetRightMargin(right_margin)
        dn .SetRightMargin(right_margin)

        # set left margins
        left_margin = 0.14
        up .SetLeftMargin(left_margin)
        mid.SetLeftMargin(left_margin)
        dn .SetLeftMargin(left_margin)

        # bottom margins
        up.SetBottomMargin(0.04)
        mid.SetBottomMargin(0.15)
        dn.SetBottomMargin(0.47)

        # set top margins
        up.SetTopMargin(0.09)
        mid.SetTopMargin(0.05)
        dn.SetTopMargin(0.02)



        up.Draw()
        mid.Draw()
        dn.Draw()
        can.Update()

        self.canvas = can
        self.upper_pad = up
        self.middle_pad = mid
        self.lower_pad = dn

