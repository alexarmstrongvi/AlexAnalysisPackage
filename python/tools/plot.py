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
import re
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
# TODO add some generic functions to plot utils

class Plot1D(object) :

    # Class defaults for axis range depending on log and normalization settings
    xmin = 0
    xmax = 50.0
    ymin = 0
    ymax =1e4

    logy_min = 1e-1
    logy_max = 1e10

    norm_ymin = 0
    norm_ymax = 1.5

    logy_norm_min = 1e-7
    logy_norm_max = 1e4

    auto_set_ylimits = True

    def __init__(self,
        region = "",
        name = "",
        variable = "",
        xlabel = "x-Label",
        xunits = "",
        ylabel = "",
        yunits = "",
        bin_range = [], # [x0, x1, y0, y1]
        bin_width = None, # Can specify to override nbins
        nbins = 20,
        xmin = None,
        xmax = None,
        ymin = None,
        ymax = None,
        logy_min = None,
        logy_max = None,
        norm_ymin = None,
        norm_ymax = None,
        logy_norm_min = None,
        logy_norm_max = None,
        is2D = False,
        doLogY = True, # change to do_logy
        doNorm = False,
        add_overflow = True,
        add_underflow = False,
        leg_is_left = False,
        leg_is_bottom_right = False,
        leg_is_bottom_left = False,
        ptype = Types.default,
        ) :
        '''
        Constructor

        Order is important so be careful when rearranging.
        '''
        # Descriptors
        self.region = region
        self.variable = variable
        self.name = name if name else self.determine_name(region, variable)

        # Objects
        self.pads = None

        # Flags
        self.is2D = is2D
        self.doLogY = doLogY
        self.doNorm = doNorm

        self.add_overflow = add_overflow
        self.add_underflow = add_underflow
        self.leg_is_left = leg_is_left
        self.leg_is_bottom_right = leg_is_bottom_right
        self.leg_is_bottom_left = leg_is_bottom_left
        self.ptype = ptype

        # Properties
        if xmin: self.xmin = xmin
        if xmax: self.xmax = xmax
        if ymin: self.ymin = ymin
        if ymax: self.ymax = ymax
        if logy_min: self.logy_min = logy_min
        if logy_max: self.logy_max = logy_max
        if norm_ymin: self.norm_ymin = norm_ymin
        if norm_ymax: self.norm_ymax = norm_ymax
        if logy_norm_min: self.logy_norm_min = logy_norm_min
        if logy_norm_max: self.logy_norm_max = logy_norm_max

        self.xmin, self.xmax, self.ymin, self.ymax = self.determine_range(bin_range)

        self.nbins = self.determine_nbins(bin_width) if bin_width else nbins

        self.xunits = xunits
        self.yunits = yunits
        self.xlabel, self.ylabel = self.determine_labels(xlabel, ylabel)

        self.bin_labels = []
        self.rebin_bins = []

    def update(self,
        region = None,
        variable = None,
        name = None,
        bin_range = None,
        bin_width = None,
        nbins = None,
        xmin = None,
        xmax = None,
        ymin = None,
        ymax = None,
        logy_min = None,
        logy_max = None,
        norm_ymin = None,
        norm_ymax = None,
        logy_norm_min = None,
        logy_norm_max = None,
        doLogY = None,
        doNorm = None,
        add_overflow = None,
        add_underflow = None):
        '''
        Update plot properties

        This is mainly for updating plot properties of a cloned plot. Therefore,
        only a minimum set of properties that one expects could be updated are
        included.
        '''
        if region: self.region = region
        if variable: self.variable = variable
        if name:
            self.name = name
        elif not name and (region or variable):
            self.name = self.determine_name(region, variable)

        # Flags
        if doLogY: self.doLogY = doLogY
        if doNorm: self.doNorm = doNorm
        if add_overflow: self.add_overflow = add_overflow
        if add_underflow: self.add_underflow = add_underflow

        # Properties
        if xmin: self.xmin = xmin
        if xmax: self.xmax = xmax
        if ymin: self.ymin = ymin
        if ymax: self.ymax = ymax
        if logy_min: self.logy_min = logy_min
        if logy_max: self.logy_max = logy_max
        if norm_ymin: self.norm_ymin = norm_ymin
        if norm_ymax: self.norm_ymax = norm_ymax
        if logy_norm_min: self.logy_norm_min = logy_norm_min
        if logy_norm_max: self.logy_norm_max = logy_norm_max

        self.xmin, self.xmax, self.ymin, self.ymax = self.determine_range(bin_range)
        if nbins or bin_width:
            self.nbins = self.determine_nbins(bin_width) if bin_width else nbins



    def determine_nbins(self, bin_width, update_range = True):
        '''
        Intelligently determine the number of bins.

        Given a desired bin width and x-axis range, the correct number of bins
        is determined for use with TH1. The x-axis range or bin width will
        likely be adjusted so that the latter divides the former. If the
        boundaries are int type values, these will be maintained.

        args:
            bin_width (float) - desired width of x-axis bins
            update_range - option to update x-axis range to be divisible by
                bin width as opposed to modying bin width to be a divisor

        returns:
            (int) - the number x-axis bins
        '''
        assert self.xmax != self.xmin, ("ERROR :: x-axis range not set")
        bin_width = float(bin_width)
        x_range = float(self.xmax - self.xmin)
        remainder = x_range % bin_width
        cutoff_range = bin_width - remainder if remainder else 0
        int_bins = (isinstance(self.xmin, (int, long))
                and isinstance(self.xmax, (int, long)))

        eps = 0.0000001
        notify = cutoff_range > eps and abs(cutoff_range - bin_width) > eps
        if update_range and notify:
            print "WARNING :: bin_width is not a divisor of x-axis range.",
            print "Variable = %s, Range = [%f,%f], bin_width = %f"%(
                    self.variable, self.xmin, self.xmax, bin_width)
            print "Expanding range by %f to fit"%cutoff_range
        if update_range and self.xmin == 0:
            self.xmax += cutoff_range
        elif update_range and self.xmin != 0:
            self.xmin -= 0.5 * cutoff_range
            self.xmax += 0.5 * cutoff_range

        if int_bins:
            self.xmin = int(round(self.xmin))
            self.xmax = int(round(self.xmax))

        x_range = (self.xmax - self.xmin)
        nbins = round( x_range / bin_width )

        return nbins

    def bin_width(self):
        return (self.xmax - self.xmin) / self.nbins

    def determine_range(self, bin_range):
        '''
        Inteligently determine bin range

        The user can provide 2 or 4 values in a list. The default range values
        will be modified. It is assumed providing two values indicates one only
        wants to set the x-axis range.

        args:
            bin_range (list(int or float)) - list of bin range values

        returns:
            tuple (4 int or float) - range values for both x- and y-axis

        '''
        # default values

        assert not bin_range or len(bin_range) in [2,4],(
            'ERROR :: Unrecognized bin range format:', bin_range)


        if self.doLogY and self.doNorm:
            ymin = self.logy_norm_min
            ymax = self.logy_norm_max
        elif self.doLogY:
            ymin = self.logy_min
            ymax = self.logy_max
        elif self.doNorm:
            ymin = self.norm_ymin
            ymax = self.norm_ymax
        else:
            ymin = self.ymin
            ymax = self.ymax

        xmin = self.xmin
        xmax = self.xmax

        if not bin_range:
            return  xmin, xmax, ymin, ymax
        elif len(bin_range) == 2:
            return bin_range[0], bin_range[1], ymin, ymax
        elif len(bin_range) == 4:
            return bin_range[0], bin_range[1], bin_range[2], bin_range[3]

    def determine_labels(self, xlabel, ylabel):

        # Set x-axis label
        if self.xunits:
            xlabel = "%s [%s]"%(xlabel, self.xunits)

        # set y-axis label defaults
        if ylabel:
            pass
        elif self.doNorm:
            ylabel = "a.u."
        else:
            ylabel = "Events"

        # set y-axis label
        width_label = str(round(self.bin_width(), 2))
        if not self.xunits and width_label == '1.0':
            pass
        elif self.xunits and width_label == '1.0':
            ylabel = "%s / %s"%(ylabel, self.xunits)
        else:
            ylabel = "%s / %s %s"%(ylabel,width_label,self.xunits)

        return xlabel, ylabel

    def determine_name(self, region, variable):
        var_stripped = re.sub(r'[(){}[\]]+','', variable)
        return "%s_%s"%(region, var_stripped)

    def set_bin_labels(self, axis_hist):
        '''
        Add arbitrary labels to TH1 bins

        Set labels for a TH1 histogram after defining bin_labels. It is
        required that the number of labels be the same as the number of bins. If
        one wants only certain bins to be labeled, then put empty strings for
        those indices and they will be skipped. Put a blank string to create
        empty labels

        args:
            axis_hist (TH1): histogram that includes plot axis

        '''
        n_labels = len(self.bin_labels)
        n_bins = axis_hist.GetNbinsX()
        assert n_labels == n_bins, (
            "ERROR :: nbins (%d) != nlabels (%d)"%(n_bins, n_labels))

        x_axis = axis_hist.GetXaxis()
        for ibin, label in zip(range(n_labels), self.bin_labels):
            if not label: continue
            x_axis.SetBinLabel(ibin+1, label)
        x_axis.SetLabelOffset(0.005)
        x_axis.SetLabelSize(1.5*x_axis.GetLabelSize())

    def setDefaultPads(self, name) :
        self.pads = Pads(name)
        self.ptype = Types.default

    def setStackPads(self, name):
        self.pads = StackPads(name)
        self.ptype = Types.stack

    def setRatioPads(self, name) :
        self.pads = RatioPads(name)
        self.ptype = Types.ratio

    def setDoubleRatioPads(self, name) :
        self.pads = DoubleRatioPads(name)
        self.ptype = Types.double_ratio

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
        self.xlabel = "x-Label"
        self.ylabel = "y-Label"
        self.x_bin_width = 1.0
        self.xmin = 0.0
        self.xmax = 50.0
        self.ymin = 0.0
        self.ymax = 50.0
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

    def setDefaultPads(self, name) :
        c = r.TCanvas("c_"+name, "c_"+name, 800, 600)
        self.canvas = c

    def labels(self, x="",y="") :
        '''
        Set the x- and y-axis titles
        '''
        self.xlabel = x
        self.ylabel = y

    def xax(self, width=1.0, min=0.0, max=50.0) :
        '''
        Set the x-axis attributes
        '''
        self.x_bin_width = width
        self.xmin = min
        self.xmax = max
        self.n_binsX = self.get_n_bins(width, min, max)

    def yax(self, width=1.0, min=0.0, max=50.0) :
        '''
        Set the y-axis attributes
        '''
        self.y_bin_width = width
        self.ymin = min
        self.ymax = max
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
class Pads :
    def __init__(self, name):
        self.name = "c_" + name
        self.canvas = r.TCanvas(self.name, self.name, 800, 600)
        self.set_pad_dimensions()

    def set_pad_dimensions(self):
        pass

class StackPads(Pads):
    def __init__(self, name):
        Pads.__init__(self, name)

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

class RatioPads :
    ylabel = 'Data / MC'
    ymax = 2
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

class DoubleRatioPads :
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

