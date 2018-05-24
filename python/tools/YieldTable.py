"""
================================================================================
Class for storing and manipulating yields

Author:
    Alex Armstrong <alarmstr@cern.ch>

License:
    Copyright: (C) <May 20th, 2018>; University of California, Irvine
================================================================================
"""

# General python
import sys, os
import re
from numbers import Number
from math import sqrt
from collections import OrderedDict

class UncFloat :
    precision = 2

    def __init__(self, value = 0 , uncertainty = 0):
        if isinstance(value, Number):
            pass # Expected
        elif isinstance(value, basestring):
            value, uncertainty = parse_uncertainty_string(value)
        else:
            print "WARNING : Unexpected input type"

        self.value = float(value)
        self.uncertainty = float(uncertainty)

    def __add__(self, other):
        value = self.value + other.value
        uncertainty = sqrt(self.uncertainty**2 + other.uncertainty**2)
        return UncFloat(value, uncertainty)
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)
    def __sub__(self, other):
        value = self.value - other.value
        uncertainty = sqrt(self.uncertainty**2 + other.uncertainty**2)
        return UncFloat(value, uncertainty)
    def __mul__(self, other):
        value = self.value * other.value
        rel_unc1 = self.uncertainty /abs(self.value)
        rel_unc2 = other.uncertainty / abs(other.value)
        rel_unc = sqrt(rel_unc1**2 + rel_unc2**2)
        uncertainty = rel_unc * abs(value)
        return UncFloat(value, uncertainty)
    def __div__(self, other):
        value = self.value / other.value
        rel_unc1 = self.uncertainty /abs(self.value)
        rel_unc2 = other.uncertainty / abs(other.value)
        rel_unc = sqrt(rel_unc1**2 + rel_unc2**2)
        uncertainty = rel_unc * abs(value)
        return UncFloat(value, uncertainty)
    def __lt__(self, other) :
        return self.value < other.value
    def __le__(self, other) :
        return self.value <= other.value
    def __eq__(self, other) :
        return (isinstance(self, other)
            and self.value == other.value
            and self.uncertainty == other.uncertainty)
    def __ne__(self, other) :
        return not self.value == other.value
    def __gt__(self, other) :
        return self.value > other.value
    def __ge__(self, other) :
        return self.value >= other.value

    def __str__(self):
        return "%.*f +/- %.*f"%(self.precision, self.value,
                                self.precision, self.uncertainty)

    def parse_uncertainty_string(string):
        if '+/-' not in string:
            print "WARNING : Unrecognized uncertainty format:", string
            return 0, 0

        nums = [n.strip() for n in string.split('+/-')]
        nums = (float(n) for n in nums if self.is_number(n))
        if len(nums) != 2:
            print "WARNING : Unrecognized uncertainty format:", string
            return 0, 0

        return nums[0], nums[1]

    def is_number(string):
        try:
            float(string)
            return True
        except ValueError:
            return False

class YieldTable :
    precision = 2

    def __init__(self):
        self.mc = OrderedDict()
        self.data = {}
        self.signals = {}
        self.data_mc_ratio = True
        self.formulas = {}
        self.region = ""
        self.variable = ""

    def all_yields(self):
        all_yields = self.mc.copy()
        all_yields.update(self.data)
        all_yields.update(self.signals)
        return all_yields

    def all_names(self):
        all_names = []
        all_names += [k for k in self.mc]
        all_names += [k for k in self.signals]
        all_names += [k for k in self.data]
        all_names += [k for k in self.formulas]
        return all_names

    def __eq__(self, other):
        ''' Check if all yields, both value and uncertainty, are identical'''
        if not isinstance(other, self.__class__): return false

        y1 = self.all_yields()
        y2 = other.all_yields()
        for (name1, yld1), (name2, yld2) in zip(y1.iteritems(),y2.iteritems()):
            if name1 != name2 or str(yld1) != str(yld2):
                return False
        return True

    def is_empty(self):
        return not (self.mc or self.data or self.signals)

    def Print(self, no_uncertainty = False):
        if no_uncertainty:
            bkg_strings = [str(v.value) for _, v in self.mc.iteritems()]
            sig_strings = [str(v.value) for _, v in self.signals.iteritems()]
            data_string = [str(v.value) for _, v in self.data.iteritems()]
            mc_total = sum([v.value for _, v in self.mc.iteritems()])
            data_total = sum([v.value for _, v in self.data.iteritems()])
        else:
            bkg_strings = [str(v) for _, v in self.mc.iteritems()]
            sig_strings = [str(v) for _, v in self.signals.iteritems()]
            data_string = [str(v) for _, v in self.data.iteritems()]
            mc_total = sum([v for _, v in self.mc.iteritems()])
            data_total = sum([v for _, v in self.data.iteritems()])

        formula_values = {}
        for key, formula in self.formulas.iteritems():
            formula_values[key] = self.apply_formula(formula)

        # Get formatting settings
        longest_name = max(self.all_names(), key=len)
        space = len(longest_name) + 2

        # Print Table
        print "==============  Yield Table =============="
        print " Variable(s): ", self.variable
        print " Region: ", self.region
        print "=========================================="
        if len(self.mc):
            print "------------------------------------------"
            print "Backgrounds:"
            for name, yield_value in zip(self.mc.keys(), bkg_strings):
                print "%*s : %s"%(space, name, yield_value)
        if len(self.data):
            print "------------------------------------------"
            print "%*s : %s"%(space, "MC", mc_total)
            print "%*s : %.*f"%(space, "Data", self.precision, data_total.value)

        if len(self.signals):
            print "------------------------------------------"
            print "Signal:"
            for name, yield_value in zip(self.signals.keys(), bkg_strings):
                print "%*s : %s"%(space, name, yield_value)
        if self.data_mc_ratio or len(self.formulas):
            print "------------------------------------------"
        if self.data_mc_ratio:
            ratio = (data_total/mc_total).value
            print "%*s : %.*f"%(space, "Data/MC", self.precision, ratio)
        if len(self.formulas):
            for name, value in formula_values.iteritems():
                print "%*s : %.*f"%(space, name, self.precision, value)
        print "=========================================="

    def apply_formula(self, formula, no_uncertainty=True):
        # Get samples from formula
        samples = re.sub("(\+|\-|\)|\(|\*|\/)"," ", formula)
        samples = [s.strip() for s in samples.split()]
        assert all(s.replace("_","").isalpha() for s in samples), (
            "ERROR :: Unacceptable formula format", formula)


        # Replace names in formula with values
        all_yields = self.all_yields()
        assert len(all_yields) == len(self.mc) + len(self.data) + len(self.signals), (
            "ERROR (YieldTable) :: Overlapping key values")
        samples.sort(key=len, reverse=True)

        for sample_name in samples:
            assert sample_name in all_yields, (
                "ERROR :: Formula sample not stored:", sample_name)
            if no_uncertainty:
                replace_str = "all_yields['%s'].value"%sample_name
            else:
                replace_str = "all_yields['%s']"%sample_name
            formula = re.sub(r"\b%s\b"%sample_name, replace_str, formula)

        # Evaluate formula
        return eval(formula)

    def reset(self):
        self.mc.clear()
        self.data.clear()
        self.signals.clear()
        self.region = ""
        self.variable = ""


