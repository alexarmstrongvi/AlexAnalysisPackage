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
from number import Numbers
from math import sqrt
from collections import OrderedDict

class UncFloat :
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
        return (self.value == other.value
            and self.uncertainty == other.uncertainty)
    def __ne__(self, other) :
        return not self.value == other.value
    def __gt__(self, other) :
        return self.value > other.value
    def __ge__(self, other) :
        return self.value >= other.value

    def __str__(self):
        return "%s +/- %s"%(self.value, self.uncertainty)

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
    def __init__(self):
        self.mc = OrderedDict()
        self.data = {}
        self.signals = {}
        self.data_mc_ratio = True
        self.formulas = {}
        self.region = ""
        self.variable = ""

    def Print(no_uncertainty = False):
        if no_uncertainty:
            bkg_strings = [str(v.value) for _, v in self.mc.iteritems()]
            sig_strings = [str(v.value) for _, v in self.signals.iteritems()]
            data_string = [str(v.value) for _, v in self.data.iteritems()]
            mc_total = sum([v.value for _, v in self.mc.iteritems()])
            data_total = sum([v.value for _, v in self.data.iteritems()])
        else:
            bkd_strings = [str(v) for _, v in self.mc.iteritems()]
            sig_strings = [str(v) for _, v in self.signals.iteritems()]
            data_string = [str(v) for _, v in self.data.iteritems()]
            mc_total = sum([v for _, v in self.mc.iteritems()])
            data_total = sum([v for _, v in self.data.iteritems()])

        formula_values = {}
        for key, formula in numerators.iteritems():
            formula_values[key] = apply_formula(formula)

        # Get formatting settings
        all_names = []
        all_names += [k for k in self.mc]
        all_names += [k for k in self.signals]
        all_names += [k for k in self.data]
        all_names += [k for k in self.formulas]
        longest_name = max(all_names, key=len)
        space = len(longest_name) + 2

        # Print Table
        print "===== Yields (Region : %s, Quantity : %s) ====="%(self.region, self.variable)
        if len(self.mc):
            print "-"*40
            print "Backgrounds:"
            for name, yield_value in zip(self.mc.keys(), bkg_strings):
                print "%*s : %s"%(space, name, yield_value)
        if len(self.data):
            print "-"*40
            print "Data:"
            for name, yield_value in zip(self.data.keys(), data_string):
                print "%*s : %s"%(space, name, yield_value)

        if len(self.signals):
            print "-"*40
            print "Signal:"
            for name, yield_value in zip(self.signals.keys(), bkg_strings):
                print "%*s : %s"%(space, name, yield_value)
        if self.data_mc_ratio or len(self.formulas):
            print "-"*40
        if self.data_mc_ratio:
            print "%*s : %s"%(space, "Data/MC", data_total/mc_total)
        if len(self.formulas):
            for name, value in formula_values:
                print "%*s : %s"%(space, name, value)

    def apply_formula(self, formula, no_uncertainty=True):
        # Get samples from formula
        samples = formula.replace("+"," ").replace("-"," ")
        samples = samples.replace("*"," ").replace("/"," ")
        samples = samples.replace(")"," ").replace("("," ")
        samples = [s.strip() for s in samples.split()]
        assert all(s.replace("_","").isalpha() for s in samples), (
            "ERROR :: Unacceptable formula format", formula)


        # Replace names in formula with values
        all_values = dict(self.mc, **self.data, **self.signals)
        assert len(all_values) == len(self.mc) + len(self.data) + len(self.signals), (
            "ERROR (YieldTable) :: Overlapping key values")

        for sample_name in samples.sort(key=len, reverse=True)
            assert sample_name in all_values, (
                "ERROR :: Formula sample not stored:", sample_name)
            if no_uncertainty
                formula.replace(sample_name, "all_values['%s'].value"%sample)
            else:
                formula.replace(sample_name, "all_values['%s']"%sample)

        # Evaluate formula
        return eval(formula)


    def reset(self):
        #TODO: Reset for new region/plot
        pass

    def PrettyPrint(no_uncertainty = False):


