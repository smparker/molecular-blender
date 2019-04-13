# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: util.py
#  Copyright (C) 2017 Shane Parker, Joshua Szekely
#
#  This file is part of Molecular Blender.
#
#  Molecular Blender is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3, or (at your option)
#  any later version.
#
#  Molecular Blender is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Library General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Molecular Blender; see COPYING.
#  If not, see <http://www.gnu.org/licenses/>.
#

"""Common utilities for Molecular Blender"""

import os
import time
import math

DEBUG = os.environ.get('MB_DEBUG', '') == '1'


class Timer(object):
    """Convenience class for measuring timing of chunks of code"""

    def __init__(self):
        """Start a new timer"""
        self.last_time = time.time()

    def tick(self):
        """Set a new reference time and return the time since the last tick"""
        lt, self.last_time = self.last_time, time.time()
        return self.last_time - lt

    def tick_print(self, label):
        """Calls tick and automatically prints the output with the given label"""
        out = self.tick()
        if DEBUG:
            print("  %40s: %.4f sec" % (label, out))
        return out


def stopwatch(routine, verbose=DEBUG):
    """Decorator to measure time in a function using blender timer"""
    def stopwatch_dec(func):
        """Specific decorator"""
        def wrapper(*args, **kwargs):
            """Wrapper for decorator to measure timer"""
            start = time.time()
            out = func(*args, **kwargs)
            end = time.time()
            if verbose:
                print("%.4f sec elapsed in routine %s" %
                      ((end - start), routine))
            return out
        return wrapper
    return stopwatch_dec


def unique_name(name, existing_names, starting_suffix=None):
    """If name is not in existing_names, returns name. Otherwise, returns name + "0", 1, 2, etc."""
    testname = name if starting_suffix is None else "%s%d" % (
        name, starting_suffix)
    if testname in existing_names:
        i = 0 if starting_suffix is None else starting_suffix + 1
        while True:
            testname = "%s%d" % (name, i)
            if testname in existing_names:
                i += 1
            else:
                return testname
    else:
        return testname

def round_sigfigs(num, sig_figs):
    """Round to specified number of sigfigs.

    >>> round_sigfigs(0, sig_figs=4)
    0
    >>> int(round_sigfigs(12345, sig_figs=2))
    12000
    >>> int(round_sigfigs(-12345, sig_figs=2))
    -12000
    >>> int(round_sigfigs(1, sig_figs=2))
    1
    >>> '{0:.3}'.format(round_sigfigs(3.1415, sig_figs=2))
    '3.1'
    >>> '{0:.3}'.format(round_sigfigs(-3.1415, sig_figs=2))
    '-3.1'
    >>> '{0:.5}'.format(round_sigfigs(0.00098765, sig_figs=2))
    '0.00099'
    >>> '{0:.6}'.format(round_sigfigs(0.00098765, sig_figs=3))
    '0.000988'
    """
    if num != 0:
        return round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0  # Can't take the log of 0
