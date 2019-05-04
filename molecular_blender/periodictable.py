# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: periodictable.py
#  Copyright (C) 2014 Shane Parker, Joshua Szekely
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

"""Collection of information for Periodic Table of Elements"""

class Element():
    """Element of the Periodic Table"""

    def __init__(self, vdw, cov, m, n, s):
        """Construct Element"""
        self.vdw = vdw
        self.covalent = cov
        self.mass = m
        self.name = n
        self.symbol = s

    def __str__(self):
        """Return string of symbol"""
        return self.symbol

elements = {
    # key: vdw, cov, mass, full name, symbol
    'q':  Element(1.00, 0.50,           0.0, "dummy", "q"),
    'h':  Element(1.20, 0.31,       1.00784, "Hydrogen", "h"),
    'he': Element(1.40, 0.28,     4.0026022, "Helium", "he"),
    'li': Element(1.82, 1.28,         6.938, "Lithium", "li"),
    'be': Element(1.00, 0.96,    9.01218315, "Beryllium", "be"),
    'b':  Element(1.00, 0.84,        10.806, "Boron", "b"),
    'c':  Element(1.70, 0.76,       12.0096, "Carbon", "c"),
    'n':  Element(1.55, 0.71,      14.00643, "Nitrogen", "n"),
    'o':  Element(1.52, 0.66,      15.99903, "Oxygen", "o"),
    'f':  Element(1.47, 0.57, 18.9984031636, "Fluorine", "f"),
    'ne': Element(1.54, 0.58,      20.17976, "Neon", "ne"),
    'na': Element(2.27, 1.66,  22.989769282, "Sodium", "na"),
    'mg': Element(1.73, 1.41,        24.304, "Magnesium", "mg"),
    'al': Element(1.84, 1.21,   26.98153857, "Aluminum", "al"),
    'si': Element(2.10, 1.11,        28.084, "Silicon", "si"),
    'p':  Element(1.80, 1.07, 30.9737619985, "Phosphorus", "p"),
    's':  Element(1.80, 1.05,        32.059, "Sulfur", "s"),
    'cl': Element(1.75, 1.02,        35.446, "Chlorine", "cl"),
    'ar': Element(1.88, 1.06,       39.9481, "Argon", "ar"),
    'k':  Element(2.75, 2.03,      39.09831, "Potassium", "k"),
    'ca': Element(2.31, 1.76,       40.0784, "Calcium", "ca"),
    'sc': Element(1.00, 1.70,    44.9559085, "Scandium", "sc"),
    'ti': Element(2.25, 1.36,       47.8671, "Titanium", "ti"),
    'v':  Element(1.00, 1.53,      50.94151, "Vanadium", "v"),
    'cr': Element(1.00, 1.39,      51.99616, "Chromium", "cr"),
    'mn': Element(1.00, 1.39,    54.9380443, "Manganese", "mn"),
    'fe': Element(1.00, 1.32,       55.8452, "Iron", "fe"),
    'co': Element(1.00, 1.26,    58.9331944, "Cobalt", "co"),
    'ni': Element(1.63, 1.24,      58.69344, "Nickel", "ni"),
    'cu': Element(1.45, 1.32,       63.5463, "Copper", "cu"),
    'zn': Element(1.42, 1.22,        65.382, "Zinc", "zn"),
    'ga': Element(1.87, 1.22,       69.7231, "Gallium", "ga"),
    'ge': Element(1.00, 1.20,       72.6308, "Germanium", "ge"),
    'as': Element(1.85, 1.19,    74.9215956, "Arsenic", "as"),
    'se': Element(1.90, 1.20,       78.9718, "Selenium", "se"),
    'br': Element(1.85, 1.20,        79.901, "Bromine", "br"),
    'kr': Element(2.02, 1.16,       83.7982, "Krypton", "kr"),
    'rb': Element(1.00, 2.20,      85.46783, "Rubidium", "rb"),
    'sr': Element(1.00, 1.95,        87.621, "Strontium", "sr"),
    'y':  Element(1.00, 1.90,     88.905842, "Yttrium", "y"),
    'zr': Element(1.00, 1.75,       91.2242, "Zirconium", "zr"),
    'nb': Element(1.00, 1.64,     92.906372, "Niobium", "nb"),
    'mo': Element(1.00, 1.54,        95.951, "Molybdenum", "mo"),
    'tc': Element(1.00, 1.47,            98, "Technetium", "tc"),
    'ru': Element(1.00, 1.46,       101.072, "Ruthenium", "ru"),
    'rh': Element(1.00, 1.42,    102.905502, "Rhodium", "rh"),
    'pd': Element(1.63, 1.39,       106.421, "Palladium", "pd"),
    'ag': Element(1.72, 1.45,     107.86822, "Silver", "ag"),
    'cd': Element(1.58, 1.44,      112.4144, "Cadmium", "cd"),
    'in': Element(1.93, 1.42,      114.8181, "Indium", "in"),
    'sn': Element(2.17, 1.39,      118.7107, "Tin", "sn"),
    'sb': Element(1.00, 1.39,      121.7601, "Antimony", "sb"),
    'te': Element(2.06, 1.38,       127.603, "Tellurium", "te"),
    'i':  Element(1.98, 1.39,    126.904473, "Iodine", "i"),
    'xe': Element(2.16, 1.40,      131.2936, "Xenon", "xe"),
    'cs': Element(1.00, 2.44, 132.905451966, "Cesium", "cs"),
    'ba': Element(2.68, 2.15,      137.3277, "Barium", "ba"),
    'la': Element(1.00, 2.07,    138.905477, "Lanthanum", "la"),
    'ce': Element(1.00, 2.04,      140.1161, "Cerium", "ce"),
    'pr': Element(1.00, 2.03,    140.907662, "Praseodymium", "pr"),
    'nd': Element(1.00, 2.01,      144.2423, "Neodymium", "nd"),
    'pm': Element(1.00, 1.99,           145, "Promethium", "pm"),
    'sm': Element(1.00, 1.98,       150.362, "Samarium", "sm"),
    'eu': Element(1.00, 1.98,      151.9641, "Europium", "eu"),
    'gd': Element(1.00, 1.96,       157.253, "Gadolinium", "gd"),
    'tb': Element(1.00, 1.94,    158.925352, "Terbium", "tb"),
    'dy': Element(1.00, 1.92,      162.5001, "Dysprosium", "dy"),
    'ho': Element(1.00, 1.92,    164.930332, "Holmium", "ho"),
    'er': Element(1.00, 1.89,      167.2593, "Erbium", "er"),
    'tm': Element(1.00, 1.90,    168.934222, "Thulium", "tm"),
    'yb': Element(1.00, 1.87,      173.0545, "Ytterbium", "yb"),
    'lu': Element(1.00, 1.87,     174.96681, "Lutetium", "lu"),
    'hf': Element(1.00, 1.75,       178.492, "Hafnium", "hf"),
    'ta': Element(1.00, 1.70,    180.947882, "Tantalum", "ta"),
    'w':  Element(1.00, 1.62,       183.841, "Tungsten", "w"),
    're': Element(1.00, 1.51,      186.2071, "Rhenium", "re"),
    'os': Element(1.00, 1.44,       190.233, "Osmium", "os"),
    'ir': Element(1.00, 1.41,      192.2173, "Iridium", "ir"),
    'pt': Element(1.75, 1.36,      195.0849, "Platinum", "pt"),
    'au': Element(1.66, 1.36,   196.9665695, "Gold", "au"),
    'hg': Element(1.55, 1.32,      200.5923, "Mercury", "hg"),
    'tl': Element(1.96, 1.45,       204.382, "Thallium", "tl"),
    'pb': Element(2.02, 1.46,        207.21, "Lead", "pb"),
    'bi': Element(1.00, 1.48,    208.980401, "Bismuth", "bi"),
    'po': Element(1.00, 1.40,           209, "Polonium", "po"),
    'at': Element(1.00, 1.50,           210, "Astatine", "at"),
    'rn': Element(2.20, 1.50,           222, "Radon", "rn"),
    'fr': Element(1.00, 2.60,           223, "Francium", "fr"),
    'ra': Element(1.00, 2.21,           226, "Radium", "ra"),
    'ac': Element(1.00, 2.15,           227, "Actinium", "ac"),
    'th': Element(1.00, 2.06,     232.03774, "Thorium", "th"),
    'pa': Element(1.00, 2.00,    231.035882, "Protactinium", "pa"),
    'u':  Element(1.00, 1.96,    238.028913, "Uranium", "u"),
    'np': Element(1.00, 1.90,           237, "Neptunium", "np"),
    'pu': Element(1.00, 1.87,           244, "Plutonium", "pu"),
    'am': Element(1.00, 1.80, 241.056829319, "Americium", "am"),
    'cm': Element(1.00, 1.69, 243.061389322, "Curium", "cm"),
    'bk': Element(1.00, 1.00, 247.070307359, "Berkelium", "bk"),
    'cf': Element(1.00, 1.00, 249.074853923, "Californium", "cf"),
    'es': Element(1.00, 1.00,  252.08298054, "Einsteinium", "es"),
    'fm': Element(1.00, 1.00, 257.095106169, "Fermium", "fm"),
    'md': Element(1.00, 1.00, 258.098431550, "Mendelevium", "md"),
    'no': Element(1.00, 1.00,   259.1010311, "Nobelium", "no"),
    'lr': Element(1.00, 1.00,   262.1096122, "Lawrencium", "lr"),
    'rf': Element(1.00, 1.00,   267.1217962, "Rutherfordium", "rf"),
    'db': Element(1.00, 1.00,   268.1256757, "Dubnium", "db"),
    'sg': Element(1.00, 1.00,   271.1339363, "Seaborgium", "sg"),
    'bh': Element(1.00, 1.00,   272.1382658, "Bohrium", "bh"),
    'hs': Element(1.00, 1.00,   270.1342927, "Hassium", "hs"),
    'mt': Element(1.00, 1.00,   276.1515959, "Meitnerium", "mt")
}


symbols = ["q",
           "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na", "mg", "al",
           "si", "p", "s", "cl", "ar", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe",
           "co", "ni", "cu", "zn", "ga", "ge", "as", "se", "br", "kr", "rb", "sr",
           "y", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn",
           "sb", "te", "i", "xe", "cs", "ba", "la", "ce", "pr", "nd", "pm", "sm",
           "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu", "hf", "ta", "w",
           "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi", "po", "at", "rn",
           "fr", "ra", "ac", "th", "pa", "u", "np", "pu", "am", "cm", "bk", "cf",
           "es", "fm", "md", "no", "lr", "rf", "db", "sg", "bh", "hs", "mt", "ds",
           "rg", "cn", "uut", "fl", "uup", "lv", "uus", "uuo"]
