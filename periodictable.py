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
    'h':  Element(1.20, 0.31,   1.0, "Hydrogen", "h"),
    'he': Element(1.40, 0.28,   4.0, "Helium", "he"),
    'li': Element(1.82, 1.28,   6.9, "Lithium", "li"),
    'be': Element(1.00, 0.96,   9.0, "Beryllium", "be"),
    'b':  Element(1.00, 0.84,  10.8, "Boron", "b"),
    'c':  Element(1.70, 0.76,  12.0, "Carbon", "c"),
    'n':  Element(1.55, 0.71,  14.0, "Nitrogen", "n"),
    'o':  Element(1.52, 0.66,  16.0, "Oxygen", "o"),
    'f':  Element(1.47, 0.57,  19.0, "Fluorine", "f"),
    'ne': Element(1.54, 0.58,  20.2, "Neon", "ne"),
    'na': Element(2.27, 1.66,  23.0, "Sodium", "na"),
    'mg': Element(1.73, 1.41,  24.3, "Magnesium", "mg"),
    'al': Element(1.84, 1.21,  27.0, "Aluminum", "al"),
    'si': Element(2.10, 1.11,  28.1, "Silicon", "si"),
    'p':  Element(1.80, 1.07,  31.0, "Phosphorus", "p"),
    's':  Element(1.80, 1.05,  32.1, "Sulfur", "s"),
    'cl': Element(1.75, 1.02,  35.5, "Chlorine", "cl"),
    'ar': Element(1.88, 1.06,  39.9, "Argon", "ar"),
    'k':  Element(2.75, 2.03,  39.1, "Potassium", "k"),
    'ca': Element(2.31, 1.76,  40.1, "Calcium", "ca"),
    'sc': Element(1.00, 1.70,  45.0, "Scandium", "sc"),
    'ti': Element(2.25, 1.36,  47.9, "Titanium", "ti"),
    'v':  Element(1.00, 1.53,  50.9, "Vanadium", "v"),
    'cr': Element(1.00, 1.39,  52.0, "Chromium", "cr"),
    'mn': Element(1.00, 1.39,  54.9, "Manganese", "mn"),
    'fe': Element(1.00, 1.32,  55.8, "Iron", "fe"),
    'co': Element(1.00, 1.26,  58.9, "Cobalt", "co"),
    'ni': Element(1.63, 1.24,  58.7, "Nickel", "ni"),
    'cu': Element(1.45, 1.32,  63.5, "Copper", "cu"),
    'zn': Element(1.42, 1.22,  65.4, "Zinc", "zn"),
    'ga': Element(1.87, 1.22,  69.7, "Gallium", "ga"),
    'ge': Element(1.00, 1.20,  72.6, "Germanium", "ge"),
    'as': Element(1.85, 1.19,  74.9, "Arsenic", "as"),
    'se': Element(1.90, 1.20,  79.0, "Selenium", "se"),
    'br': Element(1.85, 1.20,  79.9, "Bromine", "br"),
    'kr': Element(2.02, 1.16,  83.8, "Krypton", "kr"),
    'rb': Element(1.00, 2.20,  85.5, "Rubidium", "rb"),
    'sr': Element(1.00, 1.95,  87.6, "Strontium", "sr"),
    'y':  Element(1.00, 1.90,  88.9, "Yttrium", "y"),
    'zr': Element(1.00, 1.75,  91.2, "Zirconium", "zr"),
    'nb': Element(1.00, 1.64,  92.9, "Niobium", "nb"),
    'mo': Element(1.00, 1.54,  95.9, "Molybdenum", "mo"),
    'tc': Element(1.00, 1.47,  98.0, "Technetium", "tc"),
    'ru': Element(1.00, 1.46, 101.1, "Ruthenium", "ru"),
    'rh': Element(1.00, 1.42, 102.9, "Rhodium", "rh"),
    'pd': Element(1.63, 1.39, 106.4, "Palladium", "pd"),
    'ag': Element(1.72, 1.45, 107.9, "Silver", "ag"),
    'cd': Element(1.58, 1.44, 112.4, "Cadmium", "cd"),
    'in': Element(1.93, 1.42, 114.8, "Indium", "in"),
    'sn': Element(2.17, 1.39, 118.7, "Tin", "sn"),
    'sb': Element(1.00, 1.39, 121.8, "Antimony", "sb"),
    'te': Element(2.06, 1.38, 127.6, "Tellurium", "te"),
    'i':  Element(1.98, 1.39, 126.9, "Iodine", "i"),
    'xe': Element(2.16, 1.40, 131.3, "Xenon", "xe"),
    'cs': Element(1.00, 2.44, 132.9, "Cesium", "cs"),
    'ba': Element(2.68, 2.15, 137.3, "Barium", "ba"),
    'la': Element(1.00, 2.07, 138.9, "Lanthanum", "la"),
    'ce': Element(1.00, 2.04, 140.1, "Cerium", "ce"),
    'pr': Element(1.00, 2.03, 140.9, "Praseodymium", "pr"),
    'nd': Element(1.00, 2.01, 144.2, "Neodymium", "nd"),
    'pm': Element(1.00, 1.99, 145.0, "Promethium", "pm"),
    'sm': Element(1.00, 1.98, 150.4, "Samarium", "sm"),
    'eu': Element(1.00, 1.98, 152.0, "Europium", "eu"),
    'gd': Element(1.00, 1.96, 157.2, "Gadolinium", "gd"),
    'tb': Element(1.00, 1.94, 158.9, "Terbium", "tb"),
    'dy': Element(1.00, 1.92, 162.5, "Dysprosium", "dy"),
    'ho': Element(1.00, 1.92, 164.9, "Holmium", "ho"),
    'er': Element(1.00, 1.89, 167.3, "Erbium", "er"),
    'tm': Element(1.00, 1.90, 168.9, "Thulium", "tm"),
    'yb': Element(1.00, 1.87, 173.0, "Ytterbium", "yb"),
    'lu': Element(1.00, 1.87, 175.0, "Lutetium", "lu"),
    'hf': Element(1.00, 1.75, 178.5, "Hafnium", "hf"),
    'ta': Element(1.00, 1.70, 180.9, "Tantalum", "ta"),
    'w':  Element(1.00, 1.62, 183.8, "Tungsten", "w"),
    're': Element(1.00, 1.51, 186.2, "Rhenium", "re"),
    'os': Element(1.00, 1.44, 190.2, "Osmium", "os"),
    'ir': Element(1.00, 1.41, 192.2, "Iridium", "ir"),
    'pt': Element(1.75, 1.36, 195.1, "Platinum", "pt"),
    'au': Element(1.66, 1.36, 197.0, "Gold", "au"),
    'hg': Element(1.55, 1.32, 200.6, "Mercury", "hg"),
    'tl': Element(1.96, 1.45, 204.4, "Thallium", "tl"),
    'pb': Element(2.02, 1.46, 207.2, "Lead", "pb"),
    'bi': Element(1.00, 1.48, 209.0, "Bismuth", "bi"),
    'po': Element(1.00, 1.40, 209.0, "Polonium", "po"),
    'at': Element(1.00, 1.50, 210.0, "Astatine", "at"),
    'rn': Element(2.20, 1.50, 222.0, "Radon", "rn"),
    'fr': Element(1.00, 2.60, 223.0, "Francium", "fr"),
    'ra': Element(1.00, 2.21, 226.0, "Radium", "ra"),
    'ac': Element(1.00, 2.15, 227.0, "Actinium", "ac"),
    'th': Element(1.00, 2.06, 232.0, "Thorium", "th"),
    'pa': Element(1.00, 2.00, 231.0, "Protactinium", "pa"),
    'u':  Element(1.00, 1.96, 38.0, "Uranium", "u"),
    'np': Element(1.00, 1.90, 237.0, "Neptunium", "np"),
    'pu': Element(1.00, 1.87, 244.0, "Plutonium", "pu"),
    'am': Element(1.00, 1.80, 243.0, "Americium", "am"),
    'cm': Element(1.00, 1.69, 247.0, "Curium", "cm"),
    'bk': Element(1.00, 1.00, 247.0, "Berkelium", "bk"),
    'cf': Element(1.00, 1.00, 251.0, "Californium", "cf"),
    'es': Element(1.00, 1.00, 252.0, "Einsteinium", "es"),
    'fm': Element(1.00, 1.00, 257.0, "Fermium", "fm"),
    'md': Element(1.00, 1.00, 258.0, "Mendelevium", "md"),
    'no': Element(1.00, 1.00, 259.0, "Nobelium", "no"),
    'lr': Element(1.00, 1.00, 262.0, "Lawrencium", "lr"),
    'rf': Element(1.00, 1.00, 261.0, "Rutherfordium", "rf"),
    'db': Element(1.00, 1.00, 262.0, "Dubnium", "db"),
    'sg': Element(1.00, 1.00, 266.0, "Seaborgium", "sg"),
    'bh': Element(1.00, 1.00, 264.0, "Bohrium", "bh"),
    'hs': Element(1.00, 1.00, 277.0, "Hassium", "hs"),
    'mt': Element(1.00, 1.00, 268.0, "Meitnerium", "mt")
}


symbols = ["empty",
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
