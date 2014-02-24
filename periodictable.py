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

class element():
    vdw = 2.0
    covalent = 1.0
    color = (1.0, 1.0, 1.0)
    mass = 1.0
    name = ""
    symbol = ""

    def __init__(self, vdw, cov, col, m, n, s):
        self.vdw = vdw
        self.covalent = cov
        self.color = col
        self.mass = m
        self.name= n
        self.symbol = s

elements = {
  'h'  : element(1.20, 0.31, (1.000,1.000,1.000),   1.0, "Hydrogen", "h"),
  'he' : element(1.40, 0.28, (0.851,1.000,1.000),   4.0, "Helium", "he"),
  'li' : element(1.82, 1.28, (0.800,0.502,1.000),   6.9, "Lithium", "li"),
  'be' : element(1.00, 0.96, (0.761,1.000,0.000),   9.0, "Beryllium", "be"),
  'b'  : element(1.00, 0.84, (1.000,0.710,0.710),  10.8, "Boron", "b"),
  'c'  : element(1.70, 0.76, (0.565,0.565,0.565),  12.0, "Carbon", "c"),
  'n'  : element(1.55, 0.71, (0.188,0.314,0.973),  14.0, "Nitrogen", "n"),
  'o'  : element(1.52, 0.66, (1.000,0.051,0.051),  16.0, "Oxygen", "o"),
  'f'  : element(1.47, 0.57, (0.565,0.878,0.314),  19.0, "Fluorine", "f"),
  'ne' : element(1.54, 0.58, (0.702,0.890,0.961),  20.2, "Neon", "ne"),
  'na' : element(2.27, 1.66, (0.671,0.361,0.949),  23.0, "Sodium", "na"),
  'mg' : element(1.73, 1.41, (0.541,1.000,0.000),  24.3, "Magnesium", "mg"),
  'al' : element(1.84, 1.21, (0.749,0.651,0.651),  27.0, "Aluminum", "al"),
  'si' : element(2.10, 1.11, (0.941,0.784,0.627),  28.1, "Silicon", "si"),
  'p'  : element(1.80, 1.07, (1.000,0.502,0.000),  31.0, "Phosphorus", "p"),
  's'  : element(1.80, 1.05, (1.000,1.000,0.188),  32.1, "Sulfur", "s"),
  'cl' : element(1.75, 1.02, (0.122,0.941,0.122),  35.5, "Chlorine", "cl"),
  'ar' : element(1.88, 1.06, (0.502,0.820,0.890),  39.9, "Argon", "ar"),
  'k'  : element(2.75, 2.03, (0.561,0.251,0.831),  39.1, "Potassium", "k"),
  'ca' : element(2.31, 1.76, (0.239,1.000,0.000),  40.1, "Calcium", "ca"),
  'sc' : element(1.00, 1.70, (0.902,0.902,0.902),  45.0, "Scandium", "sc"),
  'ti' : element(1.00, 1.60, (0.749,0.761,0.780),  47.9, "Titanium", "ti"),
  'v'  : element(1.00, 1.53, (0.651,0.651,0.671),  50.9, "Vanadium", "v"),
  'cr' : element(1.00, 1.39, (0.541,0.600,0.780),  52.0, "Chromium", "cr"),
  'mn' : element(1.00, 1.39, (0.612,0.478,0.780),  54.9, "Manganese", "mn"),
  'fe' : element(1.00, 1.32, (0.878,0.400,0.200),  55.8, "Iron", "fe"),
  'co' : element(1.00, 1.26, (0.941,0.565,0.627),  58.9, "Cobalt", "co"),
  'ni' : element(1.63, 1.24, (0.314,0.816,0.314),  58.7, "Nickel", "ni"),
  'cu' : element(1.45, 1.32, (0.784,0.502,0.200),  63.5, "Copper", "cu"),
  'zn' : element(1.42, 1.22, (0.490,0.502,0.690),  65.4, "Zinc", "zn"),
  'ga' : element(1.87, 1.22, (0.761,0.561,0.561),  69.7, "Gallium", "ga"),
  'ge' : element(1.00, 1.20, (0.400,0.561,0.561),  72.6, "Germanium", "ge"),
  'as' : element(1.85, 1.19, (0.741,0.502,0.890),  74.9, "Arsenic", "as"),
  'se' : element(1.90, 1.20, (1.000,0.631,0.000),  79.0, "Selenium", "se"),
  'br' : element(1.85, 1.20, (0.651,0.161,0.161),  79.9, "Bromine", "br"),
  'kr' : element(2.02, 1.16, (0.361,0.722,0.820),  83.8, "Krypton", "kr"),
  'rb' : element(1.00, 2.20, (0.439,0.180,0.690),  85.5, "Rubidium", "rb"),
  'sr' : element(1.00, 1.95, (0.000,1.000,0.000),  87.6, "Strontium", "sr"),
  'y'  : element(1.00, 1.90, (0.580,1.000,1.000),  88.9, "Yttrium", "y"),
  'zr' : element(1.00, 1.75, (0.580,0.878,0.878),  91.2, "Zirconium", "zr"),
  'nb' : element(1.00, 1.64, (0.451,0.761,0.788),  92.9, "Niobium", "nb"),
  'mo' : element(1.00, 1.54, (0.329,0.710,0.710),  95.9, "Molybdenum", "mo"),
  'tc' : element(1.00, 1.47, (0.231,0.620,0.620),  98.0, "Technetium", "tc"),
  'ru' : element(1.00, 1.46, (0.141,0.561,0.561), 101.1, "Ruthenium", "ru"),
  'rh' : element(1.00, 1.42, (0.039,0.490,0.549), 102.9, "Rhodium", "rh"),
  'pd' : element(1.63, 1.39, (0.000,0.412,0.522), 106.4, "Palladium", "pd"),
  'ag' : element(1.72, 1.45, (0.753,0.753,0.753), 107.9, "Silver", "ag"),
  'cd' : element(1.58, 1.44, (1.000,0.851,0.561), 112.4, "Cadmium", "cd"),
  'in' : element(1.93, 1.42, (0.651,0.459,0.451), 114.8, "Indium", "in"),
  'sn' : element(2.17, 1.39, (0.400,0.502,0.502), 118.7, "Tin", "sn"),
  'sb' : element(1.00, 1.39, (0.620,0.388,0.710), 121.8, "Antimony", "sb"),
  'te' : element(2.06, 1.38, (0.831,0.478,0.000), 127.6, "Tellurium", "te"),
  'i'  : element(1.98, 1.39, (0.580,0.000,0.580), 126.9, "Iodine", "i"),
  'xe' : element(2.16, 1.40, (0.259,0.620,0.690), 131.3, "Xenon", "xe"),
  'cs' : element(1.00, 2.44, (0.341,0.090,0.561), 132.9, "Cesium", "cs"),
  'ba' : element(2.68, 2.15, (0.000,0.788,0.000), 137.3, "Barium", "ba"),
  'la' : element(1.00, 2.07, (0.439,0.831,1.000), 138.9, "Lanthanum", "la"),
  'ce' : element(1.00, 2.04, (1.000,1.000,0.780), 140.1, "Cerium", "ce"),
  'pr' : element(1.00, 2.03, (0.851,1.000,0.780), 140.9, "Praseodymium", "pr"),
  'nd' : element(1.00, 2.01, (0.780,1.000,0.780), 144.2, "Neodymium", "nd"),
  'pm' : element(1.00, 1.99, (0.639,1.000,0.780), 145.0, "Promethium", "pm"),
  'sm' : element(1.00, 1.98, (0.561,1.000,0.780), 150.4, "Samarium", "sm"),
  'eu' : element(1.00, 1.98, (0.380,1.000,0.780), 152.0, "Europium", "eu"),
  'gd' : element(1.00, 1.96, (0.271,1.000,0.780), 157.2, "Gadolinium", "gd"),
  'tb' : element(1.00, 1.94, (0.188,1.000,0.780), 158.9, "Terbium", "tb"),
  'dy' : element(1.00, 1.92, (0.122,1.000,0.780), 162.5, "Dysprosium", "dy"),
  'ho' : element(1.00, 1.92, (0.000,1.000,0.612), 164.9, "Holmium", "ho"),
  'er' : element(1.00, 1.89, (0.000,0.902,0.459), 167.3, "Erbium", "er"),
  'tm' : element(1.00, 1.90, (0.000,0.831,0.322), 168.9, "Thulium", "tm"),
  'yb' : element(1.00, 1.87, (0.000,0.749,0.220), 173.0, "Ytterbium", "yb"),
  'lu' : element(1.00, 1.87, (0.000,0.671,0.141), 175.0, "Lutetium", "lu"),
  'hf' : element(1.00, 1.75, (0.302,0.761,1.000), 178.5, "Hafnium", "hf"),
  'ta' : element(1.00, 1.70, (0.302,0.651,1.000), 180.9, "Tantalum", "ta"),
  'w'  : element(1.00, 1.62, (0.129,0.580,0.839), 183.8, "Tungsten", "w"),
  're' : element(1.00, 1.51, (0.149,0.490,0.671), 186.2, "Rhenium", "re"),
  'os' : element(1.00, 1.44, (0.149,0.400,0.588), 190.2, "Osmium", "os"),
  'ir' : element(1.00, 1.41, (0.090,0.329,0.529), 192.2, "Iridium", "ir"),
  'pt' : element(1.75, 1.36, (0.816,0.816,0.878), 195.1, "Platinum", "pt"),
  'au' : element(1.66, 1.36, (1.000,0.820,0.137), 197.0, "Gold", "au"),
  'hg' : element(1.55, 1.32, (0.722,0.722,0.816), 200.6, "Mercury", "hg"),
  'tl' : element(1.96, 1.45, (0.651,0.329,0.302), 204.4, "Thallium", "tl"),
  'pb' : element(2.02, 1.46, (0.341,0.349,0.380), 207.2, "Lead", "pb"),
  'bi' : element(1.00, 1.48, (0.620,0.310,0.710), 209.0, "Bismuth", "bi"),
  'po' : element(1.00, 1.40, (0.671,0.361,0.000), 209.0, "Polonium", "po"),
  'at' : element(1.00, 1.50, (0.459,0.310,0.271), 210.0, "Astatine", "at"),
  'rn' : element(2.20, 1.50, (0.259,0.510,0.588), 222.0, "Radon", "rn"),
  'fr' : element(1.00, 2.60, (0.259,0.000,0.400), 223.0, "Francium", "fr"),
  'ra' : element(1.00, 2.21, (0.000,0.490,0.000), 226.0, "Radium", "ra"),
  'ac' : element(1.00, 2.15, (0.439,0.671,0.980), 227.0, "Actinium", "ac"),
  'th' : element(1.00, 2.06, (0.000,0.729,1.000), 232.0, "Thorium", "th"),
  'pa' : element(1.00, 2.00, (0.000,0.631,1.000), 231.0, "Protactinium", "pa"),
  'u'  : element(1.00, 1.96, (0.000,0.561,1.000), 238.0, "Uranium", "u"),
  'np' : element(1.00, 1.90, (0.000,0.502,1.000), 237.0, "Neptunium", "np"),
  'pu' : element(1.00, 1.87, (0.000,0.420,1.000), 244.0, "Plutonium", "pu"),
  'am' : element(1.00, 1.80, (0.329,0.361,0.949), 243.0, "Americium", "am"),
  'cm' : element(1.00, 1.69, (0.471,0.361,0.890), 247.0, "Curium", "cm"),
  'bk' : element(1.00, 1.00, (0.541,0.310,0.890), 247.0, "Berkelium", "bk"),
  'cf' : element(1.00, 1.00, (0.631,0.212,0.831), 251.0, "Californium", "cf"),
  'es' : element(1.00, 1.00, (0.702,0.122,0.831), 252.0, "Einsteinium", "es"),
  'fm' : element(1.00, 1.00, (0.702,0.122,0.729), 257.0, "Fermium", "fm"),
  'md' : element(1.00, 1.00, (0.702,0.051,0.651), 258.0, "Mendelevium", "md"),
  'no' : element(1.00, 1.00, (0.741,0.051,0.529), 259.0, "Nobelium", "no"),
  'lr' : element(1.00, 1.00, (0.780,0.000,0.400), 262.0, "Lawrencium", "lr"),
  'rf' : element(1.00, 1.00, (0.800,0.000,0.349), 261.0, "Rutherfordium", "rf"),
  'db' : element(1.00, 1.00, (0.820,0.000,0.310), 262.0, "Dubnium", "db"),
  'sg' : element(1.00, 1.00, (0.851,0.000,0.271), 266.0, "Seaborgium", "sg"),
  'bh' : element(1.00, 1.00, (0.878,0.000,0.220), 264.0, "Bohrium", "bh"),
  'hs' : element(1.00, 1.00, (0.902,0.000,0.180), 277.0, "Hassium", "hs"),
  'mt' : element(1.00, 1.00, (0.922,0.000,0.149), 268.0, "Meitnerium", "mt")
}

symbols = [ "empty",
            "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na", "mg", "al",
            "si", "p", "s", "cl", "ar", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe",
            "co", "ni", "cu", "zn", "ga", "ge", "as", "se", "br", "kr", "rb", "sr",
            "y", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn",
            "sb", "te", "i", "xe", "cs", "ba", "la", "ce", "pr", "nd", "pm", "sm",
            "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu", "hf", "ta", "w",
            "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi", "po", "at", "rn",
            "fr", "ra", "ac", "th", "pa", "u", "np", "pu", "am", "cm", "bk", "cf",
            "es", "fm", "md", "no", "lr", "rf", "db", "sg", "bh", "hs", "mt", "ds",
            "rg", "cn", "uut", "fl", "uup", "lv", "uus", "uuo"
          ]
