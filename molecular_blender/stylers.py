# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: stylers.py
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

"""Collection of classes to construct default materials/colors for plotting"""

import bpy
import mathutils

from .nodes import arrange_nodes

class PaletteElementStyler(object):
    """Base class to color atom based on element"""

    bond_color = (0.9, 0.9, 0.9, 1.0)
    ring_color = (0.137, 0.523, 0.120, 1.0)
    charge_minus_color = (1.0, 0.0, 0.0, 1.0)
    charge_plus_color = (1.0, 0.0, 0.0, 1.0)
    iso_minus_color = (0.023153, 0.527115, 0.102242, 1.0)
    iso_plus_color = (0.40724, 0.088656, 1.0, 1.0)

    def make_principled_material(self, name, base_color=(1.0, 0.0, 0.0, 1.0),
            subsurface=0.0, subsurface_color=(1.0, 1.0, 1.0, 1.0),
            metallic=0.0, specular=0.5, specular_tint = 0.0,
            roughness=0.0, anisotropic=0.0, anisotropic_rotation=0.0,
            sheen=0.0, sheen_tint=0.0, clearcoat=0.0, clearcoat_roughness=0.0,
            IOR=1.45, transmission=0.0, transmission_roughness=0.0):

        mat = bpy.data.materials.new(name)
        mat.use_nodes = True
        mat.diffuse_color = base_color

        options = {
            'Base Color' : base_color,
            'Subsurface' : subsurface,
            'Subsurface Color' : subsurface_color,
            'Metallic' : metallic,
            'Specular' : specular,
            'Specular Tint' : specular_tint,
            'Roughness'  : roughness,
            'Anisotropic' : anisotropic,
            'Anisotropic Rotation' : anisotropic_rotation,
            'Sheen'      : sheen,
            'Sheen Tint' : sheen_tint,
            'Clearcoat'  : clearcoat,
            'Clearcoat Roughness' : clearcoat_roughness,
            'IOR' : IOR,
            'Transmission' : transmission,
            'Transmission Roughness' : transmission_roughness
        }

        nodes = mat.node_tree.nodes
        nodes.clear()
        links = mat.node_tree.links

        material_output = nodes.new('ShaderNodeOutputMaterial')

        principled = nodes.new('ShaderNodeBsdfPrincipled')
        for o in options:
            principled.inputs[o].default_value = options[o]
        links.new(principled.outputs['BSDF'], material_output.inputs['Surface'])

        arrange_nodes(nodes, "socket")

        return mat

    def make_diffuse_material(self, name, color, roughness=0.3):
        mat = bpy.data.materials.new(name)
        mat.use_nodes = True
        mat.diffuse_color = color
        nodes = mat.node_tree.nodes
        nodes.clear()
        links = mat.node_tree.links

        material_output = nodes.new('ShaderNodeOutputMaterial')

        diffuse = nodes.new('ShaderNodeBsdfDiffuse')
        diffuse.inputs['Color'].default_value = color
        diffuse.inputs['Roughness'].default_value = roughness

        links.new(diffuse.outputs['BSDF'], material_output.inputs['Surface'])

        arrange_nodes(nodes, "socket")

        return mat

    def atom_material(self, name, element):
        """Return atom material"""
        color = self.element_color(element)
        return self.make_diffuse_material(name, color)

    def bond_material(self, name, bond):
        """Return bond material"""
        return self.make_diffuse_material(name, self.bond_color, roughness=1.0)

    def ring_material(self, name):
        """Return ring material"""
        mat = self.make_principled_material(name,
            base_color=self.ring_color,
            clearcoat=1.0,
            clearcoat_roughness=0.2,
            transmission=1.0,
            transmission_roughness=0.2,
            IOR=1.25)
        mat.blend_method = 'BLEND'
        mat.use_screen_refraction = True
        mat.refraction_depth = 0.1
        mat.shadow_method = 'CLIP'

        return mat

    def charge_material(self, pname, mname, element):
        """Return charge material"""
        pmat = self.make_diffuse_material(pname, self.charge_plus_color)
        mmat = self.make_diffuse_material(mname, self.charge_minus_color)

        return pmat, mmat

    def isosurface_material(self, isoname):
        """Return isosurface material"""
        color = self.iso_plus_color if "plus" in isoname else self.iso_minus_color
        return self.make_principled_material(isoname, base_color=color,
                metallic=0.7, specular=1.0, roughness=0.25, sheen=1.0)

    def outer_isosurface_material(self, isoname):
        """Return isosurface material"""
        base_color = self.iso_plus_color if "plus" in isoname else self.iso_minus_color

        mat = bpy.data.materials.new(isoname)
        mat.use_nodes = True
        mat.blend_method = 'BLEND'
        mat.diffuse_color = base_color
        nodes = mat.node_tree.nodes
        nodes.clear() # remove defaults
        links = mat.node_tree.links

        material_output = nodes.new('ShaderNodeOutputMaterial')

        mix = nodes.new('ShaderNodeMixShader')
        links.new(mix.outputs['Shader'], material_output.inputs['Surface'])

        layerweight = nodes.new('ShaderNodeLayerWeight')
        layerweight.inputs['Blend'].default_value = 0.7

        power = nodes.new('ShaderNodeMath')
        power.operation = 'POWER'
        power.inputs[1].default_value = 8.0
        links.new(layerweight.outputs['Facing'], power.inputs[0])
        links.new(power.outputs['Value'], mix.inputs['Fac'])

        transparent = nodes.new('ShaderNodeBsdfTransparent')
        links.new(transparent.outputs['BSDF'], mix.inputs[1])

        emission = nodes.new('ShaderNodeEmission')
        emission.inputs[0].default_value = base_color
        emission.inputs[1].default_value = 1.0
        links.new(emission.outputs['Emission'], mix.inputs[2])

        # arrange nodes
        arrange_nodes(nodes, "socket")

        return mat

    def element_color(self, element):
        """Returns RGB triple for element"""
        r, g, b = self.palette[element.symbol]
        a = 1.0
        return (r, g, b, a)

class DefaultElementStyler(PaletteElementStyler):
    """Color elements with Molecular Blender defaults"""

    def __init__(self):
        """Construct empty colorizer"""
        self.palette = {
            'h': (1.000, 1.000, 1.000),
            'he': (0.851, 1.000, 1.000),
            'li': (0.800, 0.502, 1.000),
            'be': (0.761, 1.000, 0.000),
            'b': (1.000, 0.710, 0.710),
            'c': (0.210, 0.210, 0.210),
            'n': (0.188, 0.314, 0.973),
            'o': (1.000, 0.051, 0.051),
            'f': (0.565, 0.878, 0.314),
            'ne': (0.702, 0.890, 0.961),
            'na': (0.671, 0.361, 0.949),
            'mg': (0.541, 1.000, 0.000),
            'al': (0.749, 0.651, 0.651),
            'si': (0.941, 0.784, 0.627),
            'p': (1.000, 0.502, 0.000),
            's': (1.000, 1.000, 0.188),
            'cl': (0.122, 0.941, 0.122),
            'ar': (0.502, 0.820, 0.890),
            'k': (0.561, 0.251, 0.831),
            'ca': (0.239, 1.000, 0.000),
            'sc': (0.902, 0.902, 0.902),
            'ti': (0.749, 0.761, 0.780),
            'v': (0.651, 0.651, 0.671),
            'cr': (0.541, 0.600, 0.780),
            'mn': (0.612, 0.478, 0.780),
            'fe': (0.878, 0.400, 0.200),
            'co': (0.941, 0.565, 0.627),
            'ni': (0.314, 0.816, 0.314),
            'cu': (0.784, 0.502, 0.200),
            'zn': (0.490, 0.502, 0.690),
            'ga': (0.761, 0.561, 0.561),
            'ge': (0.400, 0.561, 0.561),
            'as': (0.741, 0.502, 0.890),
            'se': (1.000, 0.631, 0.000),
            'br': (0.651, 0.161, 0.161),
            'kr': (0.361, 0.722, 0.820),
            'rb': (0.439, 0.180, 0.690),
            'sr': (0.000, 1.000, 0.000),
            'y': (0.580, 1.000, 1.000),
            'zr': (0.580, 0.878, 0.878),
            'nb': (0.451, 0.761, 0.788),
            'mo': (0.329, 0.710, 0.710),
            'tc': (0.231, 0.620, 0.620),
            'ru': (0.141, 0.561, 0.561),
            'rh': (0.039, 0.490, 0.549),
            'pd': (0.000, 0.412, 0.522),
            'ag': (0.753, 0.753, 0.753),
            'cd': (1.000, 0.851, 0.561),
            'in': (0.651, 0.459, 0.451),
            'sn': (0.400, 0.502, 0.502),
            'sb': (0.620, 0.388, 0.710),
            'te': (0.831, 0.478, 0.000),
            'i': (0.580, 0.000, 0.580),
            'xe': (0.259, 0.620, 0.690),
            'cs': (0.341, 0.090, 0.561),
            'ba': (0.000, 0.788, 0.000),
            'la': (0.439, 0.831, 1.000),
            'ce': (1.000, 1.000, 0.780),
            'pr': (0.851, 1.000, 0.780),
            'nd': (0.780, 1.000, 0.780),
            'pm': (0.639, 1.000, 0.780),
            'sm': (0.561, 1.000, 0.780),
            'eu': (0.380, 1.000, 0.780),
            'gd': (0.271, 1.000, 0.780),
            'tb': (0.188, 1.000, 0.780),
            'dy': (0.122, 1.000, 0.780),
            'ho': (0.000, 1.000, 0.612),
            'er': (0.000, 0.902, 0.459),
            'tm': (0.000, 0.831, 0.322),
            'yb': (0.000, 0.749, 0.220),
            'lu': (0.000, 0.671, 0.141),
            'hf': (0.302, 0.761, 1.000),
            'ta': (0.302, 0.651, 1.000),
            'w': (0.129, 0.580, 0.839),
            're': (0.149, 0.490, 0.671),
            'os': (0.149, 0.400, 0.588),
            'ir': (0.090, 0.329, 0.529),
            'pt': (0.816, 0.816, 0.878),
            'au': (1.000, 0.820, 0.137),
            'hg': (0.722, 0.722, 0.816),
            'tl': (0.651, 0.329, 0.302),
            'pb': (0.341, 0.349, 0.380),
            'bi': (0.620, 0.310, 0.710),
            'po': (0.671, 0.361, 0.000),
            'at': (0.459, 0.310, 0.271),
            'rn': (0.259, 0.510, 0.588),
            'fr': (0.259, 0.000, 0.400),
            'ra': (0.000, 0.490, 0.000),
            'ac': (0.439, 0.671, 0.980),
            'th': (0.000, 0.729, 1.000),
            'pa': (0.000, 0.631, 1.000),
            'u': (0.000, 0.561, 1.000),
            'np': (0.000, 0.502, 1.000),
            'pu': (0.000, 0.420, 1.000),
            'am': (0.329, 0.361, 0.949),
            'cm': (0.471, 0.361, 0.890),
            'bk': (0.541, 0.310, 0.890),
            'cf': (0.631, 0.212, 0.831),
            'es': (0.702, 0.122, 0.831),
            'fm': (0.702, 0.122, 0.729),
            'md': (0.702, 0.051, 0.651),
            'no': (0.741, 0.051, 0.529),
            'lr': (0.780, 0.000, 0.400),
            'rf': (0.800, 0.000, 0.349),
            'db': (0.820, 0.000, 0.310),
            'sg': (0.851, 0.000, 0.271),
            'bh': (0.878, 0.000, 0.220),
            'hs': (0.902, 0.000, 0.180),
            'mt': (0.922, 0.000, 0.149)
        }

VMD_COLORS = {
    "blue": (0.000000, 0.000000, 1.000000),
    "red": (1.000000, 0.000000, 0.000000),
    "gray": (0.350000, 0.350000, 0.350000),
    "orange": (1.000000, 0.500000, 0.000000),
    "yellow": (1.000000, 1.000000, 0.000000),
    "tan": (0.500000, 0.500000, 0.200000),
    "silver": (0.600000, 0.600000, 0.600000),
    "green": (0.000000, 1.000000, 0.000000),
    "white": (1.000000, 1.000000, 1.000000),
    "pink": (1.000000, 0.600000, 0.600000),
    "cyan": (0.250000, 0.750000, 0.750000),
    "purple": (0.650000, 0.000000, 0.650000),
    "lime": (0.500000, 0.900000, 0.400000),
    "mauve": (0.900000, 0.400000, 0.700000),
    "ochre": (0.500000, 0.300000, 0.000000),
    "iceblue": (0.500000, 0.500000, 0.750000),
    "black": (0.000000, 0.000000, 0.000000),
    "yellow2": (0.880000, 0.970000, 0.020000),
    "yellow3": (0.550000, 0.900000, 0.020000),
    "green2": (0.000000, 0.900000, 0.040000),
    "green3": (0.000000, 0.900000, 0.500000),
    "cyan2": (0.000000, 0.880000, 1.000000),
    "cyan3": (0.000000, 0.760000, 1.000000),
    "blue2": (0.020000, 0.380000, 0.670000),
    "blue3": (0.010000, 0.040000, 0.930000),
    "violet": (0.270000, 0.000000, 0.980000),
    "violet2": (0.450000, 0.000000, 0.900000),
    "magenta": (0.900000, 0.000000, 0.900000),
    "magenta2": (1.000000, 0.000000, 0.660000),
    "red2": (0.980000, 0.000000, 0.230000),
    "red3": (0.810000, 0.000000, 0.000000),
    "orange2": (0.890000, 0.350000, 0.000000),
    "orange3": (0.960000, 0.720000, 0.000000)
}

class VMDElementStyler(PaletteElementStyler):
    """Color elements with Molecular Blender defaults"""

    def __init__(self):
        """Construct empty colorizer"""
        self.palette = {
            "ac": VMD_COLORS["ochre"],
            "ag": VMD_COLORS["ochre"],
            "al": VMD_COLORS["ochre"],
            "am": VMD_COLORS["ochre"],
            "ar": VMD_COLORS["ochre"],
            "as": VMD_COLORS["ochre"],
            "at": VMD_COLORS["ochre"],
            "au": VMD_COLORS["ochre"],
            "b": VMD_COLORS["ochre"],
            "ba": VMD_COLORS["ochre"],
            "be": VMD_COLORS["ochre"],
            "bh": VMD_COLORS["ochre"],
            "bi": VMD_COLORS["ochre"],
            "bk": VMD_COLORS["ochre"],
            "br": VMD_COLORS["ochre"],
            "c": VMD_COLORS["cyan"],
            "ca": VMD_COLORS["ochre"],
            "cd": VMD_COLORS["ochre"],
            "ce": VMD_COLORS["ochre"],
            "cf": VMD_COLORS["ochre"],
            "cl": VMD_COLORS["ochre"],
            "cm": VMD_COLORS["ochre"],
            "co": VMD_COLORS["ochre"],
            "cr": VMD_COLORS["ochre"],
            "cs": VMD_COLORS["ochre"],
            "cu": VMD_COLORS["ochre"],
            "db": VMD_COLORS["ochre"],
            "ds": VMD_COLORS["ochre"],
            "dy": VMD_COLORS["ochre"],
            "er": VMD_COLORS["ochre"],
            "es": VMD_COLORS["ochre"],
            "eu": VMD_COLORS["ochre"],
            "f": VMD_COLORS["ochre"],
            "fe": VMD_COLORS["ochre"],
            "fm": VMD_COLORS["ochre"],
            "fr": VMD_COLORS["ochre"],
            "ga": VMD_COLORS["ochre"],
            "gd": VMD_COLORS["ochre"],
            "ge": VMD_COLORS["ochre"],
            "h": VMD_COLORS["white"],
            "he": VMD_COLORS["ochre"],
            "hf": VMD_COLORS["ochre"],
            "hg": VMD_COLORS["ochre"],
            "ho": VMD_COLORS["ochre"],
            "hs": VMD_COLORS["ochre"],
            "i": VMD_COLORS["ochre"],
            "in": VMD_COLORS["ochre"],
            "ir": VMD_COLORS["ochre"],
            "k": VMD_COLORS["ochre"],
            "kr": VMD_COLORS["ochre"],
            "la": VMD_COLORS["ochre"],
            "li": VMD_COLORS["ochre"],
            "lr": VMD_COLORS["ochre"],
            "lu": VMD_COLORS["ochre"],
            "md": VMD_COLORS["ochre"],
            "mg": VMD_COLORS["ochre"],
            "mn": VMD_COLORS["ochre"],
            "mo": VMD_COLORS["ochre"],
            "mt": VMD_COLORS["ochre"],
            "n": VMD_COLORS["blue"],
            "na": VMD_COLORS["ochre"],
            "nb": VMD_COLORS["ochre"],
            "nd": VMD_COLORS["ochre"],
            "ne": VMD_COLORS["ochre"],
            "ni": VMD_COLORS["ochre"],
            "no": VMD_COLORS["ochre"],
            "np": VMD_COLORS["ochre"],
            "o": VMD_COLORS["red"],
            "os": VMD_COLORS["ochre"],
            "p": VMD_COLORS["tan"],
            "pa": VMD_COLORS["ochre"],
            "pb": VMD_COLORS["ochre"],
            "pd": VMD_COLORS["ochre"],
            "pm": VMD_COLORS["ochre"],
            "po": VMD_COLORS["ochre"],
            "pr": VMD_COLORS["ochre"],
            "pt": VMD_COLORS["ochre"],
            "pu": VMD_COLORS["ochre"],
            "ra": VMD_COLORS["ochre"],
            "rb": VMD_COLORS["ochre"],
            "re": VMD_COLORS["ochre"],
            "rf": VMD_COLORS["ochre"],
            "rg": VMD_COLORS["ochre"],
            "rh": VMD_COLORS["ochre"],
            "rn": VMD_COLORS["ochre"],
            "ru": VMD_COLORS["ochre"],
            "s": VMD_COLORS["yellow"],
            "sb": VMD_COLORS["ochre"],
            "sc": VMD_COLORS["ochre"],
            "se": VMD_COLORS["ochre"],
            "sg": VMD_COLORS["ochre"],
            "si": VMD_COLORS["ochre"],
            "sm": VMD_COLORS["ochre"],
            "sn": VMD_COLORS["ochre"],
            "sr": VMD_COLORS["ochre"],
            "ta": VMD_COLORS["ochre"],
            "tb": VMD_COLORS["ochre"],
            "tc": VMD_COLORS["ochre"],
            "te": VMD_COLORS["ochre"],
            "th": VMD_COLORS["ochre"],
            "ti": VMD_COLORS["ochre"],
            "tl": VMD_COLORS["ochre"],
            "tm": VMD_COLORS["ochre"],
            "u": VMD_COLORS["ochre"],
            "v": VMD_COLORS["ochre"],
            "w": VMD_COLORS["ochre"],
            "x": VMD_COLORS["purple"],
            "xe": VMD_COLORS["ochre"],
            "y": VMD_COLORS["ochre"],
            "yb": VMD_COLORS["ochre"],
            "zn": VMD_COLORS["silver"],
            "zr": VMD_COLORS["ochre"]
        }

def get_styler(options):
    """Returns styler given option set"""
    colors = options["colors"]
    if "vmd" in colors:
        return VMDElementStyler()
    return DefaultElementStyler()
