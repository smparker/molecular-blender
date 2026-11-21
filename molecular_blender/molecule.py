# Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
# SPDX-License-Identifier: GPL-3.0-or-later

"""Historical model imports. Scientific code can always import molecular_data directly."""
from . import _HAS_BPY
from .molecular_data import VolumeData, bonded_pairs

if _HAS_BPY:
    from .scene_model import Atom, Bond, Molecule, Snapshot
else:
    from .molecular_data import Atom, Bond, Molecule, Snapshot
