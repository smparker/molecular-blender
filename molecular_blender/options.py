"""Serializable import settings, separate from processed runtime options."""
import math
import numbers

DEFAULTS = {"bonds": True,
                "bond_thickness": 0.2,
                "draw_bond_orders": False,
                "perceive_bond_orders": False,
                "molecular_charge": 0,
                "bond_order_max_iterations": 100000,
                "hook_atoms": "auto",
                "plot_style": "sticks",
                "plot_type": "auto",
                "object_type": "mesh",
                "keystride": 1,
                "animate_bonds": "staticfirst",
                "universal_bonds": True,
                "ignore_hydrogen": True,
                "gradient": False,
                "charges": "none",
                "charge_offset": 0.0,
                "charge_factor": 1.0,
                "find_aromatic": False,
                "recycle_materials": False,
                "isovalues": "0.25,0.50",
                "cumulative": True,
                "volume": "orbital",
                "cube_dataset": 0,
                "cube_field_type": "auto",
                "orbital": "homo",
                "orbital_spin": "all",
                "remesh" : 6,
                "resolution": 0.05,
                "colors": "default"
               }


# Additional reader controls also belong to saved settings.
DEFAULTS.update(orca_fix=False, sdf_record=0)
ENUMS = {
    'hook_atoms': ('auto', 'on', 'off'),
    'plot_style': ('vdw', 'bs', 'fixedbs', 'sticks'),
    'plot_type': ('auto', 'frame', 'animate'),
    'object_type': ('mesh', 'nurbs', 'meta', 'wireframe'),
    'animate_bonds': ('staticfirst', 'staticall', 'dynamic'),
    'charges': ('none', 'scale'),
    'volume': ('orbital', 'density', 'alpha', 'beta', 'spin'),
    'cube_field_type': ('auto', 'orbital', 'density', 'signed'),
    'orbital_spin': ('all', 'alpha', 'beta'),
    'colors': ('default', 'vmd'),
}


def validated_options(values):
    """Return a fresh settings block without mutating the caller's choices.

    Resolution is in Angstrom. Isovalues and hook policy remain user-facing
    strings; process_options converts a separate copy for numerical/scene work.
    Historical Python API defaults are retained independently of UI presets.
    """
    out = dict(DEFAULTS, **values)
    for name, choices in ENUMS.items():
        if out[name] not in choices:
            raise ValueError('{} must be one of {}'.format(name, ', '.join(choices)))
    for name, default in DEFAULTS.items():
        value = out[name]
        if isinstance(default, bool):
            if not isinstance(value, bool):
                raise ValueError(name + ' must be a boolean')
        elif isinstance(default, int):
            if isinstance(value, bool) or not isinstance(value, numbers.Integral):
                raise ValueError(name + ' must be an integer')
            out[name] = int(value)
        elif isinstance(default, float):
            if isinstance(value, bool) or not isinstance(value, numbers.Real) or not math.isfinite(value):
                raise ValueError(name + ' must be a finite number')
            out[name] = float(value)
        elif not isinstance(value, str):
            raise ValueError(name + ' must be text')
    for name in ('resolution', 'keystride', 'bond_order_max_iterations'):
        if out[name] <= 0:
            raise ValueError(name + ' must be positive')
    for name in ('bond_thickness', 'cube_dataset', 'sdf_record', 'remesh'):
        if out[name] < 0:
            raise ValueError(name + ' must not be negative')
    if out['isovalues'].strip():
        try:
            levels = [float(v) for v in out['isovalues'].split(',')]
        except ValueError as exc:
            raise ValueError('Isovalues must be comma-separated numbers') from exc
        if not all(math.isfinite(v) for v in levels):
            raise ValueError('Isovalues must be finite')
        if out['cumulative'] and any(not 0 < abs(v) <= 1 for v in levels):
            raise ValueError('Cumulative isovalues must have magnitude in (0, 1]')
    return out
