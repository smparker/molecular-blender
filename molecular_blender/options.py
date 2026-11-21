"""Import option declarations shared by the Python API, scene records and Blender RNA."""
import math
import numbers
from dataclasses import dataclass


UNSET = object()  # An omitted option inherits the selected field's sign policy.


@dataclass(frozen=True)
class Option:
    default: object
    ui: dict
    choices: tuple


def option(default, *, ui_default=None, api_extra=(), **ui):
    """A UI label exposes this setting in the import operator. UI limits constrain the widget."""
    choices = tuple(item[0] for item in ui.get('items', ())) + api_extra
    if ui:
        ui['default'] = default if ui_default is None else ui_default
    return Option(default, ui, choices)


DENSITIES = (('density', 'Total density', 'Total electron density'),
             ('alpha', 'Alpha density', 'Alpha electron density'),
             ('beta', 'Beta density', 'Beta electron density'),
             ('spin', 'Spin density', 'Alpha minus beta density'))
ORBITAL_SPINS = (('all', 'All', 'All orbitals'), ('alpha', 'Alpha', 'Alpha frontier orbitals'),
                 ('beta', 'Beta', 'Beta frontier orbitals'))
MOLDEN_FIELDS = DENSITIES + (('orbital', 'Molecular orbital', 'Plot positive and negative orbital lobes'),)

OPTIONS = {
    'bonds': option(True),
    'bond_thickness': option(0.2, name='Thickness of bonds', ui_default=0.15,
        description='Determine overall thickness of bonds (0.0 turns bonds off)', min=0.0, max=50.0, step=0.05, precision=4),
    'draw_bond_orders': option(False, name='Draw bond orders',
        description='Draw supplied bond orders, or perceive them from coordinates (requires all hydrogens)'),
    'perceive_bond_orders': option(False),
    'molecular_charge': option(0, name='Molecular charge',
        description='Total integer charge for bond-order perception (not partial atomic charges)'),
    'bond_order_max_iterations': option(100000),
    'hook_atoms': option('auto', name='Hook bonds',
        description='Hook bonds to atoms to make bonds follow manual manipulation of the atoms',
        items=(
            ('auto', 'Automatic', 'Hook bonds for animated trajectories and Molden vibrational modes'),
            ('on', 'On', 'Hook bonds (slow for molecules with hundreds or more atoms)'),
            ('off', 'Off', 'Do not hook bonds to atoms (faster, but manipulating atoms in' + ' blender will not move the bonds alongside)'))),
    'plot_style': option('sticks', name='Style', ui_default='fixedbs',
        description='Plot style',
        items=(
            ('vdw', 'Van der Waals', 'Plot using VDW radii'),
            ('bs', 'Ball-and-stick', 'Ball and stick with bond radii determined automaticaly'),
            ('fixedbs', 'Fixed ball-and-stick', 'Ball and stick with a fixed bond radius'),
            ('sticks', 'Stick model', 'Just sticks'))),
    'plot_type': option('auto', api_extra=('auto', 'frame', 'animate')),
    'object_type': option('mesh', name='Object type',
        description='Object type to use to draw atoms',
        items=(
            ('meta', 'Metaballs', 'Metaballs'),
            ('nurbs', 'NURBS', 'NURBS'),
            ('mesh', 'Mesh', 'Mesh'),
            ('wireframe', 'Wireframe', 'Wireframe'))),
    'keystride': option(1, name='Keystride',
        description='Striding between keyframes in animation'),
    'trajectory_interpolation': option('blender', name='Playback',
        items=(
            ('blender', 'Blender Default', 'Use the keyframe interpolation from Blender preferences'),
            ('hold', 'Hold', 'Show each sampled geometry until the next snapshot'),
            ('linear', 'Linear', 'Interpolate atom positions and charge sizes linearly between snapshots'))),
    'frame_start': option(0, name='First input frame',
        description='Zero-based position in the XYZ trajectory or JSON entry list', min=0),
    'frame_end': option(-1, name='Last input frame',
        description='Inclusive input frame, or -1 to read through the end', min=-1),
    'frame_step': option(1, name='Read every Nth frame',
        description='Sample starting at the first selected input frame and preserve source timing', min=1),
    'animate_bonds': option('staticfirst', name='Bonds',
        description='Determine how the bonds are handled for an animation',
        items=(
            ('staticfirst', 'Static: first frame', 'Use bonds from only the first frame'),
            ('staticall', 'Static: all frames', 'Draw all bonds that are formed during any frame'),
            ('dynamic', 'Dynamically draw frames', 'Dynamically form and break bonds during animation'))),
    'bond_break_margin': option(5.0, name='Break margin (%)',
        description='Extra distance above the formation cutoff before an inferred bond disappears', min=0, soft_max=25),
    'universal_bonds': option(True, name='Universal bonds',
        description='Use one bond type for whole plot'),
    'ignore_hydrogen': option(True, name='Ignore Hydrogens', ui_default=False,
        description='Ignores Hydrogen atoms for cleaner images'),
    'gradient': option(False, name='Plot gradient',
        description='Draw arrows for data found in gradients column'),
    'charges': option('none', name='Plot charges',
        description='Style in which to plot atomic charges',
        items=(
            ('none', 'none', 'No charges plotted'),
            ('scale', 'scale', 'Charge magnitude encoded by scale of sphere'))),
    'charge_offset': option(0.0, name='chgoff', ui_default=0.9,
        description='Use chgfac*(charge + chgoff) to control visibility of charges', min=-5.0, max=5.0, step=0.05, precision=4),
    'charge_factor': option(1.0, name='chgfac',
        description='Use chgfac*(charge + chgoff) to control visibility of charges', min=-100.0, max=100.0, step=0.05, precision=2),
    'find_aromatic': option(False, name='Plot Aromatics',
        description='Find closed rings and if planar, fill in with object'),
    'recycle_materials': option(False, name='Recycle materials', ui_default=True,
        description='Re-use materials generated for previously imported molecules'),
    'isovalues': option('0.25,0.50', name='Isovalues to plot',
        description='Comma-separated surface levels. Leave empty to import just the structure'),
    'cumulative': option(True, name='Cumulative isovalues',
        description='Set isovalues so they contain given proportion of the total density'),
    'volume': option('orbital', name='Field',
        description='Field to calculate from the Molden wavefunction',
        items=MOLDEN_FIELDS, api_extra=('response',)),
    'cube_both_signs': option(UNSET, name='Include both signs', ui_default=True,
        description='Generate both signs for each cube level. Turn off to use exactly the signs entered'),
    'cube_dataset': option(0, name='Cube dataset index',
        description='Zero-based dataset index within a cube file', min=0),
    'cube_field_type': option('auto', name='Cube field type',
        description='How cumulative cube surfaces integrate values',
        items=(
            ('auto', 'Auto / scalar magnitude', 'Use file metadata; untyped fields integrate absolute magnitude'),
            ('orbital', 'Orbital amplitude', 'Integrate squared amplitude'),
            ('density', 'Density', 'Integrate absolute density'),
            ('signed', 'Signed scalar', 'Integrate absolute field magnitude'))),
    'orbital': option('homo', name='Orbital',
        description='Plot this orbital from molden file. Accepts comma separated list, homo-n, lumo+m. (ignored for other inputs)'),
    'orbital_spin': option('all', name='Orbital Spin',
        description='Channel for frontier-orbital names; explicit numbers always use file indices',
        items=ORBITAL_SPINS),
    'response_eigenpair': option(0, name='DIPL eigenpair',
        description='Eigenpair ID in the DIPL file, or zero for the first listed pair', min=0),
    'align_orbital_phase': option(False, name='Align orbital phase across frames',
        description='Use the overlap with the preceding displayed orbital to keep lobe colors consistent'),
    'remesh': option(6, name='Remesh Octree Depth', ui_default=8,
        description='Octree depth for remesh modifier'),
    'resolution': option(0.05, name='Resolution',
        description='Desired spacing between points in isosurface (Angstrom)', min=0.001, max=0.5, step=0.01, precision=3),
    'colors': option('default', name='Colors',
        description='Palette of colors to use',
        items=(
            ('default', 'default', 'Colors defined by molecular blender'),
            ('vmd', 'vmd', 'Same colors defined by VMD Element'))),
    'molden_convention': option('auto', name='Molden convention',
        items=(
            ('auto', 'Auto', 'Use ORCA for an orca_2mkl producer title, otherwise standard Molden'),
            ('standard', 'Standard Molden', 'Use standard Molden normalization and phases'),
            ('orca', 'ORCA', 'Use orca_2mkl normalization and phases'))),
    'sdf_record': option(0),
    'atom_scale': option(1.0),
    'surface_opacity': option(1.0),
}

DEFAULTS = {key: spec.default for key, spec in OPTIONS.items() if spec.default is not UNSET}
ENUMS = {key: spec.choices for key, spec in OPTIONS.items() if spec.choices}
IMPORT_PROPERTIES = tuple(key for key, spec in OPTIONS.items() if spec.ui and spec.default is not UNSET)


def molden_convention(options):
    """Translate the legacy boolean in scripts and saved scene/surface settings."""
    if 'orca_fix' in options and not isinstance(options['orca_fix'], bool):
        raise ValueError('orca_fix must be a boolean')
    if 'molden_convention' not in options and 'orca_fix' in options:
        return 'orca' if options['orca_fix'] else 'standard'
    mode = options.get('molden_convention', 'auto')
    if mode not in ENUMS['molden_convention']:
        raise ValueError('molden_convention must be auto, standard, or orca')
    return mode


def trajectory_range(options):
    start, end, step = (options.get("frame_start", 0), options.get("frame_end", -1),
                        options.get("frame_step", 1))
    if any(type(v) is not int for v in (start, end, step)) or start < 0 or step < 1 or (end != -1 and end < start):
        raise ValueError("Frame range needs a nonnegative start, an inclusive end >= start (or -1), and a positive step")
    return start, end, step


def validated_options(values):
    """Return a fresh settings block without mutating the caller's choices.

    Resolution is in Angstrom. Isovalues and hook policy remain user-facing
    strings; process_options converts a separate copy for numerical/scene work.
    Historical Python API defaults are retained independently of UI presets.
    """
    out = dict(DEFAULTS, **values)
    out['molden_convention'] = molden_convention(values)
    out.pop('orca_fix', None)
    if 'cube_both_signs' in out and not isinstance(out['cube_both_signs'], bool):
        raise ValueError('cube_both_signs must be a boolean')
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
    for name in ('bond_thickness', 'bond_break_margin', 'cube_dataset', 'sdf_record', 'remesh', 'response_eigenpair'):
        if out[name] < 0:
            raise ValueError(name + ' must not be negative')
    trajectory_range(out)
    if out['atom_scale'] <= 0:
        raise ValueError('Atom size multiplier must be positive')
    if not 0 <= out['surface_opacity'] <= 1:
        raise ValueError('Surface opacity must be between zero and one')
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


def default_options(options):
    """add defaults to options in case called outside of blender"""
    options['molden_convention'] = molden_convention(options)
    for name, value in DEFAULTS.items():
        options.setdefault(name, value)
    return options


def process_options(options, *, cube=False):
    """postprocess choices that might interact with each other"""

    hooking = {"on": True, "off": False,
               "auto": options["plot_type"] == "animate" }
    options["hook_atoms"] = hooking[options["hook_atoms"]]

    options["isosurfaces"] = "isovalues" in options and options["isovalues"] != ""
    if options["isosurfaces"]:
        isovals = options["isovalues"].split(',')
        isovals = [float(x) for x in isovals]
        if not isovals:
            raise Exception("Isovalues must be specified")

        both_signs = options["volume"] in ("orbital", "spin", "response")
        if cube:
            both_signs = options.get("cube_both_signs", both_signs)
        if both_signs:
            # automatically include positive and negative
            iso_set = set()
            for iso in isovals:
                iso_set.add(iso)
                iso_set.add(-iso)
            isovals = list(iso_set)

        options["isovalues"] = isovals

    return options


def get_frame_mapper(options, frame_numbers=None):
    """Map a snapshot position to a Blender frame, including explicit manifest timing."""
    stride = options["keystride"]
    def f(i):
        return (frame_numbers[i] if frame_numbers is not None else i) * stride + 1
    return f
