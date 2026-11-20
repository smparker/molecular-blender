"""Render the saved surface groups of a molecule from Blender's command line."""
import argparse
from pathlib import Path
import sys

import bpy

from .scene_state import RECORD, read_record
from .scene_surfaces import surface_groups
from .surface_export import export_surface_pngs

_command = None


def render_series(argv):
    parser = argparse.ArgumentParser(prog='blender -b figure.blend -c molecular_render',
                                     description='Render saved surface groups with the active scene camera.')
    parser.add_argument('--molecule', help='Exact molecule root object name. Optional for a single molecule.')
    parser.add_argument('--output', required=True, type=Path, help='Destination directory for transparent PNGs.')
    parser.add_argument('--surface', action='append', default=[], help='Exact sidebar label or group ID. Repeat to select several.')
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        return exc.code
    try:
        if not bpy.data.filepath:
            raise ValueError('Pass a saved .blend file before -c molecular_render')
        owners = {obj.name: obj for obj in bpy.context.scene.objects if RECORD in obj}
        if args.molecule is None and len(owners) == 1:
            owner, = owners.values()
        else:
            owner = owners.get(args.molecule)
        if owner is None:
            raise ValueError('Choose --molecule from the active scene: ' + (', '.join(sorted(owners)) or '(none)'))
        groups = surface_groups(owner, read_record(owner))
        keys = []
        for name in args.surface:
            matches = [k for k, g in groups.items() if name == k or name == g['label']]
            if len(matches) != 1:
                choices = ', '.join(f"{g['label']} ({k})" for k, g in groups.items())
                raise ValueError(f'Surface {name!r} needs one match. Available surfaces: {choices}')
            if matches[0] not in keys:
                keys.append(matches[0])
        directory = args.output.expanduser().resolve()
        directory.mkdir(parents=True, exist_ok=True)
        keys = keys or list(groups)
        export_surface_pngs(bpy.context, owner, groups, keys, directory)
        print(f'Exported {len(keys)} PNG(s) to {directory}')
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        print(f'molecular_render: {exc}', file=sys.stderr)
        return 1
    return 0


def register():
    global _command
    if hasattr(bpy.utils, 'register_cli_command'):
        _command = bpy.utils.register_cli_command('molecular_render', render_series)


def unregister():
    global _command
    if _command is not None:
        bpy.utils.unregister_cli_command(_command)
        _command = None
