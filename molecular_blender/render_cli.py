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
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--output', type=Path, help='Destination directory for transparent PNGs.')
    action.add_argument('--list', action='store_true', help='List the camera, molecule names, surface labels and group IDs.')
    parser.add_argument('--surface', action='append', default=[], help='Exact sidebar label or group ID. Repeat to select several.')
    parser.add_argument('--legends', action='store_true', help='Write matching SVG field legends for mapped PNGs.')
    try:
        args = parser.parse_args(argv)
        if args.list and args.surface:
            parser.error('--surface selects panels to render with --output')
        if args.list and args.legends:
            parser.error('--legends requires --output')
    except SystemExit as exc:
        return exc.code
    try:
        if not bpy.data.filepath:
            raise ValueError('Pass a saved .blend file before -c molecular_render')
        scene = bpy.context.scene
        owners = {obj.name: obj for obj in scene.objects if RECORD in obj}
        if (args.molecule is not None and args.molecule not in owners
                or args.molecule is None and not args.list and len(owners) != 1):
            raise ValueError('Choose --molecule from the active scene: ' + (', '.join(sorted(owners)) or '(none)'))
        if args.molecule is not None:
            owners = {args.molecule: owners[args.molecule]}
        if args.list:
            print('Camera:', scene.camera.name if scene.camera else '(none)')
            if not owners:
                print('No saved molecules in this scene.')
            for name, owner in sorted(owners.items()):
                print('Molecule:', name)
                groups = surface_groups(owner, read_record(owner))
                for key, group in groups.items():
                    print(f"  {group['label']}  ({key})")
                if not groups:
                    print('  (no saved surfaces)')
            return 0
        owner, = owners.values()
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
        export_surface_pngs(bpy.context, owner, groups, keys, directory, legends=args.legends)
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
