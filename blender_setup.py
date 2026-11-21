'''
Setup originally created by JacquesLucke
Modified for molecular_blender by Shane Parker 2019-2024
'''

import os
import sys
import subprocess
import argparse
import json
import platform
import sysconfig
import shutil

currentDirectory = os.path.dirname(os.path.abspath(__file__))

if not os.path.samefile(currentDirectory, os.getcwd()):
    print("You are not in the correct directory.")
    print("Expected:", currentDirectory)
    print("Got:     ", os.getcwd())
    raise RuntimeError("blender_setup.py should be run from the directory it is in.")

from _setuputils.generic import currentOS, find_user_blender_dirs, getAddonVersion
from _setuputils.export import execute_Export

addonName = "molecular_blender"
addonDirectory = os.path.join(currentDirectory, addonName)

addonVersion = getAddonVersion(os.path.join(addonDirectory, "__init__.py"))
exportName = "{}_v{}_{}_{}_py{}_{}".format(
    addonName, *addonVersion[:2], currentOS, *sys.version_info[:2])

exportPath = os.path.join(currentDirectory, exportName + ".zip")

# Main
####################################################
def main():
    parser = argparse.ArgumentParser("Setup molecular_blender")
    parser.add_argument("command", choices=['build', 'install', 'export'])
    parser.add_argument("--wheel-dir", help="export native modules from these wheels")
    parser.add_argument("--output", help="export ZIP destination")
    parser.add_argument("--mode", "-m", choices=['link', 'copy'], default='copy',
        help='install by copying or linking the add-on')
    parser.add_argument("--check-version", action='store_true', help='verify build ABI against BLENDER_BIN or blender on PATH')
    parser.add_argument("--blender", metavar="EXECUTABLE",
        help="verify the build interpreter against this Blender executable")
    parser.add_argument("--link-target", default='all',
        help='Blender versions to install into, comma-separated (default: all found)')
    parser.add_argument("--dry-run", action='store_true', help='preview install destinations without changing files')

    args = parser.parse_args()

    if args.command == 'build':
        build(args.blender, check_environment=args.check_version)
    elif args.command == 'install':
        install(install_mode=args.mode, link_target=args.link_target, dry_run=args.dry_run)
    elif args.command == 'export':
        export(args.wheel_dir, args.output)

# Build
####################################################

def build(blender=None, check_environment=False):
    if blender or check_environment:
        check_build_environment(blender or os.environ.get('BLENDER_BIN', 'blender'))
    subprocess.run([sys.executable, 'setup.py', 'build_ext', '--inplace', '--force'],
                   cwd=currentDirectory, check=True)


def check_build_environment(blender):
    marker = 'MB_BUILD_RUNTIME '
    probe = ("import json, platform, sysconfig, numpy, bpy; "
             "print('MB_BUILD_RUNTIME ' + json.dumps([bpy.app.version_string, numpy.__version__, "
             "sysconfig.get_config_var('EXT_SUFFIX'), platform.machine().lower()]))")
    output = subprocess.check_output([blender, '-b', '--factory-startup', '--python-exit-code', '1',
                                      '--python-expr', probe], text=True)
    record = next((line[len(marker):] for line in output.splitlines() if line.startswith(marker)), None)
    if record is None:
        raise ValueError(f'Could not identify the Python runtime in {blender}: {output.strip()}')
    version, numpy_version, suffix, machine = json.loads(record)
    print(f'Blender {version} uses {suffix} on {machine} with NumPy {numpy_version}', flush=True)
    build_suffix = sysconfig.get_config_var('EXT_SUFFIX')
    build_machine = platform.machine().lower()
    if (build_suffix, build_machine) != (suffix, machine):
        raise ValueError(f'Build interpreter uses {build_suffix} on {build_machine}. '
                         'Select a Python with the same ABI and architecture as Blender.')


def install(install_mode='copy', link_target='all', dry_run=False):
    if install_mode not in ('link', 'copy'):
        raise ValueError("install_mode must be either 'link' or 'copy'")
    paths = find_user_blender_dirs()
    if link_target != 'all':
        requested = {version.strip() for version in link_target.split(',')}
        available = {os.path.basename(path) for path in paths}
        missing = requested - available
        if missing:
            raise ValueError('Unknown Blender install target(s): {}. Available: {}'.format(
                ', '.join(sorted(missing)), ', '.join(sorted(available)) or 'none'))
        paths = [path for path in paths if os.path.basename(path) in requested]

    source = os.path.join(currentDirectory, addonName)
    targets = [os.path.join(path, 'scripts', 'addons', addonName) for path in paths]
    for target in targets:
        print('  - {} {} -> {}{}'.format(install_mode, source, target,
              ' (replaces existing installation)' if os.path.lexists(target) else ''))
        if os.path.lexists(target) and not (os.path.islink(target) or os.path.isdir(target)):
            raise FileExistsError('Install destination is not a directory or symlink: ' + target)
    if dry_run:
        return targets
    for target in targets:
        addon = os.path.dirname(target)
        os.makedirs(addon, exist_ok=True)

        if os.path.islink(target):
            os.unlink(target)
        elif os.path.isdir(target):
            print("  - removing directory {}".format(target))
            shutil.rmtree(target)

        if install_mode == 'copy':
            print("  - copying directory {} to {}".format(source, target))
            shutil.copytree(source, target)
        else:
            print("  - creating symlink from {} to {}".format(source, target))
            os.symlink(source, target)
    return targets

def export(wheel_dir=None, output=None):
    if wheel_dir:
        from _setuputils.release import export_wheels
        export_wheels(addonDirectory, output or exportPath, wheel_dir)
    else:
        execute_Export(addonDirectory, output or exportPath, addonName)

if __name__ == "__main__":
    main()

