'''
Setup originally created by JacquesLucke
Modified for molecular_blender by Shane Parker 2019-2024
'''

import os
import sys
import subprocess
import argparse
import importlib
import shutil

currentDirectory = os.path.dirname(os.path.abspath(__file__))

if not os.path.samefile(currentDirectory, os.getcwd()):
    print("You are not in the correct directory.")
    print("Expected:", currentDirectory)
    print("Got:     ", os.getcwd())
    raise RuntimeError("setup.py should be run from the directory it is in.")

from _setuputils.generic import *
from _setuputils.setup_info_files import getSetupInfoList
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
    parser.add_argument("--mode", "-m", choices=['link', 'copy'], default='copy',
        help='install by either linking of copying')
    parser.add_argument("--check-version", action='store_true', help='check python version')
    parser.add_argument("--link-target", default='all',
        help='Blender versions to install into, comma-separated (default: all found)')
    parser.add_argument("--dry-run", action='store_true', help='preview install destinations without changing files')

    args = parser.parse_args()

    if args.command == 'build':
        build(check_environment=args.check_version)
    elif args.command == 'install':
        install(install_mode=args.mode, link_target=args.link_target, dry_run=args.dry_run)
    elif args.command == 'export':
        export()

# Build
####################################################

def build(check_environment=False):
    if check_environment:
        check_build_environment()

    changedFileStates = build_impl()
    printChangedFileStates(changedFileStates, currentDirectory)

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

def export():
    execute_Export(addonDirectory, exportPath, addonName)

def printChangedFileStates(states, basepath):
    print()
    print('Change Summary')
    for p in states['new']:
        print(' + {}'.format(p))
    for p in states['removed']:
        print(' - {}'.format(p))
    for p in states['changed']:
        print(' ~ {}'.format(p))

@returnChangedFileStates(currentDirectory)
def build_impl():
    from _setuputils.compilation import compile_cxx, compile_cython
    compile_cython(addonDirectory)
    setupInfoList = getSetupInfoList(addonDirectory)

    compile_cxx(addonDirectory)

def check_build_environment(cython=True, python=True):
    fail = False
    if cython and (importlib.util.find_spec("Cython") is None):
        print("Cython is not installed for this Python version.")
        print(sys.version)
        fail = True
    if python and not (sys.version_info.major, sys.version_info.minor) == (3, 10):
        print("Blender 3.6 uses Python 3.10")
        print(f"You are using {sys.version}")
        fail = True
    if fail:
        raise Exception("Build environment is not correct.")

if __name__ == "__main__":
    main()

