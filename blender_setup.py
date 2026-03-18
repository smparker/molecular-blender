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
from _setuputils.compilation import compile_cxx, compile_cython
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
    parser.add_argument("--link-target", default='all', help='which available blender versions to symlink (default: all found)')

    args = parser.parse_args()

    if args.command == 'build':
        build(check_environment=args.check_version)
    elif args.command == 'install':
        install(install_mode=args.mode)
    elif args.command == 'export':
        export()

# Build
####################################################

def build(check_environment=False):
    if check_environment:
        check_build_environment()

    changedFileStates = build_impl()
    printChangedFileStates(changedFileStates, currentDirectory)

def install(install_mode='copy'):
    if install_mode not in ('link', 'copy'):
        raise ValueError("install_mode must be either 'link' or 'copy'")
    paths = find_user_blender_dirs()

    source = os.path.join(currentDirectory, addonName)
    for p in paths:
        addon = os.path.join(p, 'scripts', 'addons')
        os.makedirs(addon, exist_ok=True)
        target = os.path.join(addon, addonName)

        if os.path.isdir(target):
            print("  - removing directory {}".format(target))
            if os.path.islink(target):
                os.unlink(target)
            else:
                shutil.rmtree(target)

        if install_mode == 'copy':
            print("  - copying directory {} to {}".format(source, target))
            shutil.copytree(source, target)
        else:
            print("  - creating symlink from {} to {}".format(source, target))
            os.symlink(source, target)

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

