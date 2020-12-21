'''
Setup originally created by JacquesLucke
Modified for molecular_blender by Shane Parker 2019-2020
'''

import os
import sys
import json
import textwrap
import subprocess
from pprint import pprint

import argparse

currentDirectory = os.path.dirname(os.path.abspath(__file__))

if not os.path.samefile(currentDirectory, os.getcwd()):
    print("You are not in the correct directory.")
    print("Expected:", currentDirectory)
    print("Got:     ", os.getcwd())
    sys.exit()

if currentDirectory not in sys.path:
    sys.path.append(currentDirectory)

from _setuputils.generic import *
from _setuputils.addon_files import *
from _setuputils.compilation import compile_cxx, compile_cython
from _setuputils.setup_info_files import getSetupInfoList
from _setuputils.export import execute_Export

addonName = "molecular_blender"
addonDirectory = os.path.join(currentDirectory, addonName)

addonVersion = getAddonVersion(os.path.join(addonDirectory, "__init__.py"))
exportName = "{}_v{}_{}_{}_py{}{}".format(
    addonName, *addonVersion[:2], currentOS, *sys.version_info[:2])

exportPath = os.path.join(currentDirectory, exportName + ".zip")

# Main
####################################################
def main():
    parser = argparse.ArgumentParser("Setup molecular_blender")
    parser.add_argument("command", choices=['build', 'link', 'export'])
    parser.add_argument("--link", "-c", action='store_true', help='link build to available blender installations')
    parser.add_argument("--export", "-e", action='store_true', help='create installable .zip file')
    parser.add_argument("--no-check", action='store_true', help='skip python version check')
    parser.add_argument("--link-target", default='all', help='which available blender versions to symlink (default: all found)')

    args = parser.parse_args()

    if args.command == 'build':
        build()
    elif args.command == 'link':
        link(args.link_target)
    elif args.command == 'export':
        export()

# Build
####################################################

def build():
    check_build_environment()

    changedFileStates = build_impl()
    printChangedFileStates(changedFileStates, currentDirectory)

def link(link_target='all'):
    paths = find_user_blender_dirs()

    source = os.path.join(currentDirectory, addonName)
    for p in paths:
        addon = os.path.join(p, 'scripts', 'addons')
        os.makedirs(addon, exist_ok=True)
        target = os.path.join(addon, addonName)
        if not os.path.isdir(target):
            print("  - creating symlink at {}".format(target))
            os.symlink(source, target)
        else:
            print("  - skipping directory {}".format(target))

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

def printIndentedPathList(paths, basepath):
    if len(paths) == 0:
        print("  <none>")
    else:
        for path in sorted(paths):
            print("  {}".format(os.path.relpath(path, basepath)))

@returnChangedFileStates(currentDirectory)
def build_impl():
    compile_cython(addonDirectory)
    setupInfoList = getSetupInfoList(addonDirectory)

    compile_cxx(addonDirectory)

def have_cython():
    try:
        import Cython
        return True
    except:
        return False

def have_python_37():
    v = sys.version_info
    if v.major != 3 or v.minor != 7:
        return False
    return True

def check_build_environment(cython=True, python=True):
    fail = False
    if cython and not have_cython():
        print("Cython is not installed for this Python version.")
        print(sys.version)
        fail = True
    if python and not have_python_37():
        print(textwrap.dedent('''\
        Blender 2.8 officially uses Python 3.7.x.
        You are using: {}

        Use the --no-check option to disable this check.\
        '''.format(sys.version)))
        fail = True
    if fail:
        sys.exit()

# Run Main
###############################################

main()

