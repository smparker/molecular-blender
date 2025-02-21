import os
import re
import sys
import time
from .generic import splitPath, changeFileExtension
import numpy as np

from Cython.Build import cythonize

def compile_cython(addon_dir):
    cython_files = [ os.path.join(addon_dir, f) for f in os.listdir(addon_dir) if f.endswith('.pyx') ]
    for f in cython_files:
        cythonize(f, compiler_directives = {"language_level" : "3"}, include_path = [ np.get_include() ])

def compile_cxx(addon_dir):
    print("Compile:\n{:40s}".format('-'))

    pyx_files = [ os.path.join(addon_dir, f) for f in os.listdir(addon_dir) if f.endswith('.pyx') ]
    cxx_files = [ f.replace('.pyx', '.cpp') for f in pyx_files ]

    for cxx in cxx_files:
        ext = getExtensionFromPath(cxx, addon_dir)
        buildExtensionInplace(ext)

def getExtensionFromPath(path, addonDirectory, includeDirs = [np.get_include()]):
    from distutils.core import Extension
    moduleName = getModuleNameOfPath(path, addonDirectory)
    print("moduleName: ", moduleName)

    kwargs = {
        "sources" : [path],
        "include_dirs" : includeDirs,
        "define_macros" : [],
        "undef_macros" : [],
        "library_dirs" : [],
        "libraries" : [],
        "runtime_library_dirs" : [],
        "extra_objects" : [],
        "extra_compile_args" : [ "-march=native", "-ffast-math" ],
        "extra_link_args" : [],
        "export_symbols" : [],
        "depends" : []
    }

    for key, values in getExtensionArgsFromSetupOptions(getSetupOptions(path)).items():
        kwargs[key].extend(values)

    infoFile = changeFileExtension(path, "_setup_info.py")
    for key, values in getExtensionsArgsFromInfoFile(infoFile).items():
        kwargs[key].extend(values)

    return Extension(moduleName, **kwargs)

def getModuleNameOfPath(path, basePath):
    relativePath = os.path.relpath(os.path.splitext(path)[0], os.path.dirname(basePath))
    return ".".join(splitPath(relativePath))

def getExtensionsArgsFromInfoFile(infoFilePath):
    if not os.path.isfile(infoFilePath):
        return {}

    data = executePythonFile(infoFilePath)
    fName = "getExtensionArgs"
    if fName not in data:
        return {}

    return data[fName](Utils)

def buildExtensionInplace(extension):
    from distutils.core import setup
    oldArgs = sys.argv
    sys.argv = [oldArgs[0], "build_ext", "--inplace" ]
    setup(ext_modules = [extension])
    sys.argv = oldArgs

def getSetupOptions(path):
    pyxPath = changeFileExtension(path, ".pyx")
    if not os.path.isfile(pyxPath):
        return set()

    options = set()
    with open(pyxPath, "rt") as f:
        text = f.read()
    for match in re.finditer(r"^#\s*setup\s*:\s*options\s*=(.*)$", text, flags = re.MULTILINE):
        options.update(match.group(1).split())
    return options

def getExtensionArgsFromSetupOptions(options):
    args = {}
    if "c++11" in options:
        if onLinux or onMacOS:
            args["extra_compile_args"] = ["-std=c++11"]
    return args
