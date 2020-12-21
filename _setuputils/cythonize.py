import os

import numpy as np
from Cython.Build import cythonize

def compile_cython(addon_dir):
    cython_files = [ os.path.join(addon_dir, f) for f in os.listdir(addon_dir) if f.endswith('.pyx') ]
    for f in cython_files:
        cythonize(f, compiler_directives = {"language_level" : "3"}, include_path = [ np.get_include() ])
