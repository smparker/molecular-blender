try:
    from .orbitals_cy import *
except:
    print("Warning: using pure python isosurface implementation.")
    print("Building the cython version could give up to factor of 50 speedup.")
    from .orbitals_py import *
