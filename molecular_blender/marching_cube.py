CYTHON_ENABLED = False

try:
    from .marching_cube_cy import *
    CYTHON_ENABLED = True
except:
    print("Warning: using pure python isosurface implementation.")
    print("Building the cython version could give up to factor of 50 speedup.")
    from .marching_cube_py import *
