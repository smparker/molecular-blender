#distutils: language = c++
#cython: profile=True
#cython: cdivision=True

cimport cython

import numpy as np
cimport numpy as np
from numpy cimport ndarray

from libcpp.vector cimport vector as vector
from libc.math cimport sqrt, exp, fabs

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

DEF CYDEBUG = False

cdef int *edgetable=[0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
                                0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
                                0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
                                0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
                                0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
                                0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
                                0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
                                0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
                                0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
                                0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
                                0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
                                0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
                                0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
                                0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
                                0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
                                0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
                                0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
                                0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
                                0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
                                0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
                                0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
                                0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
                                0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
                                0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
                                0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
                                0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
                                0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
                                0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
                                0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
                                0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
                                0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
                                0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0]


cdef int tritable[256][16]
tritable[0][:] = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[1][:] = [0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[2][:] = [0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[3][:] = [1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[4][:] = [1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[5][:] = [0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[6][:] = [9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[7][:] = [2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1]
tritable[8][:] = [3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[9][:] = [0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[10][:] = [1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[11][:] = [1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1]
tritable[12][:] = [3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[13][:] = [0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1]
tritable[14][:] = [3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1]
tritable[15][:] = [9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[16][:] = [4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[17][:] = [4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[18][:] = [0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[19][:] = [4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1]
tritable[20][:] = [1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[21][:] = [3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1]
tritable[22][:] = [9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1]
tritable[23][:] = [2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1]
tritable[24][:] = [8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[25][:] = [11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1]
tritable[26][:] = [9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1]
tritable[27][:] = [4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1]
tritable[28][:] = [3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1]
tritable[29][:] = [1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1]
tritable[30][:] = [4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1]
tritable[31][:] = [4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1]
tritable[32][:] = [9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[33][:] = [9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[34][:] = [0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[35][:] = [8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1]
tritable[36][:] = [1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[37][:] = [3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1]
tritable[38][:] = [5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1]
tritable[39][:] = [2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1]
tritable[40][:] = [9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[41][:] = [0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1]
tritable[42][:] = [0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1]
tritable[43][:] = [2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1]
tritable[44][:] = [10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1]
tritable[45][:] = [4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1]
tritable[46][:] = [5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1]
tritable[47][:] = [5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1]
tritable[48][:] = [9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[49][:] = [9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1]
tritable[50][:] = [0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1]
tritable[51][:] = [1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[52][:] = [9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1]
tritable[53][:] = [10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1]
tritable[54][:] = [8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1]
tritable[55][:] = [2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1]
tritable[56][:] = [7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1]
tritable[57][:] = [9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1]
tritable[58][:] = [2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1]
tritable[59][:] = [11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1]
tritable[60][:] = [9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1]
tritable[61][:] = [5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1]
tritable[62][:] = [11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1]
tritable[63][:] = [11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[64][:] = [10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[65][:] = [0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[66][:] = [9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[67][:] = [1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1]
tritable[68][:] = [1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[69][:] = [1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1]
tritable[70][:] = [9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1]
tritable[71][:] = [5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1]
tritable[72][:] = [2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[73][:] = [11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1]
tritable[74][:] = [0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1]
tritable[75][:] = [5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1]
tritable[76][:] = [6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1]
tritable[77][:] = [0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1]
tritable[78][:] = [3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1]
tritable[79][:] = [6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1]
tritable[80][:] = [5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[81][:] = [4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1]
tritable[82][:] = [1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1]
tritable[83][:] = [10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1]
tritable[84][:] = [6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1]
tritable[85][:] = [1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1]
tritable[86][:] = [8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1]
tritable[87][:] = [7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1]
tritable[88][:] = [3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1]
tritable[89][:] = [5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1]
tritable[90][:] = [0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1]
tritable[91][:] = [9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1]
tritable[92][:] = [8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1]
tritable[93][:] = [5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1]
tritable[94][:] = [0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1]
tritable[95][:] = [6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1]
tritable[96][:] = [10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[97][:] = [4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1]
tritable[98][:] = [10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1]
tritable[99][:] = [8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1]
tritable[100][:] = [1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1]
tritable[101][:] = [3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1]
tritable[102][:] = [0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[103][:] = [8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1]
tritable[104][:] = [10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1]
tritable[105][:] = [0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1]
tritable[106][:] = [3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1]
tritable[107][:] = [6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1]
tritable[108][:] = [9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1]
tritable[109][:] = [8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1]
tritable[110][:] = [3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1]
tritable[111][:] = [6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[112][:] = [7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1]
tritable[113][:] = [0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1]
tritable[114][:] = [10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1]
tritable[115][:] = [10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1]
tritable[116][:] = [1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1]
tritable[117][:] = [2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1]
tritable[118][:] = [7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1]
tritable[119][:] = [7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[120][:] = [2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1]
tritable[121][:] = [2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1]
tritable[122][:] = [1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1]
tritable[123][:] = [11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1]
tritable[124][:] = [8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1]
tritable[125][:] = [0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[126][:] = [7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1]
tritable[127][:] = [7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[128][:] = [7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[129][:] = [3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[130][:] = [0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[131][:] = [8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1]
tritable[132][:] = [10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[133][:] = [1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1]
tritable[134][:] = [2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1]
tritable[135][:] = [6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1]
tritable[136][:] = [7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[137][:] = [7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1]
tritable[138][:] = [2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1]
tritable[139][:] = [1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1]
tritable[140][:] = [10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1]
tritable[141][:] = [10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1]
tritable[142][:] = [0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1]
tritable[143][:] = [7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1]
tritable[144][:] = [6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[145][:] = [3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1]
tritable[146][:] = [8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1]
tritable[147][:] = [9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1]
tritable[148][:] = [6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1]
tritable[149][:] = [1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1]
tritable[150][:] = [4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1]
tritable[151][:] = [10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1]
tritable[152][:] = [8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1]
tritable[153][:] = [0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[154][:] = [1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1]
tritable[155][:] = [1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1]
tritable[156][:] = [8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1]
tritable[157][:] = [10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1]
tritable[158][:] = [4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1]
tritable[159][:] = [10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[160][:] = [4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[161][:] = [0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1]
tritable[162][:] = [5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1]
tritable[163][:] = [11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1]
tritable[164][:] = [9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1]
tritable[165][:] = [6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1]
tritable[166][:] = [7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1]
tritable[167][:] = [3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1]
tritable[168][:] = [7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1]
tritable[169][:] = [9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1]
tritable[170][:] = [3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1]
tritable[171][:] = [6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1]
tritable[172][:] = [9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1]
tritable[173][:] = [1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1]
tritable[174][:] = [4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1]
tritable[175][:] = [7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1]
tritable[176][:] = [6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1]
tritable[177][:] = [3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1]
tritable[178][:] = [0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1]
tritable[179][:] = [6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1]
tritable[180][:] = [1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1]
tritable[181][:] = [0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1]
tritable[182][:] = [11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1]
tritable[183][:] = [6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1]
tritable[184][:] = [5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1]
tritable[185][:] = [9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1]
tritable[186][:] = [1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1]
tritable[187][:] = [1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[188][:] = [1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1]
tritable[189][:] = [10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1]
tritable[190][:] = [0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[191][:] = [10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[192][:] = [11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[193][:] = [11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1]
tritable[194][:] = [5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1]
tritable[195][:] = [10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1]
tritable[196][:] = [11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1]
tritable[197][:] = [0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1]
tritable[198][:] = [9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1]
tritable[199][:] = [7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1]
tritable[200][:] = [2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1]
tritable[201][:] = [8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1]
tritable[202][:] = [9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1]
tritable[203][:] = [9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1]
tritable[204][:] = [1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[205][:] = [0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1]
tritable[206][:] = [9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1]
tritable[207][:] = [9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[208][:] = [5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1]
tritable[209][:] = [5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1]
tritable[210][:] = [0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1]
tritable[211][:] = [10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1]
tritable[212][:] = [2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1]
tritable[213][:] = [0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1]
tritable[214][:] = [0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1]
tritable[215][:] = [9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[216][:] = [2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1]
tritable[217][:] = [5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1]
tritable[218][:] = [3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1]
tritable[219][:] = [5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1]
tritable[220][:] = [8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1]
tritable[221][:] = [0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[222][:] = [8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1]
tritable[223][:] = [9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[224][:] = [4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1]
tritable[225][:] = [0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1]
tritable[226][:] = [1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1]
tritable[227][:] = [3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1]
tritable[228][:] = [4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1]
tritable[229][:] = [9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1]
tritable[230][:] = [11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1]
tritable[231][:] = [11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1]
tritable[232][:] = [2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1]
tritable[233][:] = [9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1]
tritable[234][:] = [3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1]
tritable[235][:] = [1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[236][:] = [4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1]
tritable[237][:] = [4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1]
tritable[238][:] = [4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[239][:] = [4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[240][:] = [9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[241][:] = [3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1]
tritable[242][:] = [0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1]
tritable[243][:] = [3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[244][:] = [1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1]
tritable[245][:] = [3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1]
tritable[246][:] = [0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[247][:] = [3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[248][:] = [2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1]
tritable[249][:] = [9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[250][:] = [2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1]
tritable[251][:] = [1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[252][:] = [1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[253][:] = [0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[254][:] = [0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
tritable[255][:] = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]

cdef struct point:
    float x
    float y
    float z

cdef struct triangle:
    point p0
    point p1
    point p2

cdef struct box:
    point p0
    point p1

ctypedef vector[triangle] triangle_vec
ctypedef vector[triangle_vec] iso_triangles

@cython.boundscheck(CYDEBUG)
@cython.wraparound(CYDEBUG)
@cython.initializedcheck(CYDEBUG)
@cython.nonecheck(CYDEBUG)
cdef void polygonise(triangle_vec* tri_list, const DTYPE_t* cornervalues,
        const DTYPE_t* cornerpos,
        float isolevel) nogil:
    global edgetable, tritable

    cdef point[12] vertlist

    #   Determine the index into the edge table which
    #   tells us which vertices are inside of the surface
    cdef int cubeindex = 0
    if cornervalues[0] < isolevel: cubeindex = cubeindex | 1
    if cornervalues[1] < isolevel: cubeindex = cubeindex | 2
    if cornervalues[2] < isolevel: cubeindex = cubeindex | 4
    if cornervalues[3] < isolevel: cubeindex = cubeindex | 8
    if cornervalues[4] < isolevel: cubeindex = cubeindex | 16
    if cornervalues[5] < isolevel: cubeindex = cubeindex | 32
    if cornervalues[6] < isolevel: cubeindex = cubeindex | 64
    if cornervalues[7] < isolevel: cubeindex = cubeindex | 128

    # Cube is entirely in/out of the surface
    if edgetable[cubeindex] == 0: return

    if (edgetable[cubeindex] & 1):    vertex_interp(&vertlist[ 0], isolevel, &cornerpos[0*3], &cornerpos[1*3], cornervalues[0], cornervalues[1])
    if (edgetable[cubeindex] & 2):    vertex_interp(&vertlist[ 1], isolevel, &cornerpos[1*3], &cornerpos[2*3], cornervalues[1], cornervalues[2])
    if (edgetable[cubeindex] & 4):    vertex_interp(&vertlist[ 2], isolevel, &cornerpos[2*3], &cornerpos[3*3], cornervalues[2], cornervalues[3])
    if (edgetable[cubeindex] & 8):    vertex_interp(&vertlist[ 3], isolevel, &cornerpos[3*3], &cornerpos[0*3], cornervalues[3], cornervalues[0])
    if (edgetable[cubeindex] & 16):   vertex_interp(&vertlist[ 4], isolevel, &cornerpos[4*3], &cornerpos[5*3], cornervalues[4], cornervalues[5])
    if (edgetable[cubeindex] & 32):   vertex_interp(&vertlist[ 5], isolevel, &cornerpos[5*3], &cornerpos[6*3], cornervalues[5], cornervalues[6])
    if (edgetable[cubeindex] & 64):   vertex_interp(&vertlist[ 6], isolevel, &cornerpos[6*3], &cornerpos[7*3], cornervalues[6], cornervalues[7])
    if (edgetable[cubeindex] & 128):  vertex_interp(&vertlist[ 7], isolevel, &cornerpos[7*3], &cornerpos[4*3], cornervalues[7], cornervalues[4])
    if (edgetable[cubeindex] & 256):  vertex_interp(&vertlist[ 8], isolevel, &cornerpos[0*3], &cornerpos[4*3], cornervalues[0], cornervalues[4])
    if (edgetable[cubeindex] & 512):  vertex_interp(&vertlist[ 9], isolevel, &cornerpos[1*3], &cornerpos[5*3], cornervalues[1], cornervalues[5])
    if (edgetable[cubeindex] & 1024): vertex_interp(&vertlist[10], isolevel, &cornerpos[2*3], &cornerpos[6*3], cornervalues[2], cornervalues[6])
    if (edgetable[cubeindex] & 2048): vertex_interp(&vertlist[11], isolevel, &cornerpos[3*3], &cornerpos[7*3], cornervalues[3], cornervalues[7])

    #Create the triangle
    cdef int i=0
    cdef triangle tri_dummy
    while tritable[cubeindex][i] != -1:
        tri_dummy.p0 = vertlist[tritable[cubeindex][i  ]]
        tri_dummy.p1 = vertlist[tritable[cubeindex][i+1]]
        tri_dummy.p2 = vertlist[tritable[cubeindex][i+2]]
        tri_list.push_back(tri_dummy)
        i+=3


@cython.boundscheck(CYDEBUG)
@cython.wraparound(CYDEBUG)
@cython.initializedcheck(CYDEBUG)
@cython.nonecheck(CYDEBUG)
cdef void vertex_interp(point* vert_out, float isolevel, const DTYPE_t* p1, const DTYPE_t* p2, const float v1, const float v2) nogil:
    if (fabs(isolevel-v1) < 0.00001):
        vert_out.x = p1[0]
        vert_out.y = p1[1]
        vert_out.z = p1[2]
        return
    if (fabs(isolevel-v2) < 0.00001):
        vert_out.x = p2[0]
        vert_out.y = p2[1]
        vert_out.z = p2[2]
        return
    if (fabs(v1-v2) < 0.00001):
        vert_out.x = p1[0]
        vert_out.y = p1[1]
        vert_out.z = p1[2]
        return

    cdef float mu = (isolevel - v1) / (v2 - v1);
    vert_out.x = p1[0] + mu * (p2[0] - p1[0]);
    vert_out.y = p1[1] + mu * (p2[1] - p1[1]);
    vert_out.z = p1[2] + mu * (p2[2] - p1[2]);

    return


@cython.boundscheck(CYDEBUG)
@cython.wraparound(CYDEBUG)
@cython.initializedcheck(CYDEBUG)
@cython.nonecheck(CYDEBUG)
cdef bint active_box(const DTYPE_t* cornervalues, const int niso, const DTYPE_t* isolevels) nogil:
    """Determines whether any isovalues are crossed based on given corner values"""
    global edgetable

    cdef int iiso
    cdef int cubeindex
    cdef float isolevel
    for iiso in range(niso):
        isolevel = isolevels[iiso]
        cubeindex = 0

        if cornervalues[0] < isolevel: cubeindex |= 1
        if cornervalues[1] < isolevel: cubeindex |= 2
        if cornervalues[2] < isolevel: cubeindex |= 4
        if cornervalues[3] < isolevel: cubeindex |= 8
        if cornervalues[4] < isolevel: cubeindex |= 16
        if cornervalues[5] < isolevel: cubeindex |= 32
        if cornervalues[6] < isolevel: cubeindex |= 64
        if cornervalues[7] < isolevel: cubeindex |= 128

        if edgetable[cubeindex] != 0:
            return True

    return False


@cython.boundscheck(CYDEBUG)
@cython.wraparound(CYDEBUG)
@cython.initializedcheck(CYDEBUG)
@cython.nonecheck(CYDEBUG)
cpdef marching_cube_box(ndarray[DTYPE_t, ndim=3, mode="c"] data, object isovalues):
    cdef int nx = data.shape[0]
    cdef int ny = data.shape[1]
    cdef int nz = data.shape[2]

    cdef int niso = len(isovalues)
    cdef int iiso
    cdef vector[float] isovals
    for iiso in range(niso):
        isovals.push_back(isovalues[iiso])

    cdef ndarray[DTYPE_t, ndim=1, mode="c"] cornervalues = np.zeros(8, dtype=DTYPE)
    cdef ndarray[DTYPE_t, ndim=2, mode="c"] cornerpos = np.zeros([8,3], dtype=DTYPE)

    cdef iso_triangles c_triangles = iso_triangles(niso)

    cdef int ix, iy, iz
    cdef int ix2, iy2, iz2
    cdef float x, y, z
    cdef float x2, y2, z2

    for iz in range(nz-1):
        iz2 = iz + 1
        for ix in range(nx-1):
            ix2 = ix + 1
            for iy in range(ny-1):
                iy2 = iy + 1

                cornervalues[0] = data[ix ,iy ,iz ]
                cornervalues[1] = data[ix ,iy2,iz ]
                cornervalues[2] = data[ix2,iy2,iz ]
                cornervalues[3] = data[ix2,iy ,iz ]
                cornervalues[4] = data[ix ,iy ,iz2]
                cornervalues[5] = data[ix ,iy2,iz2]
                cornervalues[6] = data[ix2,iy2,iz2]
                cornervalues[7] = data[ix2,iy ,iz2]

                x = <float>(ix)
                y = <float>(iy)
                z = <float>(iz)

                x2 = <float>(ix2)
                y2 = <float>(iy2)
                z2 = <float>(iz2)

                cornerpos[0,0] =  x
                cornerpos[0,1] =  y
                cornerpos[0,2] =  z

                cornerpos[1,0] =  x
                cornerpos[1,1] = y2
                cornerpos[1,2] =  z

                cornerpos[2,0] = x2
                cornerpos[2,1] = y2
                cornerpos[2,2] =  z

                cornerpos[3,0] = x2
                cornerpos[3,1] =  y
                cornerpos[3,2] =  z

                cornerpos[4,0] =  x
                cornerpos[4,1] =  y
                cornerpos[4,2] = z2

                cornerpos[5,0] =  x
                cornerpos[5,1] = y2
                cornerpos[5,2] = z2

                cornerpos[6,0] = x2
                cornerpos[6,1] = y2
                cornerpos[6,2] = z2

                cornerpos[7,0] = x2
                cornerpos[7,1] =  y
                cornerpos[7,2] = z2

                for iiso in range(niso):
                    polygonise(&c_triangles[iiso], &cornervalues[0], &cornerpos[0,0], isovals[iiso])

    # now repackage in a python friendly format
    tri_list = [[] for iso in isovalues]
    cdef int ntriangles
    cdef int itri
    for iiso in range(niso):
        ntriangles = c_triangles[iiso].size()
        for itri in range(ntriangles):
            tri_tmp = c_triangles[iiso][itri]
            tri_list[iiso].append( [
                [ tri_tmp.p0.x, tri_tmp.p0.y, tri_tmp.p0.z ],
                [ tri_tmp.p1.x, tri_tmp.p1.y, tri_tmp.p1.z ],
                [ tri_tmp.p2.x, tri_tmp.p2.y, tri_tmp.p2.z ] ] )

    return tri_list


@cython.boundscheck(CYDEBUG)
@cython.wraparound(CYDEBUG)
@cython.initializedcheck(CYDEBUG)
@cython.nonecheck(CYDEBUG)
def marching_cube_outline(ndarray[DTYPE_t, ndim=3, mode="c"] data,
        ndarray[DTYPE_t, ndim=1, mode="c"] xvals,
        ndarray[DTYPE_t, ndim=1, mode="c"] yvals,
        ndarray[DTYPE_t, ndim=1, mode="c"] zvals,
        object isovalues):
    """
    Return [ (pp0, pp1) ] where pp0 and pp1 are lower and upper bounds
    of active boxes.
    """

    cdef int nx = data.shape[0]
    cdef int ny = data.shape[1]
    cdef int nz = data.shape[2]

    cdef int niso = len(isovalues)
    cdef int iiso
    cdef vector[float] isovals
    for iiso in range(niso):
        isovals.push_back(isovalues[iiso])

    cdef ndarray[DTYPE_t, ndim=1, mode="c"] cornervalues = np.zeros(8, dtype=DTYPE)

    outlines = []

    cdef int ix, iy, iz, ix2, iy2, iz2
    cdef float x, y, z, x2, y2, z2

    cdef vector[box] active_boxes
    cdef box box_tmp

    for iz in range(nz-1):
        iz2 = iz + 1
        for ix in range(nx-1):
            ix2 = ix + 1
            for iy in range(ny-1):
                iy2 = iy + 1

                cornervalues[0] = data[ix ,iy ,iz ]
                cornervalues[1] = data[ix ,iy2,iz ]
                cornervalues[2] = data[ix2,iy2,iz ]
                cornervalues[3] = data[ix2,iy ,iz ]
                cornervalues[4] = data[ix ,iy ,iz2]
                cornervalues[5] = data[ix ,iy2,iz2]
                cornervalues[6] = data[ix2,iy2,iz2]
                cornervalues[7] = data[ix2,iy ,iz2]

                x = xvals[ix]
                y = yvals[iy]
                z = zvals[iz]

                x2 = xvals[ix2]
                y2 = yvals[iy2]
                z2 = zvals[iz2]

                if active_box(&cornervalues[0], niso, &isovals[0]):
                    box_tmp.p0.x = x
                    box_tmp.p0.y = y
                    box_tmp.p0.z = z
                    box_tmp.p1.x = x2
                    box_tmp.p1.y = y2
                    box_tmp.p1.z = z2

                    active_boxes.push_back(box_tmp)

    outlines = []
    cdef int nboxes = active_boxes.size()
    cdef int ibox
    for ibox in range(nboxes):
        outlines.append( (
            (
                active_boxes[ibox].p0.x,
                active_boxes[ibox].p0.y,
                active_boxes[ibox].p0.z
            ),
            (
                active_boxes[ibox].p1.x,
                active_boxes[ibox].p1.y,
                active_boxes[ibox].p1.z
            )
            ) )
    return outlines
