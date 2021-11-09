import pdb
import math
import numpy as np
from read_fort_cnt import read_fort_cnt

nx, ny, nz, nv, nm = (6, 3, 4, 6, 3)
npy, npz, npv, npm, nps = (2, 2, 2, 2, 2)
# dims = (2*nx+1, ny+1, 2*nz, 2*nv, nm+1)
# gdims = (2*nx+1, gny+1, 2*gnz, 2*gnv, gnm+1)

nrank = 32
file_lst = ['gkvp.%06d.cnt.004' % r for r in range(nrank)]

dir1='../cnt/'
dir2='../run/chgres_cnt_fort/'

for file in file_lst:
    print('check ' + file)
    time1, cnt1 = read_fort_cnt(dir1+file, nx, ny, nz, nv, nm)
    time2, cnt2 = read_fort_cnt(dir2+file, nx, ny, nz, nv, nm)
    if time1 != time2:
        print(' time mismatch.')
        print('  time1=%f, time2=%d' %(time1, time2))
    cnt_rd = cnt1.real - cnt2.real
    cnt_id = cnt1.imag - cnt2.imag
    if abs(cnt_rd.min()) > 1e-17 or abs(cnt_rd.max()) > 1e-17 or \
       abs(cnt_id.min()) > 1e-17 or abs(cnt_id.max()) > 1e-17 :
        #pdb.set_trace()
        print(' cnt data mismatch.')
        print('  diff real min=%g max=%g' % (cnt_rd.min(), cnt_rd.max()))
        print('  diff imag min=%g max=%g' % (cnt_id.min(), cnt_id.max()))


