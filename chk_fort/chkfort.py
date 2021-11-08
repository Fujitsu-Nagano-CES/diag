import numpy as np
from read_fort_cnt import read_fort_cnt

nx, ny, nz, nv, nm = (6, 3, 4, 6, 3)
npy, npz, npv, npm, nps = (2, 2, 2, 2, 2)


path1='../cnt/gkvp.000001.cnt.004'
time1, cnt1 = read_fort_cnt(path1, nx, ny, nz, nv, nm)
print(path1)
print('time=', time1)
print('shape=', cnt1.shape)
print('cnt0=', cnt1[0,0,0,0,0])
print('min=', cnt1.min())
print('max=', cnt1.max())

path2='../run/chgres_cnt/gkvp.000001.cnt.004'
time2, cnt2 = read_fort_cnt(path2, nx, ny, nz, nv, nm)
print(path2)
print('time=', time2)
print('cnt0=', cnt2[0,0,0,0,0])
print('shape=', cnt2.shape)
print('min=', cnt2.min())
print('max=', cnt2.max())

