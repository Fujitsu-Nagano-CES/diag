import pdb
import numpy as np
import netCDF4

nx, ny, nz, nv, nm = (6, 3, 4, 6, 3)
npy, npz, npv, npm, nps = (2, 2, 2, 2, 2)
# dims = (2*nx+1, ny+1, 2*nz, 2*nv, nm+1)
# gdims = (2*nx+1, gny+1, 2*gnz, 2*gnv, gnm+1)

file1 = '../run/data/cnt.004.nc'
file2 = '../run/chgres_cnt/gkvp.cnt.004.nc'

cnt1 = netCDF4.Dataset(file1, 'r')
cnt2 = netCDF4.Dataset(file2, 'r')

cnt1r = cnt1['recnt'][0]
cnt2r = cnt2['recnt'][0]
cnt_rd = cnt1r - cnt2r
print('diff real min=%g, max=%g' %(cnt_rd.min(), cnt_rd.max()))
cnt1i = cnt1['imcnt'][0]
cnt2i = cnt2['imcnt'][0]
cnt_id = cnt1i - cnt2i
print('diff imag min=%g, max=%g' %(cnt_id.min(), cnt_id.max()))
