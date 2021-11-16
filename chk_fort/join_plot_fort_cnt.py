import numpy as np
import matplotlib.pyplot as plt
import sys, os, glob
from read_fort_cnt import read_fort_cnt

nx, ny, nz, nv, nm = (6, 3, 4, 6, 3)
npy, npz, npv, npm, nps = (2, 2, 2, 2, 2)
dir = '../cnt'


def join_fort_cnt(dir,  nx, ny, nz, nv, nm,  npy, npz, npv, npm, nps):
    path_lst = glob.glob(os.path.join(dir, '*'))
    file_lst = [os.path.basename(x) for x in path_lst]
    sfx_set = {x[-3:] for x in file_lst}
    stp_lst = []
    for s in sfx_set:
        if s.isdigit():
            stp_lst.append(s)
    last_stp_files = []
    for fn in file_lst:
        if fn.endswith(stp_lst[-1]):
            last_stp_files.append(fn)
    last_stp_files.sort()
    n_proc = npy * npz * npv * npm * nps
    if n_proc != len(last_stp_files):
        print('Error: invalid process num.')
        return
    print('the last step is {}'.format(stp_lst[-1]))

    # check box size
    gny, gnz, gnv, gnm = (ny*npy, nz*npz, nv*npv, nm*(npm+1)-1)
    wcnt = np.zeros((nx*2+1, gny+1, 2*gnz, 2*gnv, gnm+1), dtype=np.complex128)
    hx, hy, hz, hv, hm \
        = (int(nx/2), int(gny/2), int(gnz/2), int(gnv/2), int(gnm/2))
    irank = 0
    for ips in range(nps):
        for ipm in range(npm):
            for ipv in range(npv):
                for ipz in range(npz):
                    for ipy in range(npy):
                        tm, pcnt = read_fort_cnt(os.path.join(dir,
                                                      last_stp_files[irank]),
                                                 nx, ny, nz, nv, nm)
                        if ipy == npy-1:
                            yedg = gny+1
                            yedl = gny+1 - (ny+1)*(npy-1)
                        else:
                            yedg = (ny+1)*(ipy+1)
                            yedl = ny+1
                        wcnt[:,
                             (ny+1)*ipy:yedg,
                             2*nz*ipz:2*nz*(ipz+1),
                             2*nv*ipv:2*nv*(ipv+1),
                             (nm+1)*ipm:(nm+1)*(ipm+1)] = pcnt[:,:yedl,:,:,:]
                        irank = irank + 1

        # plot real part
        fig, axes = plt.subplots(nrows=4, ncols=4, sharex=False)
        # row#0
        axes[0,0].imshow(wcnt[:, :, hz, hv, hm].real);axes[0,0].set_title("X-Y")
        axes[0,1].imshow(wcnt[:, hy, :, hv, hm].real);axes[0,1].set_title("X-Z")
        axes[0,2].imshow(wcnt[:, hy, hz, :, hm].real);axes[0,2].set_title("X-V")
        axes[0,3].imshow(wcnt[:, hy, hz, hv, :].real);axes[0,3].set_title("X-M")
        # row#1
        axes[1,0].axis('off')
        axes[1,1].imshow(wcnt[hx, :, :, hv, hm].real);axes[1,1].set_title("Y-Z")
        axes[1,2].imshow(wcnt[hx, :, hz, :, hm].real);axes[1,2].set_title("Y-V")
        axes[1,3].imshow(wcnt[hx, :, hz, hv, :].real);axes[1,3].set_title("Y-M")
        # row#2
        axes[2,0].axis('off')
        axes[2,1].axis('off')
        axes[2,2].imshow(wcnt[hx, hy, :, :, hm].real);axes[2,2].set_title("Z-V")
        axes[2,3].imshow(wcnt[hx, hy, :, hv, :].real);axes[2,3].set_title("Z-M")
        # row#3
        axes[3,0].axis('off'); axes[3,0].set_title("REAL")
        axes[3,1].axis('off'); axes[3,1].set_title("S={}".format(ips))
        axes[3,2].axis('off')
        axes[3,3].imshow(wcnt[hx, hy, hz, :, :].real);axes[3,3].set_title("V-M")
        fig.tight_layout()
        plt.show()

        # plot imag part
        fig, axes = plt.subplots(nrows=4, ncols=4, sharex=False)
        # row#0
        axes[0,0].imshow(wcnt[:, :, hz, hv, hm].imag);axes[0,0].set_title("X-Y")
        axes[0,1].imshow(wcnt[:, hy, :, hv, hm].imag);axes[0,1].set_title("X-Z")
        axes[0,2].imshow(wcnt[:, hy, hz, :, hm].imag);axes[0,2].set_title("X-V")
        axes[0,3].imshow(wcnt[:, hy, hz, hv, :].imag);axes[0,3].set_title("X-M")
        # row#1
        axes[1,0].axis('off')
        axes[1,1].imshow(wcnt[hx, :, :, hv, hm].imag);axes[1,1].set_title("Y-Z")
        axes[1,2].imshow(wcnt[hx, :, hz, :, hm].imag);axes[1,2].set_title("Y-V")
        axes[1,3].imshow(wcnt[hx, :, hz, hv, :].imag);axes[1,3].set_title("Y-M")
        # row#2
        axes[2,0].axis('off')
        axes[2,1].axis('off')
        axes[2,2].imshow(wcnt[hx, hy, :, :, hm].imag);axes[2,2].set_title("Z-V")
        axes[2,3].imshow(wcnt[hx, hy, :, hv, :].imag);axes[2,3].set_title("Z-M")
        # row#3
        axes[3,0].axis('off'); axes[3,0].set_title("IMAG")
        axes[3,1].axis('off'); axes[3,1].set_title("S={}".format(ips))
        axes[3,2].axis('off')
        axes[3,3].imshow(wcnt[hx, hy, hz, :, :].imag);axes[3,3].set_title("V-M")
        fig.tight_layout()
        plt.show()
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='join gkvp cnt fort files')
    parser.add_argument('-d', type=str, default='../cnt')
    parser.add_argument('--nx', type=int, default=nx)
    parser.add_argument('--ny', type=int, default=ny)
    parser.add_argument('--nz', type=int, default=nz)
    parser.add_argument('--nv', type=int, default=nv)
    parser.add_argument('--nm', type=int, default=nm)
    parser.add_argument('--npy', type=int, default=npy)
    parser.add_argument('--npz', type=int, default=npz)
    parser.add_argument('--npv', type=int, default=npv)
    parser.add_argument('--npm', type=int, default=npm)
    parser.add_argument('--nps', type=int, default=nps)
    args = parser.parse_args()

    if args.d: dir = args.d
    if args.nx: nx = args.nx
    if args.ny: ny = args.ny
    if args.nz: nz = args.nz
    if args.nv: nv = args.nv
    if args.nm: nm = args.nm
    if args.ny: npy = args.npy
    if args.nz: npz = args.npz
    if args.nv: npv = args.npv
    if args.nm: npm = args.npm
    if args.nm: nps = args.nps
    
    print('chek {}'.format(dir))
    join_fort_cnt(dir, nx, ny, nz, nv, nm, npy, npz, npv, npm, nps)
    sys.exit(0)
    
    
