### GNU Fortran compiler ###
FC = gfortran
FFLAGS = -mcmodel=medium -m64 -march=native -mtune=native -O3 -ffast-math # Optimization
FFLAGS += -frecord-marker=4 # 4-byte header/footer
FFLAGS += -fdefault-integer-8 # Change default integer type as integer(kind=8).
#FFLAGS += -fopenmp # OpenMP
#FFLAGS += -Wall -Wextra -pedantic -fbacktrace \
#         -fbounds-check -Wuninitialized -ffpe-trap=invalid,zero,overflow # Debug
### Usage of FFTW (spack load fftw)
##FFTW_DIR=`spack location -i fftw`
FFTW_DIR=/usr
INC += -I$(FFTW_DIR)/include
LIB += -L$(FFTW_DIR)/lib -lfftw3 -lm
### Usage of NetCDF (spack load netcdf-fortran)
##NETCDF_DIR=`spack location -i netcdf-fortran`
NETCDF_DIR=/usr
INC += -I$(NETCDF_DIR)/include
##LIB += -L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5
LIB += `nc-config --libs` `nc-config --flibs`


### Read Fortran binary GKV output
FILEIO = diag_rb_fortran
### Read NetCDF GKV output
#FILEIO = diag_rb_netcdf


SRC = ./src
PROG = diag.exe

diag:	${SRC}/diag_header.f90\
	${SRC}/diag_functions.f90\
	${SRC}/diag_clock.f90\
	${SRC}/gkvp_igs.f90\
	${SRC}/gkvp_vmecbzx.f90\
	${SRC}/diag_geom.f90\
	${SRC}/diag_intgrl.f90\
	${SRC}/diag_parity.f90\
	${SRC}/${FILEIO}.f90\
	${SRC}/diag_fft.f90\
	${SRC}/out_textfile.f90\
	${SRC}/out_netcdf.f90\
	${SRC}/out_mominxy.f90\
	${SRC}/out_mominxz.f90\
	${SRC}/out_mominkxky.f90\
	${SRC}/out_trninkxky.f90\
	${SRC}/out_triinkxky.f90\
	${SRC}/out_mominz.f90\
	${SRC}/out_momintky.f90\
	${SRC}/out_mominavs.f90\
	${SRC}/out_mominxmf.f90\
	${SRC}/out_mominvtk.f90\
	${SRC}/out_mominrz.f90\
	${SRC}/out_ffinkxky.f90\
	${SRC}/out_ffinvm.f90\
	${SRC}/out_ffinzv.f90\
	${SRC}/out_ffinzvm_vtk.f90\
	${SRC}/out_mominfreq.f90\
	${SRC}/out_linfreq.f90\
	${SRC}/out_correlate.f90\
	${SRC}/out_realsp.f90\
	${SRC}/out_bicoherence.f90\
	${SRC}/out_zfshearing.f90\
	${SRC}/out_zfdensity.f90\
	${SRC}/out_fluidtotaltrans.f90\
	${SRC}/out_fluidsubsptrans.f90\
	${SRC}/out_fluiddetailtrans.f90\
	${SRC}/diag_main.f90

	${FC} ${FFLAGS} -c ${SRC}/diag_header.f90
	${FC} ${FFLAGS} -c ${SRC}/diag_functions.f90
	${FC} ${FFLAGS} -c ${SRC}/diag_clock.f90
	${FC} ${FFLAGS} -c ${SRC}/gkvp_igs.f90
	${FC} ${FFLAGS} -c ${SRC}/gkvp_vmecbzx.f90
	${FC} ${FFLAGS} -c ${SRC}/diag_geom.f90
	${FC} ${FFLAGS} -c ${SRC}/diag_intgrl.f90
	${FC} ${FFLAGS} -c ${SRC}/diag_parity.f90
	${FC} ${FFLAGS} -c ${SRC}/${FILEIO}.f90 ${INC}
	${FC} ${FFLAGS} -c ${SRC}/diag_fft.f90 ${INC}
	${FC} ${FFLAGS} -c ${SRC}/out_textfile.f90
	${FC} ${FFLAGS} -c ${SRC}/out_netcdf.f90 ${INC}
	${FC} ${FFLAGS} -c ${SRC}/out_mominxy.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_mominxz.f90
	${FC} ${FFLAGS} -c ${SRC}/out_mominkxky.f90
	${FC} ${FFLAGS} -c ${SRC}/out_trninkxky.f90
	${FC} ${FFLAGS} -c ${SRC}/out_triinkxky.f90
	${FC} ${FFLAGS} -c ${SRC}/out_mominz.f90
	${FC} ${FFLAGS} -c ${SRC}/out_momintky.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_mominavs.f90
	${FC} ${FFLAGS} -c ${SRC}/out_mominxmf.f90
	${FC} ${FFLAGS} -c ${SRC}/out_mominvtk.f90
	${FC} ${FFLAGS} -c ${SRC}/out_mominrz.f90
	${FC} ${FFLAGS} -c ${SRC}/out_ffinkxky.f90
	${FC} ${FFLAGS} -c ${SRC}/out_ffinvm.f90
	${FC} ${FFLAGS} -c ${SRC}/out_ffinzv.f90
	${FC} ${FFLAGS} -c ${SRC}/out_ffinzvm_vtk.f90
	${FC} ${FFLAGS} -c ${SRC}/out_mominfreq.f90
	${FC} ${FFLAGS} -c ${SRC}/out_linfreq.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_correlate.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_realsp.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_bicoherence.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_zfshearing.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_zfdensity.f90
	${FC} ${FFLAGS} -c ${SRC}/out_fluidtotaltrans.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_fluidsubsptrans.f90
	${FC} ${FFLAGS} -c ${SRC}/out_fluiddetailtrans.f90
	${FC} ${FFLAGS} -c ${SRC}/diag_main.f90

	${FC} ${FFLAGS} *.o -o ${PROG} ${LIB}

	rm -f *.o *.mod

clean:
	rm -f *.o *.mod *.exe go.diag*.o* go.diag*.i*
