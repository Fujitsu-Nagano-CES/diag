MODULE diag_header
!-------------------------------------------------------------------------------
!
!     Header for data diagnostics
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  implicit none

  public

  integer, parameter :: DP = selected_real_kind(14)
  integer, parameter :: nhead = 4  ! Byte sizes of header and footer
                                   ! in sequential binary files
  integer, parameter :: nmom = 6   ! Number of output moments
  integer, parameter :: ntrn = 12  ! Number of output total transfer diagnostics
  integer, parameter :: ntri = 6   ! Number of output triad transfer diagnostics


!%%% DIAG parameters %%%
  integer, parameter :: snum = 1      ! begining of simulation runs
  integer, parameter :: enum = 1      ! end of simulation runs
!%%%%%%%%%%%%%%%%%%%%%%

!%%% GKV parameters %%%
  integer, parameter :: nxw = 10, nyw = 20
  integer, parameter :: nx = 0, global_ny = 12 ! 2/3 de-aliasing rule
  integer, parameter :: global_nz = 48, global_nv = 24, global_nm = 15

  integer, parameter :: nzb = 2, &  ! the number of ghost grids in z
                        nvb = 2     ! the number of ghost grids in v and m

!--------------------------------------
!  Data distribution for MPI
!--------------------------------------

  integer, parameter :: nprocw = 2, nprocz = 2, nprocv = 2, nprocm = 2, nprocs = 2
!%%%%%%%%%%%%%%%%%%%%%


  integer, parameter :: nxw_size = (2*nxw-1)/nprocw     ! local allocation size (0:nxw_size)
  integer, parameter :: ny       = global_ny / nprocw   ! local allocation size (0:ny)
  integer, parameter :: nz = global_nz / nprocz,          &
                        nv = global_nv / nprocv,          &
                        nm = (global_nm + 1) / nprocm - 1,&
                        ns = nprocs
  integer, parameter :: nproc = nprocw * nprocz * nprocv * nprocm * nprocs

! ---- y dimension -------
  integer :: ist_y(0:nprocw-1)               ! local start index of y
  integer :: iend_y(0:nprocw-1)              ! local end   index of y
  integer :: nsize_y(0:nprocw-1)             ! local size of y
  integer :: ist1_y(0:nprocw-1)              ! local start index of y for global start index 1 

  integer :: ist_y_g(0:nprocw-1)             ! global start index of y
  integer :: iend_y_g(0:nprocw-1)            ! global end   index of y

  real(kind=DP) :: tend                            ! end time
  real(kind=DP) :: dtout_fxv, dtout_ptn, dtout_eng ! time-spacing for output
  real(kind=DP) :: dtout_dtc                       ! time-spacing for dt control

  real(kind=DP),    parameter :: pi  = 3.141592653589793_DP, twopi = pi * 2._DP
  real(kind=DP),    parameter :: eps = 0.0000000001_DP
  complex(kind=DP), parameter :: ui  = ( 0._DP, 1._DP )


!--- unit number for file I/O ---
  integer, parameter ::  &
              inml = 5,  &
              olog = 10, &
              inbz = 11, &
              ivmc = 15, &
              ovmc = olog,&
         otextfile = 19, &
           onetcdf = 18, &
         offinkxky = 20, &
           offinvm = 21, &
           offinzv = 22, &
          offinzvm = 23, &
        omominkxky = 40, &
         omomintky = 42, &
          omominxy = 44, &
          omominxz = 46, &
           omominz = 48, &
         omominxyz = 50, &
          omominrz = 52, &
        omominfreq = 54, &
          olinfreq = 56, &
        ocorrelate = 58, &
      obicoherence = 60, &
          ophireal = 62, &
       ozfshearing = 64, &
           ozfsave = 66, &
        ozfdensity = 68, &
        otrninkxky = 70, &
        otriinkxky = 72, &
     oflttrninkxky = 74, &
     oflstrninkxky = 50000000, &
     ofldtrninkxky = 60000000


END MODULE diag_header
