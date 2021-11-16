PROGRAM diag
!-------------------------------------------------------------------------------
!
!     Data diagnostics from binary output
!     Test for chgres_cnt
!
!-------------------------------------------------------------------------------
  use diag_header
  use diag_clock, only : clock_sta, clock_end, elt
  use diag_geom, only : geom_init
  use diag_rb, only : loop_cnt_sta, loop_cnt_end

  use out_netcdf, only : cntinnetcdf
  use diag_chgres_cnt

  implicit none

  real(kind=8) :: trange
  integer :: mx, gmy, giz, giv, gim, imom, is, &
             loop, flag, loop_sta, loop_end, loop_skp, it, mxt, myt, rankz

  integer :: stpn, nnx, ngy, ngz, ngv, ngm, nnpw, nnpz, nnpv, nnpm, nnps
  !character(len=*) :: outdir


!--- Initialize ---
    call initialize
    call geom_init( snum )  ! It can be re-called when namelist is changed in runs.
    trange=0.d0;mx=0;gmy=0;giz=0;giv=0;gim=0;imom=0;is=0
    loop=0;flag=0;loop_sta=0;loop_end=0;loop_skp=0;it=0;mxt=0;myt=0;rankz=nprocz/2+1
!------------------

                                                        call clock_sta(10)
!= example to write NetCDF data format =
!    write(*,*) "OUTPUT : NetCDF data"
!    call cntinnetcdf

    
!--- Change resolution ---
    nnx = nx
    ngy = global_ny
    ngz = global_nz
    ngv = global_nv
    ngm = global_nm
    nnpw= nprocw
    nnpz= nprocz
    nnpv= nprocv
    nnpm= nprocm
    nnps= nprocs
    
    call chgres_cnt_netcdf(nnx=nnx, ngy=ngy, ngz=ngz, ngv=ngv, ngm=ngm, &
        nnpw=nnpw, nnpz=nnpz, nnpv=nnpv, nnpm=nnpm, nnps=nnps)
!------------------
                                                        call clock_end(10)


!--- Finalization ---
    call finalize
    write( *, * ) ""
    write( *, * ) "########## elapsed time (sec) #########"
    write( *, * ) "# total = ", elt(10)
    write( *, * ) "#######################################"
!---------------------


END PROGRAM diag


SUBROUTINE initialize
!-------------------------------------------------------------------------------
!
!     Initialize diag
!                                                   (S. Maeyama, 9 Dec. 2016)
!
!-------------------------------------------------------------------------------
  use diag_clock, only : clock_init
  use diag_rb, only : rb_cnt_fileopen, rb_cnt_check, rb_cnt_setnloop
  implicit none

    call clock_init
    call rb_cnt_fileopen
    call rb_cnt_check
    call rb_cnt_setnloop

END SUBROUTINE initialize


SUBROUTINE finalize
!-------------------------------------------------------------------------------
!
!     Finalize diag
!                                                   (S. Maeyama, 9 Dec. 2016)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_cnt_fileclose
  implicit none

    call rb_cnt_fileclose

END SUBROUTINE finalize
