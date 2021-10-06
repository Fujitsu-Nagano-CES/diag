!-------------------------------------------------------------------------------
!
!    diag_chgres_cnt: change the resolution and process division number for cnt
!                                                   (FUJITSU LTD, November 2021)
!
!-------------------------------------------------------------------------------
MODULE diag_chgres_cnt
  use diag_header
  use diag_rb, only : rb_cnt_gettime, rb_cnt_mxmyimisloop
  use netcdf
  use out_netcdf, only : check_nf90err
  implicit none

  private renew_dir, check_params, open_fortran, close_fortran
  integer, parameter :: cntfos = 900000000, stpfos = 100000

  public  chgres_cnt_fortran, chgres_cnt_netcdf

    
CONTAINS

  !-------------------------------------------------------------------------
  ! renew_dir: renewal output directory
  !-------------------------------------------------------------------------
  SUBROUTINE renew_dir( dir )
    character(len=*), intent(in) :: dir
    character(len=512) :: comm
    integer :: status

    write(comm, *) 'if [ -d ', trim(dir), ' ]; then rm -rf ', trim(dir), '; fi'
    call system(comm, status)
    if ( status /= 0 ) then
       write(*,*) "chgres_cnt: system exec failed: ", comm
       stop
    end if

    write(comm, *) 'mkdir -p ', trim(dir)
    call system(comm, status)
    if ( status /= 0 ) then
       write(*,*) "chgres_cnt: mkdir failed: ", trim(dir)
       stop
    end if

    return
  end SUBROUTINE renew_dir

  !-------------------------------------------------------------------------
  ! check_params: check parameters
  !-------------------------------------------------------------------------
  SUBROUTINE check_params(nnx, ngy, ngz, ngv, ngm, nnpw, nnpz, nnpv, nnpm, nnps)
    integer, intent(in) :: nnx, ngy, ngz, ngv, ngm, nnpw, nnpz, nnpv, nnpm, nnps

    if ( nnx < 1 .or. ngy < 1 .or. ngz < 1 .or. ngv < 1 .or. ngm < 1 .or. &
         nnpw < 1 .or. nnpz < 1 .or. nnpv < 1 .or. nnpm < 1 .or. nnps < 1 ) then
       write(*,*) "chgres_cnt: negative or zero value has specified."
       stop
    end if
    if ( real(ngy)/real(nnpw) /= real(ngy/nnpw) ) then
       write(*,*) "chgres_cnt: invalid ngy or nnpw has specified."
       stop
    end if
    if ( real(ngz)/real(nnpz) /= real(ngz/nnpz) ) then
       write(*,*) "chgres_cnt: invalid ngz or nnpz has specified."
       stop
    end if
    if ( real(ngv)/real(nnpv) /= real(ngv/nnpv) ) then
       write(*,*) "chgres_cnt: invalid ngv or nnpv has specified."
       stop
    end if
    if ( real(ngm+1)/real(nnpm) /= real((ngm+1)/nnpm) .or. &
         (ngm+1)/nnpm -1 < 1 ) then
       write(*,*) "chgres_cnt: invalid ngm or nnpm has specified."
       stop
    end if
    return
  end SUBROUTINE check_params

  !-------------------------------------------------------------------------
  ! open_fortran: open files for chgres_cnt_fortran
  !-------------------------------------------------------------------------
  SUBROUTINE open_fortran( stpnum, nnpw, nnpz, nnpv, nnpm, nnps, outdir )
    integer, intent(in) :: stpnum, nnpw, nnpz, nnpv, nnpm, nnps
    character(len=*), intent(in) :: outdir
    character :: crank(6), cnum(3)
    integer :: inum, ranks, rankm, rankv, rankz, rankw, ir

    ! renew 'outdir'
    call renew_dir( outdir )

    ! open files
    inum = stpnum
    write( cnum, fmt="(i3.3)" ) inum
    do ranks = 0, nnpcs-1
       do rankm = 0, nnpm-1
          do rankv = 0, nnpv-1
             do rankz = 0, nnpz-1
                do rankw = 0, nnpw-1
                   ir = rankw + nnpw*rankz + nnpw*nnpz*rankv &
                        + nnpw*nnpz*nnpv*rankm + nnpw*nnpz*nnpv*nnpm*ranks
                   write( crank, fmt="(i6.6)" ) ir
                   open( unit=cntfos+stpfos*inum+ir,                       &
                        file=trim(outdir)//"/gkvp."//crank//".cnt."//cnum, &
                        status="new", action="write",                      &
                        form="unformatted", access="stream" )
                end do
             end do
          end do
       end do
    end do
  END SUBROUTINE open_fortran

  !-------------------------------------------------------------------------
  ! close_fortran: close files for chgres_cnt_fortran
  !-------------------------------------------------------------------------
  SUBROUTINE close_fortran( stpnum, nnpw, nnpz, nnpv, nnpm, nnps )
    integer, intent(in) :: stpnum, nnpw, nnpz, nnpv, nnpm, nnps
    integer :: inum, ranks, rankm, rankv, rankz, rankw, ir

    inum = stpnum
    do ranks = 0, nnps-1
       do rankm = 0, nnpm-1
          do rankv = 0, nnpv-1
             do rankz = 0, nnpz-1
                do rankw = 0, nnpw-1
                   ir = rankw + nnpw*rankz + nnpw*nnpz*rankv  &
                        + nnpw*nnpz*nnpv*rankm + nnpw*nnpz*nnpv*nnpm*ranks
                   close( unit=cntfos+stpfos*inum+ir )
                end do
             end do
          end do
       end do
    end do
  END SUBROUTINE close_fortran
    

!-------------------------------------------------------------------------------
!
!    change the resolution and process division number for cnt
!    and write into Fortran I/O files
!                                                   (FUJITSU LTD, November 2021)
!
!-------------------------------------------------------------------------------
  SUBROUTINE chgres_cnt_fortran( stpn, nnx, ngy, ngz, ngv, ngm, &
       nnpw, nnpz, nnpv, nnpm, nnps, outdir )
    integer, optional, intent(in) :: stpn, nnx, ngy, ngz, ngv, ngm, &
         nnpw, nnpz, nnpv, nnpm, nnps
    character(len=*), optional, intent(in) :: outdir

    integer :: n_nx, n_gy, n_gz, n_gv, n_gm, n_npw, n_npz, n_npv, n_npm, n_nps
    integer :: n_ny, n_nz, n_nv, n_nm
    integer :: stpnum, ips, ipm, ipv, ipz, ipw
    ! buffer for write to file
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: nff
    ! buffer for rb_cnt_ivimisloop
    complex(kind=DP) :: ff(-nx:nx, 0:global_ny, -global_nz:global_nz-1)
    character(len=*), parameter :: default_odir = "./chgres_cnt"
    character(len=512) :: odir

    ! check stpnum
    stpnum = merge(stpn, enum, present(stpn))
    if (stpnum < snum .or. stpnum > enum) then
       write(*,*) "chgres_cnt_fortran: invalid stpnum specified(out of range)"
       stop
    end if

    ! check new resolution and process division number
    n_nx = merge(nnx, nx, present(nnx))
    n_gy = merge(ngy, gy, present(ngy))
    n_gz = merge(ngz, gz, present(ngz))
    n_gv = merge(ngv, gv, present(ngv))
    n_gm = merge(ngm, gm, present(ngm))
    n_npw = merge(nnpw, npw, present(nnpw))
    n_npz = merge(nnpz, npz, present(nnpz))
    n_npv = merge(nnpv, npv, present(nnpv))
    n_npm = merge(nnpm, npm, present(nnpm))
    n_nps = merge(nnps, nps, present(nnps))
    call check_params(n_nx,n_gy,n_gz,n_gv,n_gm, n_npw,n_npz,n_npv,n_npm,n_nps)
    n_ny = n_gy / n_npw
    n_nz = n_gz / n_npz
    n_nz = n_gv / n_npv
    n_nm = (n_gm + 1) / n_npm - 1

    ! open fortran files
    if ( present(outdir) ) then
       odir = outdir
    else
       odir = default_odir
    end if
    call open_fortran(stpnum, n_npw, n_npz, n_npv, n_npm, n_nps, odir)

    ! allocate work for new cnt
    allocate( nff(-n_nx:n_nx, 0:n_ny, -n_nz:n_nz-1, 1:2*n_nv, 0:n_nm) )

    ! main loop (in new process division)
    do ips = 0, n_nps-1
       do ipm = 0, n_npm-1
          do ipv = 0, n_npv-1
             do ipz = 0, n_npz-1
                do ipw = 0, n_npw-1
                   ! fill in nff
                   
