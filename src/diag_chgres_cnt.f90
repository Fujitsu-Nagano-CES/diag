MODULE diag_chgres_cnt
!-------------------------------------------------------------------------------
!
!    diag_chgres_cnt: change the resolution and process division number for cnt
!                                                   (FUJITSU LTD, November 2021)
!
!-------------------------------------------------------------------------------
  use diag_header
  use diag_rb, only : rb_cnt_mxmyimisloop
  use netcdf
  use out_netcdf, only : check_nf90err
  implicit none

  private renew_dir, check_params, open_fortran, close_fortran
  integer, parameter :: cntfos = 900000000, stpfos = 100000

  public  chgres_cnt_fortran, chgres_cnt_netcdf

CONTAINS

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
    if ( real(ngm+1)/real(nnpm) /= real((ngm+1)/nnpm) ) then
       write(*,*) "chgres_cnt: invalid ngm or nnpm has specified."
       stop
    end if
    return
  end SUBROUTINE check_params

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
                   open( unit=cntfos+stpfos*inum+ir,                         &
                        file=trim(outdir)//"/gkvp."//crank//".cnt."//cnum, &
                        status="new", action="write",                      &
                        form="unformatted", access="stream" )
                end do
             end do
          end do
       end do
    end do
  END SUBROUTINE open_fortran

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
    

  SUBROUTINE chgres_cnt_fortran( stpnum, nnx, ngy, ngz, ngv, ngm, &
       nnpw, nnpz, nnpv, nnpm, nnps, outdir )
!-------------------------------------------------------------------------------
!
!    change the resolution and process division number for cnt
!    and write into Fortran I/O files
!                                                   (FUJITSU LTD, November 2021)
!
!-------------------------------------------------------------------------------
    integer, optional, intent(in) :: stpnum, nnx, ngy, ngz, ngv, ngm, &
         nnpw, nnpz, nnpv, nnpm, nnps
    character(len=*), optional, intent(in) :: outdir

    integer :: n_nx, n_gy, n_gz, n_gv, n_gm, n_npw, n_npz, n_npv, n_npm, n_nps
    integer :: ips, ipm, ipv, ipz, ipw, ir, iunit
    complex(kind=DP), dimension(:,:,:), allocatable :: off, nff

    ! check stpnum
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

    ! allocate work for org cnt
    allocate( off(-nx:nx, 0:gy, -gz:gz-1) )

    ! allocate work for new cnt
    allocate( nff(-n_nx:n_nx, 0:n_gy, -n_gz:n_gz-1) )

    ! open fortran files
    call open_fortran(stpnum, n_npw, n_npz, n_npv, n_npm, n_nps, outdir)

    ! main loop (in new process division)
    do ips = 0, n_nps-1
       do ipm = 0, n_npm-1
          do ipv = 0, n_npv-1
             do ipz = 0, n_npz-1
                do ipw = 0, n_npw-1
                   ir = ipw + n_npw*ipz + n_npw*n_npz*ipv  &
                        + n_npw*n_npz*n_npv*ipm + n_npw*n_npz*n_npv*n_npm*ips
                   iunit = cntfos + stpfos*stpnum + ir
                   
