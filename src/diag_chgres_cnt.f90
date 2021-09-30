MODULE diag_chgres_cnt
!-------------------------------------------------------------------------------
!
!    diag_chgres_cnt: change the resolution and process division number for cnt
!                                                   (FUJITSU LTD, November 2021)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private renew_dir, check_params
  public  chgres_cnt_fortran, chgres_cnt_netcdf

CONTAINS

  SUBROUTINE renew_dir( dir )
    character(len=*), intent(in) :: dir
    character(len=256) :: comm
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

    
