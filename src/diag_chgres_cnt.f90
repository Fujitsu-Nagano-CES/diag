!-------------------------------------------------------------------------------
!
!    diag_chgres_cnt: change the resolution and process division number for cnt
!                                                   (FUJITSU LTD, November 2021)
!
!-------------------------------------------------------------------------------
MODULE diag_chgres_cnt
  use diag_header
  use diag_rb, only : rb_cnt_gettime, rb_cnt_ivimisloop, &
       loop_cnt_sta, loop_cnt_end
  use diag_geom, only : lz, vmax, mmax
  use netcdf
  use out_netcdf, only : check_nf90err
  implicit none

  private renew_dir, check_params, open_fortfiles, close_fortfiles
  integer, parameter :: cntfos = 900000000

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
    integer :: wny

    if ( nnx < 1 .or. ngy < 1 .or. ngz < 1 .or. ngv < 1 .or. ngm < 1 .or. &
         nnpw < 1 .or. nnpz < 1 .or. nnpv < 1 .or. nnpm < 1 .or. nnps < 1 ) then
       write(*,*) "chgres_cnt: negative or zero value has specified."
       stop
    end if
    wny = ngy / nnpw
    if ( wny < 1 .or. ((wny+1)*nnpw - ngy -1)/(wny+1) > 1 ) then
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

  SUBROUTINE get_org_ivim(v, m, oiv, oim)
    real(kind=DP), intent(in) :: v, m
    integer, dimension(2), intent(out) :: oiv, oim

  end SUBROUTINE get_org_ivim

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
    integer :: stpnum, ips, ipm, ipv, ipz, ipw, loop
    integer :: igx0, igy0, igz0, igv0, igm0, igx1, igy1, igz1, igv1, igm1
    real(kind=DP) :: time, z0, v0, m0, z1, v1, m1
    ! buffer for write to file
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: nff
    ! buffer for rb_cnt_ivimisloop (x2 x2)
    complex(kind=DP) :: off(-nx:nx, 0:global_ny, -global_nz:global_nz-1, 2, 2)
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

    ! prepare directory for fortran files
    if ( present(outdir) ) then
       odir = outdir
    else
       odir = default_odir
    end if
    call renew_dir( odir )

    ! allocate work for new cnt
    allocate( nff(-n_nx:n_nx, 0:n_ny, -n_nz:n_nz-1, 1:2*n_nv, 0:n_nm) )

    ! new delta (z, v, m)
    n_dz = lz / real(ngz, kind=DP)
    n_dv = 2._DP * vmax / real(2 * n_nv * n_npv -1, kind=DP)
    n_dm = mmax / real(n_npm * (n_nm+1) -1, kind=DP)
    
    ! main loop (in new process division)
    do loop = loop_cnt_sta(stpnum), loop_cnt_end(stpnum)
       ! get time
       call rb_cnt_gettime(loop, time)
       
       do ips = 0, n_nps-1
       do ipm = 0, n_npm-1
       do ipv = 0, n_npv-1
       do ipz = 0, n_npz-1
       do ipw = 0, n_npw-1
          ! new index range in global
          igx0 = -n_nx; igx1 = n_nx
          igy0 = ipw*(n_ny+1); igy1 = min(igy+n_ny, n_gy)
          igz0 = -n_gz + ipz*2*n_nz; igz1 = igz0 + 2*n_nz -1
          igv0 = 1 + ipv*2*n_nv; igv1 = igv0 + 2*n_nv -1
          igm0 = ipm*(n_nm+1); igm1 = igm0 + n_nm

          ! new coordinate range (z, v, m)
          z0 = igz0 * n_dz; z1 = igz1 * n_dz
          v0 = -vmax + (igv0-1)*n_dv; v1 = -vmax + (igv1-1)*n_dv
          m0 = igm0 * n_dm; m1 = igm1 * n_dm

          ! in process loop
          do igm = igm0, igm1
             m = m0 + igm*n_dm
             do igv = igv0, igv1
                v = v0 + igv*n_dv
                
                ! get off(:, :, :, 1:2, 1:2) around (v, m)
                
