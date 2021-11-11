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
  use diag_geom, only : lz, vmax, mmax, dz, dv, dm, kxmin, kymin
  use diag_interp
  use diag_cache_cnt
  use netcdf
  implicit none

  private renew_dir, check_params, get_org_ivim, check_nf90err
  integer, parameter :: cntfos = 900000000

  public  chgres_cnt_fortran , chgres_cnt_netcdf

    
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

  !-------------------------------------------------------------------------
  ! get_org_ivim: get indices range of v and m in original mesh
  !-------------------------------------------------------------------------
  SUBROUTINE get_org_ivim(v, m, oiv, oim, vflag, mflag)
    real(kind=DP), intent(in) :: v, m
    integer, dimension(2), intent(out) :: oiv, oim
    integer, optional, intent(out) :: vflag, mflag
    integer :: i

    ! indices of v
    if ( v < -vmax ) then ! equivalent to (i == 1)
       oiv(1) = 1; oiv(2) = 2
       if ( present(vflag) ) vflag = -1
    else if ( v > vmax ) then
       oiv(2) = 2 * global_nv; oiv(1) = oiv(2) -1
       if ( present(vflag) ) vflag = 1
    else
       if ( present(vflag) ) vflag = 0
       do i = 2, 2*global_nv
          if ( v < -vmax + (i-1)*dv ) then
             oiv(1) = i -1; oiv(2) = i
             exit
          end if
       end do
    end if

    ! indices of m
    if ( m < 0 ) then ! equivalent to (i == 0)
       oim(1) = 0; oim(2) = 1
       if ( present(mflag) ) mflag = -1
    else if ( m > mmax ) then
       oim(2) = global_nm; oim(1) = oim(2) -1
       if ( present(mflag) ) mflag = 1
    else
       if ( present(mflag) ) mflag = 0
       do i = 1, global_nm
          if ( m < i*dm ) then
             oim(1) = i -1; oim(2) = i
             exit
          end if
       end do
    end if

    return
  end SUBROUTINE get_org_ivim

  !-------------------------------------------------------------------------
  ! check_nf90err: check error message of nf90 (copied from out_netcdf.f90)
  !-------------------------------------------------------------------------
  SUBROUTINE check_nf90err(werr, comment)
    integer(kind=4), intent(in) :: werr
    character(len=*), intent(in) :: comment
    if(werr /= nf90_noerr) then
       write(*,*) comment//" "//trim(nf90_strerror(werr))
       stop
    end if
  END SUBROUTINE check_nf90err


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
    integer :: igx, igy, igz, igv, igm, ir, oiz, oiv(2), oim(2)
    real(kind=DP) :: time, z0, v0, m0, z1, v1, m1, n_dm, n_dv, n_dz
    real(kind=DP) :: xx, yy, zz, vv, mm
    complex(kind=DP) :: f
    ! buffer for write to file
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: nff
    ! buffer for rb_cnt_ivimisloop (v=2, m=2)
    complex(kind=DP), target :: &
         off(2*nx +1, global_ny +1, 2*global_nz, 2, 2)
    !     off(-nx:nx, 0:global_ny, -global_nz:global_nz-1, 2, 2)
    complex(kind=DP) :: woff(2*nx +1, global_ny +1, 2*global_nz)
    !complex(kind=DP) :: woff(-nx:nx, 0:global_ny, -global_nz:global_nz-1)
    character(len=*), parameter :: default_odir = "./chgres_cnt"
    character(len=512) :: odir
    character(len=6) :: crank
    character(len=3) :: cnum
    ! 5d complex interpolator
    type(interp_5d) :: intp5d

    ! check stpnum
    stpnum = merge(stpn, enum, present(stpn))
    if (stpnum < snum .or. stpnum > enum) then
       write(*,*) "chgres_cnt_fortran: invalid stpnum specified(out of range)"
       stop
    end if

    ! check new resolution and process division number
    n_nx = merge(nnx, nx, present(nnx))
    n_gy = merge(ngy, global_ny, present(ngy))
    n_gz = merge(ngz, global_nz, present(ngz))
    n_gv = merge(ngv, global_nv, present(ngv))
    n_gm = merge(ngm, global_nm, present(ngm))
    n_npw = merge(nnpw, nprocw, present(nnpw))
    n_npz = merge(nnpz, nprocz, present(nnpz))
    n_npv = merge(nnpv, nprocv, present(nnpv))
    n_npm = merge(nnpm, nprocm, present(nnpm))
    n_nps = merge(nnps, nprocs, present(nnps))
    call check_params(n_nx,n_gy,n_gz,n_gv,n_gm, n_npw,n_npz,n_npv,n_npm,n_nps)
    n_ny = n_gy / n_npw
    n_nz = n_gz / n_npz
    n_nv = n_gv / n_npv
    n_nm = (n_gm + 1) / n_npm - 1

    ! prepare directory for fortran files
    if ( present(outdir) ) then
       odir = outdir
    else
       odir = default_odir
    end if
    call renew_dir( odir )

    ! allocate work for new cnt
    !allocate( nff(-n_nx:n_nx, 0:n_ny, -n_nz:n_nz-1, 1:2*n_nv, 0:n_nm) )
    allocate( nff(2*n_nx+1, n_ny+1, 2*n_nz, 2*n_nv, n_nm+1) )

    ! new delta (z, v, m)
    n_dz = lz / real(ngz, kind=DP)
    n_dv = 2._DP * vmax / real(2 * n_nv * n_npv -1, kind=DP)
    n_dm = mmax / real(n_npm * (n_nm+1) -1, kind=DP)

    ! setup interpolator with original mesh
    call intp5d%initialize(nx*2+1, global_ny+1, global_nz*2, 2, 2)
    do oiz = 0, 2*global_nz-1
       intp5d%z(oiz+1) = -lz + dz*oiz
    end do
    
    ! main loop (in new process division)
    do loop = loop_cnt_sta(stpnum), loop_cnt_end(stpnum)
       ! get time
       call rb_cnt_gettime(loop, time)
       
       do ips = 0, n_nps-1
          
       ! initialize buffer cache
       call initialize_cache(ips, loop, 2*n_nv *3)
          
       do ipm = 0, n_npm-1
       do ipv = 0, n_npv-1
       do ipz = 0, n_npz-1
       do ipw = 0, n_npw-1
          nff(:, :, :, :, :) = (0.0, 0.0)

          ! setup rank and loop number string
          ir = ipw + n_npw*ipz + n_npw*n_npz*ipv &
               + n_npw*n_npz*n_npv*ipm + n_npw*n_npz*n_npv*n_npm*ips
          write(crank, fmt="(i6.6)") ir
          write(cnum, fmt="(i3.3)" ) stpnum

          ! case of extends 's'
          if ( ips >= nprocs ) then
             ! open FortranI/O file
             open(unit=cntfos, &
                  file=trim(odir)//"/gkvp."//crank//".cnt."//cnum, &
                  form="unformatted")

             ! write into FortranI/O file
             rewind cntfos
             write(unit=cntfos) time, nff

             ! close the file
             call flush(cntfos)
             close(cntfos)
             cycle
          end if

          ! new index range in global
          igx0 = -n_nx; igx1 = n_nx
          igy0 = ipw*(n_ny+1); igy1 = min(igy0+n_ny, n_gy)
          igz0 = -n_gz + ipz*2*n_nz; igz1 = igz0 + 2*n_nz -1
          igv0 = 1 + ipv*2*n_nv; igv1 = igv0 + 2*n_nv -1
          igm0 = ipm*(n_nm+1); igm1 = igm0 + n_nm

          ! new coordinate range (z, v, m)
          z0 = igz0 * n_dz; z1 = igz1 * n_dz
          v0 = -vmax + (igv0-1)*n_dv; v1 = -vmax + (igv1-1)*n_dv
          m0 = igm0 * n_dm; m1 = igm1 * n_dm

          ! in process loop
          do igm = igm0, igm1
             mm = m0 + (igm-igm0)*n_dm
             do igv = igv0, igv1
                vv = v0 + (igv-igv0)*n_dv
                
                ! get original indices around (v, m)
                call get_org_ivim(vv, mm, oiv, oim)

                ! get off(:, :, :, 1:2, 1:2) around (v, m)
                ! (v, m) = (1, 1)
                call get_blk(oiv(1), oim(1), woff)
                off(:, :, :, 1, 1) = woff
                ! (v, m) = (2, 1)
                call get_blk(oiv(2), oim(1), woff)
                off(:, :, :, 2, 1) = woff
                ! (v, m) = (1, 2)
                call get_blk(oiv(1), oim(2), woff)
                off(:, :, :, 1, 2) = woff
                ! (v, m) = (2, 2)
                call get_blk(oiv(2), oim(2), woff)
                off(:, :, :, 2, 2) = woff

                ! setup interpolator
                intp5d%v(1) = -vmax + dv * (oiv(1)-1)
                intp5d%v(2) = -vmax + dv * (oiv(2)-1)
                intp5d%m(1) = dm * oim(1)
                intp5d%m(2) = dm * oim(2)
                intp5d%f => off

                ! interplation loop
                do igz = igz0, igz1
                   zz = z0 + (igz-igz0)*n_dz
                   do igy = igy0, igy1
                      ! case of extends 'y'
                      if ( igy > global_nz ) then
                         do igx = igx0, igx1
                            nff(igx-igx0+1, igy-igy0+1, igz-igz0+1, &
                                 igv-igv0+1, igm-igm0+1) = (0.0, 0.0)
                         end do
                         cycle
                      end if
                      yy = real(igy)
                      do igx = igx0, igx1
                         ! case of extends 'x'
                         if ( igx < -nx .or. igx > nx ) then
                            f = (0.0, 0.0)
                         else
                            xx = real(igx)
                            call intp5d%interpolate(xx, yy, zz, vv, mm, f)
                         end if
                         nff(igx-igx0+1, igy-igy0+1, igz-igz0+1, &
                              igv-igv0+1, igm-igm0+1) = f
                      end do
                   end do
                end do

             end do ! igv
          end do ! igm

          ! open FortranI/O file
          open(unit=cntfos, &
               file=trim(odir)//"/gkvp."//crank//".cnt."//cnum, &
               form="unformatted")

          ! write into FortranI/O file
          rewind cntfos
          write(unit=cntfos) time, nff

          ! close the file
          call flush(cntfos)
          close(cntfos)

       end do ! ipw
       end do ! ipz
       end do ! ipv
       end do ! ipm
       end do ! ips

    end do ! loop

    ! finalize
    call intp5d%finalize
    call finalize_cache

    ! dealocate array
    deallocate(nff)

    return
  end SUBROUTINE chgres_cnt_fortran

!-------------------------------------------------------------------------------
!
!    change the resolution and process division number for cnt
!    and write into netCDF file
!                                                   (FUJITSU LTD, November 2021)
!
!-------------------------------------------------------------------------------
  SUBROUTINE chgres_cnt_netcdf( stpn, nnx, ngy, ngz, ngv, ngm, &
       nnpw, nnpz, nnpv, nnpm, nnps, outdir )
    integer, optional, intent(in) :: stpn, nnx, ngy, ngz, ngv, ngm, &
         nnpw, nnpz, nnpv, nnpm, nnps
    character(len=*), optional, intent(in) :: outdir

    integer :: n_nx, n_gy, n_gz, n_gv, n_gm, n_npw, n_npz, n_npv, n_npm, n_nps
    integer :: n_ny, n_nz, n_nv, n_nm
    real(kind=DP), dimension(:), allocatable :: cxl, cyl, czl, cvl, cml
    integer :: stpnum, ips, ipm, ipv, ipz, ipw, loop
    integer :: igx0, igy0, igz0, igv0, igm0, igx1, igy1, igz1, igv1, igm1
    integer :: igx, igy, igz, igv, igm, oiz, oiv(2), oim(2)
    real(kind=DP) :: time, z0, v0, m0, z1, v1, m1, n_dm, n_dv, n_dz
    real(kind=DP) :: xx, yy, zz, vv, mm
    complex(kind=DP) :: f
    ! buffer for write to file
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: nff
    ! buffer for rb_cnt_ivimisloop (v=2, m=2)
    complex(kind=DP), target :: &
         off(2*nx +1, global_ny +1, 2*global_nz, 2, 2)
    !     off(-nx:nx, 0:global_ny, -global_nz:global_nz-1, 2, 2)
    complex(kind=DP) :: woff(2*nx +1, global_ny +1, 2*global_nz)
    !complex(kind=DP) :: woff(-nx:nx, 0:global_ny, -global_nz:global_nz-1)
    character(len=*), parameter :: default_odir = "./chgres_cnt"
    character(len=512) :: odir
    character(len=3) :: cnum
    integer(kind=4) :: ncid, dimids(1:7), ierr_nf90
    integer(kind=4) :: varid_kx, varid_ky, varid_zz, varid_vl, varid_mu, &
         varid_is, varid_tt, varid_recnt, varid_imcnt
    integer(kind=4) :: start_time(1), count_time(1), start_cnt(7)
    integer(kind=4), target :: count_cnt0(7), count_cnt1(7)
    integer(kind=4), dimension(:), pointer :: count_cnt
    ! 5d complex interpolator
    type(interp_5d) :: intp5d

    ! check stpnum
    stpnum = merge(stpn, enum, present(stpn))
    if (stpnum < snum .or. stpnum > enum) then
       write(*,*) "chgres_cnt_netcdf: invalid stpnum specified(out of range)"
       stop
    end if

    ! check new resolution and process division number
    n_nx = merge(nnx, nx, present(nnx))
    n_gy = merge(ngy, global_ny, present(ngy))
    n_gz = merge(ngz, global_nz, present(ngz))
    n_gv = merge(ngv, global_nv, present(ngv))
    n_gm = merge(ngm, global_nm, present(ngm))
    n_npw = merge(nnpw, nprocw, present(nnpw))
    n_npz = merge(nnpz, nprocz, present(nnpz))
    n_npv = merge(nnpv, nprocv, present(nnpv))
    n_npm = merge(nnpm, nprocm, present(nnpm))
    n_nps = merge(nnps, nprocs, present(nnps))
    call check_params(n_nx,n_gy,n_gz,n_gv,n_gm, n_npw,n_npz,n_npv,n_npm,n_nps)
    n_ny = n_gy / n_npw
    n_nz = n_gz / n_npz
    n_nv = n_gv / n_npv
    n_nm = (n_gm + 1) / n_npm - 1

    ! prepare directory for netCDF file
    if ( present(outdir) ) then
       odir = outdir
    else
       odir = default_odir
    end if
    call renew_dir( odir )

    ! allocate work for new cnt
    !allocate( nff(-n_nx:n_nx, 0:n_ny, -n_nz:n_nz-1, 1:2*n_nv, 0:n_nm) )
    allocate( nff(2*n_nx+1, n_ny+1, 2*n_nz, 2*n_nv, n_nm+1) )

    ! new delta (z, v, m)
    n_dz = lz / real(ngz, kind=DP)
    n_dv = 2._DP * vmax / real(2 * n_nv * n_npv -1, kind=DP)
    n_dm = mmax / real(n_npm * (n_nm+1) -1, kind=DP)

    ! setup new coordinates list
    allocate( cxl(2*n_nx+1) )
    do igx = -n_nx, n_nx
       cxl(igx + n_nx + 1) = kxmin * real(igx, kind=DP)
    end do
    allocate( cyl(n_gy+1) )
    do igy = 0, n_gy
       cyl(igy + 1) = kymin * real(igy, kind=DP)
    end do
    allocate( czl(2*n_gz) )
    do igz = -n_gz, n_gz - 1
       czl(igz + n_gz + 1) = n_dz * real(igz, kind=DP)
    end do
    allocate( cvl(2*n_gv) )
    do igv = 1, 2*n_gv
       cvl(igv) = -vmax + n_dv * real(igv-1, kind=DP)
    end do
    allocate( cml(n_gm+1) )
    do igm = 0, n_gm
       cml(igm + 1) = 0.5 * (n_dm * real(igm, kind=DP))**2
    end do

    ! setup interpolator with original mesh
    call intp5d%initialize(nx*2+1, global_ny+1, global_nz*2, 2, 2)
    do oiz = 0, 2*global_nz-1
       intp5d%z(oiz+1) = -lz + dz*oiz
    end do

    ! open netCDF file
    write(cnum, fmt="(i3.3)" ) stpnum
    ierr_nf90 = nf90_create(path=trim(odir)//"/gkvp.cnt."//cnum//".nc", &
         cmode=NF90_NETCDF4, ncid=ncid)      
    call check_nf90err(ierr_nf90, "nf90_create")

    ! define dimensions in file
    ierr_nf90 = nf90_def_dim(ncid=ncid, name="kx", &
         len=int(2*n_nx+1, kind=4), dimid=dimids(1))
    ierr_nf90 = nf90_def_dim(ncid=ncid, name="ky", &
         len=int(n_gy+1, kind=4), dimid=dimids(2))
    ierr_nf90 = nf90_def_dim(ncid=ncid, name="zz", &
         len=int(2*n_gz, kind=4), dimid=dimids(3))
    ierr_nf90 = nf90_def_dim(ncid=ncid, name="vl", &
         len=int(2*n_gv, kind=4), dimid=dimids(4))
    ierr_nf90 = nf90_def_dim(ncid=ncid, name="mu", &
         len=int(n_gm+1, kind=4), dimid=dimids(5))
    ierr_nf90 = nf90_def_dim(ncid=ncid, name="is", &
         len=int(n_nps, kind=4),  dimid=dimids(6))
    ierr_nf90 = nf90_def_dim(ncid=ncid, name="t",  &
         len=NF90_UNLIMITED, dimid=dimids(7))
    call check_nf90err(ierr_nf90, "nf90_def_dim")

    ! define variables in file
    ierr_nf90 = nf90_def_var(ncid=ncid, name="kx", xtype=NF90_DOUBLE, &
         dimids=dimids(1), varid=varid_kx)
    ierr_nf90 = nf90_def_var(ncid=ncid, name="ky", xtype=NF90_DOUBLE, &
         dimids=dimids(2), varid=varid_ky)
    ierr_nf90 = nf90_def_var(ncid=ncid, name="zz", xtype=NF90_DOUBLE, &
         dimids=dimids(3), varid=varid_zz)
    ierr_nf90 = nf90_def_var(ncid=ncid, name="vl", xtype=NF90_DOUBLE, &
         dimids=dimids(4), varid=varid_vl)
    ierr_nf90 = nf90_def_var(ncid=ncid, name="mu", xtype=NF90_DOUBLE, &
         dimids=dimids(5), varid=varid_mu)
    ierr_nf90 = nf90_def_var(ncid=ncid, name="is", xtype=NF90_INT, &
         dimids=dimids(6), varid=varid_is)
    ierr_nf90 = nf90_def_var(ncid=ncid, name="t",  xtype=NF90_DOUBLE, &
         dimids=dimids(7), varid=varid_tt)
    ierr_nf90 = nf90_def_var(ncid=ncid, name="recnt", xtype=NF90_DOUBLE, &
         dimids=dimids(1:7), varid=varid_recnt)
    ierr_nf90 = nf90_def_var(ncid=ncid, name="imcnt", xtype=NF90_DOUBLE, &
         dimids=dimids(1:7), varid=varid_imcnt)
    call check_nf90err(ierr_nf90, "nf90_def_var")

    ! end of definition of file
    ierr_nf90 = nf90_enddef(ncid=ncid)
    call check_nf90err(ierr_nf90, "nf90_enddef")

    ! write variables: static coordinates
    ierr_nf90 = nf90_put_var(ncid=ncid, varid=varid_kx, values=cxl)
    ierr_nf90 = nf90_put_var(ncid=ncid, varid=varid_ky, values=cyl)
    ierr_nf90 = nf90_put_var(ncid=ncid, varid=varid_zz, values=czl)
    ierr_nf90 = nf90_put_var(ncid=ncid, varid=varid_vl, values=cvl)
    ierr_nf90 = nf90_put_var(ncid=ncid, varid=varid_mu, values=cml)
    ierr_nf90 = nf90_put_var(ncid=ncid, varid=varid_is, &
         values=(/ (ips, ips=0,n_nps-1) /))
    call check_nf90err(ierr_nf90, "nf90_putvar")

    count_time(:) = 1
    count_cnt0(:) = int((/ 2*n_nx+1, n_ny+1, 2*n_nz, 2*n_nv, n_nm+1, 1, 1 /))
    count_cnt1 = count_cnt0
    count_cnt1(2) = n_gy + 1 - (n_ny+1)*(n_npw-1)

    ! main loop (in new process division)
    do loop = loop_cnt_sta(stpnum), loop_cnt_end(stpnum)
       ! get time
       call rb_cnt_gettime(loop, time)
       start_time(:) = 1 + loop - loop_cnt_sta(stpnum)
       ierr_nf90 = nf90_put_var(ncid=ncid, varid=varid_tt, values=(/time/), &
            start=start_time, count=count_time)
       call check_nf90err(ierr_nf90, "nf90_putvar")
       
       do ips = 0, n_nps-1
          
       ! initialize buffer cache
       call initialize_cache(ips, loop, 2*n_nv *3)

       do ipm = 0, n_npm-1
       do ipv = 0, n_npv-1
       do ipz = 0, n_npz-1
       do ipw = 0, n_npw-1
          nff(:, :, :, :, :) = (0.0, 0.0)
          start_cnt(:) = int((/ 1, ipw*(n_ny+1)+1, ipz*(2*n_nz)+1, &
               ipv*(2*n_nv)+1, ipm*(n_nm+1)+1, ips+1, &
               loop-loop_cnt_sta(stpnum)+1 /))
          if ( ipw == n_npw-1 ) then
             count_cnt => count_cnt1
          else
             count_cnt => count_cnt0
          end if

          ! case of extends 's'
          if ( ips >= nprocs ) then
             ierr_nf90 = nf90_put_var(ncid=ncid, varid=varid_recnt, &
                  values=dble(nff), start=start_cnt, count=count_cnt)
             ierr_nf90 = nf90_put_var(ncid=ncid, varid=varid_imcnt, &
                  values=aimag(nff), start=start_cnt, count=count_cnt)
             call check_nf90err(ierr_nf90, "nf90_putvar")
             cycle
          end if

          ! new index range in global
          igx0 = -n_nx; igx1 = n_nx
          igy0 = ipw*(n_ny+1); igy1 = min(igy0+n_ny, n_gy)
          igz0 = -n_gz + ipz*2*n_nz; igz1 = igz0 + 2*n_nz -1
          igv0 = 1 + ipv*2*n_nv; igv1 = igv0 + 2*n_nv -1
          igm0 = ipm*(n_nm+1); igm1 = igm0 + n_nm

          ! new coordinate range (z, v, m)
          z0 = igz0 * n_dz; z1 = igz1 * n_dz
          v0 = -vmax + (igv0-1)*n_dv; v1 = -vmax + (igv1-1)*n_dv
          m0 = igm0 * n_dm; m1 = igm1 * n_dm

          ! in process loop
          do igm = igm0, igm1
             mm = m0 + (igm-igm0)*n_dm
             do igv = igv0, igv1
                vv = v0 + (igv-igv0)*n_dv
                
                ! get original indices around (v, m)
                call get_org_ivim(vv, mm, oiv, oim)

                ! get off(:, :, :, 1:2, 1:2) around (v, m)
                ! (v, m) = (1, 1)
                call get_blk(oiv(1), oim(1), woff)
                off(:, :, :, 1, 1) = woff
                ! (v, m) = (2, 1)
                call get_blk(oiv(2), oim(1), woff)
                off(:, :, :, 2, 1) = woff
                ! (v, m) = (1, 2)
                call get_blk(oiv(1), oim(2), woff)
                off(:, :, :, 1, 2) = woff
                ! (v, m) = (2, 2)
                call get_blk(oiv(2), oim(2), woff)
                off(:, :, :, 2, 2) = woff

                ! setup interpolator
                intp5d%v(1) = -vmax + dv * (oiv(1)-1)
                intp5d%v(2) = -vmax + dv * (oiv(2)-1)
                intp5d%m(1) = dm * oim(1)
                intp5d%m(2) = dm * oim(2)
                intp5d%f => off

                ! interplation loop
                do igz = igz0, igz1
                   zz = z0 + (igz-igz0)*n_dz
                   do igy = igy0, igy1
                      ! case of extends 'y'
                      if ( igy > global_nz ) then
                         do igx = igx0, igx1
                            nff(igx-igx0+1, igy-igy0+1, igz-igz0+1, &
                                 igv-igv0+1, igm-igm0+1) = (0.0, 0.0)
                         end do
                         cycle
                      end if
                      yy = real(igy)
                      do igx = igx0, igx1
                         ! case of extends 'x'
                         if ( igx < -nx .or. igx > nx ) then
                            f = (0.0, 0.0)
                         else
                            xx = real(igx)
                            call intp5d%interpolate(xx, yy, zz, vv, mm, f)
                         end if
                         nff(igx-igx0+1, igy-igy0+1, igz-igz0+1, &
                              igv-igv0+1, igm-igm0+1) = f
                      end do
                   end do
                end do

             end do ! igv
          end do ! igm

          ! set values onto netCDF
          ierr_nf90 = nf90_put_var(ncid=ncid, varid=varid_recnt, &
               values=dble(nff(:,1:count_cnt(2),:,:,:)), &
               start=start_cnt, count=count_cnt)
          ierr_nf90 = nf90_put_var(ncid=ncid, varid=varid_imcnt, &
               values=aimag(nff(:,1:count_cnt(2),:,:,:)), &
               start=start_cnt, count=count_cnt)
          call check_nf90err(ierr_nf90, "nf90_putvar")

       end do ! ipw
       end do ! ipz
       end do ! ipv
       end do ! ipm
       end do ! ips

    end do ! loop

    ! close netCDF file
    ierr_nf90 = nf90_close(ncid)
    call check_nf90err(ierr_nf90, "nf90_close")
    
    ! finalize
    call intp5d%finalize()
    call finalize_cache

    ! dealocate arrays
    deallocate(nff)
    deallocate(cxl, cyl, czl, cvl, cml)

    return
  end SUBROUTINE chgres_cnt_netcdf
  
end MODULE diag_chgres_cnt
