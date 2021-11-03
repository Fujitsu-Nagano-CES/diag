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
  use diag_geom, only : lz, vmax, mmax, dz, dv, dm
  use diag_interp
  use netcdf
  !use out_netcdf, only : check_nf90err
  implicit none

  private renew_dir, check_params, get_org_ivim
  integer, parameter :: cntfos = 900000000

  public  chgres_cnt_fortran !, chgres_cnt_netcdf

    
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
       do ipm = 0, n_npm-1
       do ipv = 0, n_npv-1
       do ipz = 0, n_npz-1
       do ipw = 0, n_npw-1
          nff(:, :, :, :, :) = (0.0, 0.0)

          ! setup rank and loop number string
          ir = ipw + n_npw*ipz + n_npw*n_npz*ipv &
               + n_npw*n_npz*n_npv*ipm + n_npw*n_npz*n_npv*n_npm*ips
          write(crank, fmt="(i6.6)") ir
          write(cnum, fmt="(i3.3)" ) loop

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
                call rb_cnt_ivimisloop(oiv(1), oim(1), ips, loop, woff)
                off(:, :, :, 1, 1) = woff
                ! (v, m) = (2, 1)
                call rb_cnt_ivimisloop(oiv(2), oim(1), ips, loop, woff)
                off(:, :, :, 2, 1) = woff
                ! (v, m) = (1, 2)
                call rb_cnt_ivimisloop(oiv(1), oim(2), ips, loop, woff)
                off(:, :, :, 1, 2) = woff
                ! (v, m) = (2, 2)
                call rb_cnt_ivimisloop(oiv(2), oim(2), ips, loop, woff)
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

    ! finalize interpolator
    call intp5d%finalize()

    ! dealocate array
    deallocate(nff)

    return
  end SUBROUTINE chgres_cnt_fortran

end MODULE diag_chgres_cnt
