PROGRAM diag
!-------------------------------------------------------------------------------
!
!     Data diagnostics from binary output
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  use diag_clock, only : clock_sta, clock_end, elt
  use diag_geom, only : geom_init
  use diag_rb, only : loop_phi_sta, loop_phi_end, &
                      loop_fxv_sta, loop_fxv_end, &
                      loop_trn_sta, loop_trn_end, &
                      loop_tri_sta, loop_tri_end, &
                      loop_cnt_sta, loop_cnt_end, &
                      num_triad_diag, triad_diag_mxt, triad_diag_myt

  use out_textfile, only : phiintext
  use out_netcdf, only : phiinnetcdf,  Alinnetcdf, mominnetcdf, &
                         fxvinnetcdf, cntinnetcdf, trninnetcdf, triinnetcdf
  use out_linfreq, only : linfreqintime, linfreqinkxky
  use out_mominfreq, only : phiinfreq
  use out_mominkxky, only : phiinkxky, Alinkxky, mominkxky
  use out_trninkxky, only : trninkxky
  use out_mominxy, only : phiinxy, Alinxy, mominxy,  &
                          phiinxy_parity, Alinxy_parity, mominxy_parity
  use out_mominz, only : phiinz, Alinz, mominz, phiinz_freq,                            &
                         phiinz_connect, Alinz_connect, mominz_connect, eneinz_connect, &
                         phiinz_parity, phiinz_parity_freq, &
                         Alinz_parity, Alinz_parity_freq
  use out_momintky, only : momintky_open, momintky_close,          &
                           wesintky, wemintky, engintky, menintky, &
                           wesintky_parity, wemintky_parity, engintky_parity, menintky_parity

  use out_ffinkxky, only : fkinkxky_fxv, fkinkxky_cnt
  use out_ffinvm, only : fkinvm_fxv, fkinvm_cnt, fluxinvm_fxv, fluxinvm_cnt
  use out_ffinzv, only : fkinzv
  use out_ffinzvm_vtk, only : fkinzvm_vtk, fkinzvm_connect_vtk

  use out_mominvtk, only : phiinvtk
  use out_mominxmf, only : mominxmf_coord, &
                           mominxmf_var_phi, mominxmf_header_phi, &
                           mominxmf_var_Al,  mominxmf_header_Al,  &
                           mominxmf_var_mom, mominxmf_header_mom
  use out_mominrz, only : phiinrz

  use out_triinkxky, only : triinkxky

  use out_fluidtotaltrans, only : fluidtotaltrans_isloop
  use out_fluiddetailtrans, only : fluiddetailtrans_mxmyisloop

!!!use out_mominxz, only : phiinxz, Alinxz, mominxz
!!!use out_triinkxky, only : triinkxky
!!!use out_mominavs, only : phiinavs
!!!use out_bicoherence, only : bicoh_phicheck, bicoh_freqcheck, bicoh_bicoherence
!!!use out_realsp, only : realsp_create
!!!use out_fluidsubsptrans, only : fluidsubsptrans_izloop, fluidsubsptrans_izloop_open,  &
!!!                                fluidsubsptrans_izloop_close, fluidsubsptrans_loop,  &
!!!                                fluidsubsptrans_loop_open, fluidsubsptrans_loop_close
!!!use out_zfshearing, only : zfs_open, zfs_close, zfshearing, zfshearing_timeaverage
!!!use out_zfdensity, only : zfd_open, zfd_close, zfdensity
  implicit none

  real(kind=8) :: trange
  integer :: mx, gmy, giz, giv, gim, imom, is, &
             loop, flag, loop_sta, loop_end, loop_skp, it, mxt, myt, rankz


!--- Initialize ---
    call initialize
    call geom_init( snum )  ! It can be re-called when namelist is changed in runs.
    trange=0.d0;mx=0;gmy=0;giz=0;giv=0;gim=0;imom=0;is=0
    loop=0;flag=0;loop_sta=0;loop_end=0;loop_skp=0;it=0;mxt=0;myt=0;rankz=nprocz/2+1
!------------------

                                                        call clock_sta(10)
!= example to write NetCDF data format =
    write(*,*) "OUTPUT : NetCDF data"
    call phiinnetcdf
    call  Alinnetcdf
    call mominnetcdf
    call fxvinnetcdf
    call cntinnetcdf
    call trninnetcdf
    do it = 0, num_triad_diag-1
      mxt = triad_diag_mxt(it)
      myt = triad_diag_myt(it)
      call triinnetcdf( mxt, myt )
    end do

!!= example to write linear frequency =
!    write(*,*) "OUTPUT : linfreqintime* linfreqinkxky*"
!    mx = 0
!    gmy = 3
!    loop_skp = 1
!    loop_sta = loop_phi_end(enum)/2
!    loop_end = loop_phi_end(enum)
!    call linfreqintime( mx, gmy, loop_sta, loop_end, loop_skp )
!    call linfreqinkxky( loop_end ) ! Evaluate linear dispersion omega(kx,ky)
!
!
!!= example to write moments in freq =
!    write(*,*) "OUTPUT : phiinfreq*"
!    mx = 0
!    gmy = 3
!    giz = 0
!    trange = 5._DP
!    call phiinfreq( mx, gmy, giz, trange )
!
!
!= example to write moments in kxky =
    write(*,*) "OUTPUT : phiinkxky* Alinkxky* mominkxky*"
    loop_skp = 10
    loop_sta = (floor(dble(loop_phi_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_phi_sta(snum)
    loop_end = loop_phi_end(enum)
    do loop = loop_sta, loop_end, loop_skp
      call phiinkxky( loop )
      call Alinkxky( loop )
      do is = 0, ns-1
        call mominkxky( is, loop )
      end do
    end do


!= example to write entropy balance diagnostics in kxky =
    write(*,*) "OUTPUT : trninkxky* "
    loop_skp = 10
    loop_sta = (floor(dble(loop_trn_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_trn_sta(snum)
    loop_end = loop_trn_end(enum)
    do loop = loop_sta, loop_end, loop_skp
      do is = 0, ns-1
        call trninkxky( is, loop )
      end do
    end do


!= example to write moments in xy =
    write(*,*) "OUTPUT : phiinxy* Alinxy* mominxy*"
    giz = 0
    loop_skp = 10
    loop_sta = (floor(dble(loop_phi_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_phi_sta(snum)
    loop_end = loop_phi_end(enum)
    do loop = loop_sta, loop_end, loop_skp
      call phiinxy( giz, loop )
      call Alinxy( giz, loop )
      do is = 0, ns-1
        call mominxy( giz, is, loop )
      end do
    end do
    !do loop = loop_sta, loop_end, loop_skp
    !  call phiinxy_parity( giz, loop )
    !  call Alinxy_parity( giz, loop )
    !  do is = 0, ns-1
    !    call mominxy_parity( giz, is, loop )
    !  end do
    !end do


!= example to write moments in zz =
    write(*,*) "OUTPUT : phiinz* Alinz* mominz*"
    mx = 0
    gmy = 6
    loop_skp = 10
    loop_sta = (floor(dble(loop_phi_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_phi_sta(snum)
    loop_end = loop_phi_end(enum)
    do loop = loop_sta, loop_end, loop_skp
!      call phiinz( mx, gmy, loop )
      call phiinz_connect( mx, gmy, loop )
!      call Alinz( mx, gmy, loop )
      call Alinz_connect( mx, gmy, loop )
      do is = 0, ns-1
!        call mominz( mx, gmy, is, loop )
        call mominz_connect( mx, gmy, is, loop )
      end do
!      call eneinz_connect( mx, gmy, loop )
    end do
!    call phiinz_freq( mx, gmy )
!    do loop = loop_sta, loop_end, loop_skp
!      call phiinz_parity( mx, gmy, loop )
!      call Alinz_parity( mx, gmy, loop )
!    end do
!    call phiinz_parity_freq( mx, gmy )
!    call Alinz_parity_freq( mx, gmy )


!!= example to write moments in tky =
!    write(*,*) "OUTPUT : momintky* "
!    loop_skp = 10
!    loop_sta = (floor(dble(loop_phi_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_phi_sta(snum)
!    loop_end = loop_phi_end(enum)
!    call momintky_open( "./data/wesintky_parity.dat" )
!    do loop = loop_sta, loop_end, loop_skp
!      call wesintky_parity( loop )
!    end do
!    call momintky_close
!    call momintky_open( "./data/wemintky_parity.dat" )
!    do loop = loop_sta, loop_end, loop_skp
!      call wemintky_parity( loop )
!    end do
!    call momintky_close
!    call momintky_open( "./data/engintky_parity.dat" )
!    do loop = loop_sta, loop_end, loop_skp
!      call engintky_parity( loop )
!    end do
!    call momintky_close
!    call momintky_open( "./data/menintky_parity.dat" )
!    do loop = loop_sta, loop_end, loop_skp
!      call menintky_parity( loop )
!    end do
!    call momintky_close
!
!
!!= example to write ff in kxky =
!    write(*,*) "OUTPUT : fkinkxky*"
!    != read from fxv/*fxv* =
!    rankz = nprocz/2
!    giv = 8
!    gim = 4
!    loop_skp = 10
!    loop_sta = (floor(dble(loop_fxv_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_fxv_sta(snum)
!    loop_end = loop_fxv_end(enum)
!    do loop = loop_sta, loop_end, loop_skp
!      do is = 0, ns-1
!        call fkinkxky_fxv( rankz, giv, gim, is, loop )
!      end do
!    end do
!    != read from cnt/*cnt* =
!    !giz = 0
!    !giv = 8
!    !gim = 4
!    !loop_skp = 10
!    !loop_sta = (floor(dble(loop_cnt_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_cnt_sta(snum)
!    !loop_end = loop_cnt_end(enum)
!    !do loop = loop_sta, loop_end, loop_skp
!    !  do is = 0, ns-1
!    !    call fkinkxky_cnt( giz, giv, gim, is, loop )
!    !  end do
!    !end do
!    
!
!!= example to write ff in vm =
!    write(*,*) "OUTPUT : fkinvm*"
!    != read from fxv/*fxv* =
!    mx = 0
!    gmy = 1
!    rankz = nprocz/2
!    loop_skp = 10
!    loop_sta = (floor(dble(loop_fxv_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_fxv_sta(snum)
!    loop_end = loop_fxv_end(enum)
!    do loop = loop_sta, loop_end, loop_skp
!      do is = 0, ns-1
!        call fkinvm_fxv( mx, gmy, rankz, is, loop )
!        call fluxinvm_fxv( rankz, is, loop )
!      end do
!    end do
!    != read from cnt/*cnt* =
!    !mx = 0
!    !gmy = 3
!    !giz = 0
!    !loop_skp = 10
!    !loop_sta = (floor(dble(loop_cnt_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_cnt_sta(snum)
!    !loop_end = loop_cnt_end(enum)
!    !do loop = loop_sta, loop_end, loop_skp
!    !  do is = 0, ns-1
!    !    call fkinvm_cnt( mx, gmy, giz, is, loop )
!    !    call fluxinvm_cnt( giz, is, loop )
!    !  end do
!    !end do
!
!
!!= example to write ff in zv =
!    write(*,*) "OUTPUT : fkinzv*"
!    != read from cnt/*cnt* =
!    mx = 0
!    gmy = 3
!    gim = 4
!    loop_skp = 10
!    loop_sta = (floor(dble(loop_cnt_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_cnt_sta(snum)
!    loop_end = loop_cnt_end(enum)
!    do loop = loop_sta, loop_end, loop_skp
!      do is = 0, ns-1
!        !do gim = 0, global_nm, (global_nm+1)/4
!          call fkinzv( mx, gmy, gim, is, loop )
!        !end do
!      end do
!    end do
!
!
!!= example to write ff in zvm =
!    write(*,*) "OUTPUT : fkinzvm*"
!    != read from cnt/*cnt* =
!    mx = 0
!    gmy = 1
!    loop_skp = 10
!    loop_sta = (floor(dble(loop_cnt_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_cnt_sta(snum)
!    loop_end = loop_cnt_end(enum)
!    do is = 0, ns-1
!      call fkinzvm_vtk( mx, gmy, is, loop_sta, loop_end, loop_skp )
!      !call fkinzvm_connect_vtk( mx, gmy, is, loop_sta, loop_end, loop_skp )
!    end do
!
!
!= example to write moments in xyz =
    write(*,*) "OUTPUT : phiinxmf* "
    loop_skp = max(1, (loop_phi_end(enum) - loop_phi_sta(snum))/100)
    loop_sta = (floor(dble(loop_phi_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_phi_sta(snum)
    loop_end = loop_phi_end(enum)
    !- XDFM (.xmf) file describes treatment of binary data, in extensible Markup Language (XML).
    call mominxmf_coord                                      ! output coordinate binary
    call mominxmf_var_phi( loop_sta, loop_end, loop_skp )    ! output phi binary
    call mominxmf_header_phi( loop_sta, loop_end, loop_skp ) ! output phiinxmf_header_*.xmf
    !call mominxmf_var_Al( loop_sta, loop_end, loop_skp )     ! output Al binary
    !call mominxmf_header_Al( loop_sta, loop_end, loop_skp )  ! output Alinxmf_header_*.xmf
    !do is = 0, ns-1
    !  do imom = 0, nmom-1
    !    call mominxmf_var_mom( imom, is, loop_sta, loop_end, loop_skp ) ! output mom binary
    !    call mominxmf_header_mom( imom, is, loop_sta, loop_end, loop_skp ) 
    !                                                         ! output mominxmf_header_*.xmf
    !  end do
    !end do
!!
!!    !- Output XML-based VTK format is also available -
!!    flag=1; call phiinvtk( flag, loop_sta, loop_end, loop_skp )  ! output coord&var(fluxtube)
!!    flag=3; call phiinvtk( flag, loop_sta, loop_end, loop_skp )  ! output coord&var(full torus)
!!    flag=5; call phiinvtk( flag, loop_sta, loop_end, loop_skp )  ! output coord&var(field aligned)
!!
!
!= example to write moments in xy =
    write(*,*) "OUTPUT : phiinrz* "
    loop_skp = 10
    loop_sta = (floor(dble(loop_phi_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_phi_sta(snum)
    loop_end = loop_phi_end(enum)
    do loop = loop_sta, loop_end, loop_skp
      call phiinrz( loop )
    end do


!!= example to write triad transfer diagnostics in kxky =
!    write(*,*) "OUTPUT : triinkxky* "
!    loop_skp = 10
!    loop_sta = (floor(dble(loop_tri_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_tri_sta(snum)
!    loop_end = loop_tri_end(enum)
!    do loop = loop_sta, loop_end, loop_skp
!      do is = 0, ns-1
!        do it = 0, num_triad_diag-1
!          mxt = triad_diag_mxt(it)
!          myt = triad_diag_myt(it)
!          call triinkxky( mxt, myt, is, loop )
!        end do
!      end do
!    end do
!
!
!!= example to write total triad transfer in fluid approximation =
!    write(*,*) "OUTPUT : fluidtotaltransinkxky* "
!    loop_skp = max(1, (loop_phi_end(enum) - loop_phi_sta(snum))/100)
!    loop_sta = (floor(dble(loop_phi_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_phi_sta(snum)
!    loop_end = loop_phi_end(enum)
!    do loop = loop_sta, loop_end, loop_skp
!      do is = 0, ns-1
!        call fluidtotaltrans_isloop( is, loop )
!      end do
!    end do
!
!
!!= example to write detailed triad transfer in fluid approximation =
!    write(*,*) "OUTPUT : fluiddetailtransinkxky* "
!    loop_skp = max(1, (loop_phi_end(enum) - loop_phi_sta(snum))/100)
!    loop_sta = (floor(dble(loop_phi_sta(snum)-1)/loop_skp)+1)*loop_skp ! loop_phi_sta(snum)
!    loop_end = loop_phi_end(enum)
!    do loop = loop_sta, loop_end, loop_skp
!      do is = 0, ns-1
!        call fluiddetailtrans_mxmyisloop(diag_mx=0, diag_my=2, is=is, loop=loop )
!!        call fluiddetailtrans_mxmyisloop(diag_mx=1, diag_my=0, is=is, loop=loop )
!      end do
!    end do




!!!
!!!  != example to write moments in xz =
!!!    gmy=3
!!!    do loop = loop_phi_sta(snum), loop_phi_end(enum), 1
!!!      call phiinxz( gmy, loop )
!!!      call Alinxz( gmy, loop )
!!!      do is = 0, ns-1
!!!        do imom = 0, nmom-1
!!!          call mominxz( gmy, imom, is, loop )
!!!        end do
!!!      end do
!!!    end do
!!!
!!!
!!!  != example to write auto bicoherence in frequency space =
!!!  ! call realsp_create( 0, 1999, 2002 )  ! create data for bicoherence analysis
!!!  ! call bicoh_phicheck( 2000 )
!!!  ! call bicoh_freqcheck( 0, 0, 10._DP )
!!!  ! call bicoh_bicoherence( 5._DP, 10._DP )
!!!
!!!
!!!  != example to write subspace transfer in fluid approximation =
!!!  ! giz = 0
!!!  ! loop_sta = loop_phi_sta(snum); loop_end = loop_phi_end(enum);
!!!  ! loop_skp = max(1, (loop_end - loop_sta)/100)
!!!  ! call fluidsubsptrans_loop_open
!!!  ! do loop = loop_sta, loop_end, loop_skp
!!!  !   call fluidsubsptrans_loop( loop )
!!!  ! end do
!!!  ! call fluidsubsptrans_loop_close
!!!    
!!!
!!!  != example to write zonal flow shearing rate =
!!!  ! call zfs_open
!!!  ! loop_sta = loop_phi_sta(snum); loop_end = loop_phi_end(enum);
!!!  ! loop_skp = max(1, (loop_end - loop_sta)/100)
!!!  ! do loop = loop_sta, loop_end, loop_skp
!!!  !   call zfshearing( loop )
!!!  ! end do
!!!  ! call zfs_close
!!!  ! call zfshearing_timeaverage( 5000, 9000, 100 )
!!!  ! call zfshearing_timeaverage( 880, 1280, 10 )
!!!
!!!
!!!  != example to write zonal fluctuations =
!!!  ! call zfd_open
!!!  ! loop_sta = loop_phi_sta(snum); loop_end = loop_phi_end(enum);
!!!  ! loop_skp = max(1, (loop_end - loop_sta)/100)
!!!  ! do loop = loop_sta, loop_end, loop_skp
!!!  !   call zfdensity( loop )
!!!  ! end do
!!!  ! call zfd_close


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
  use diag_rb, only : rb_phi_fileopen, rb_phi_check, rb_phi_setnloop,  &
                      rb_Al_fileopen,  rb_Al_check,  rb_Al_setnloop,   &
                      rb_mom_fileopen, rb_mom_check, rb_mom_setnloop,  &
                      rb_trn_fileopen, rb_trn_check, rb_trn_setnloop,  &
                      rb_tri_fileopen, rb_tri_check, rb_tri_setnloop,  &
                      rb_fxv_fileopen, rb_fxv_check, rb_fxv_setnloop,  &
                      rb_cnt_fileopen, rb_cnt_check, rb_cnt_setnloop
  use diag_fft, only : fft_pre
  implicit none

    call clock_init
    call rb_phi_fileopen
    call rb_phi_check
    call rb_phi_setnloop
    call rb_Al_fileopen
    call rb_Al_check
    call rb_Al_setnloop
    call rb_mom_fileopen
    call rb_mom_check
    call rb_mom_setnloop
    call rb_trn_fileopen
    call rb_trn_check
    call rb_trn_setnloop
    call rb_tri_fileopen
    call rb_tri_check
    call rb_tri_setnloop
    call rb_fxv_fileopen
    call rb_fxv_check
    call rb_fxv_setnloop
    call rb_cnt_fileopen
    call rb_cnt_check
    call rb_cnt_setnloop
    call fft_pre

END SUBROUTINE initialize


SUBROUTINE finalize
!-------------------------------------------------------------------------------
!
!     Finalize diag
!                                                   (S. Maeyama, 9 Dec. 2016)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_fileclose,  &
                      rb_Al_fileclose,   &
                      rb_mom_fileclose,  &
                      rb_trn_fileclose,  &
                      rb_tri_fileclose,  &
                      rb_fxv_fileclose,  &
                      rb_cnt_fileclose
  implicit none

    call rb_phi_fileclose
    call rb_Al_fileclose
    call rb_mom_fileclose
    call rb_trn_fileclose
    call rb_tri_fileclose
    call rb_fxv_fileclose
    call rb_cnt_fileclose

END SUBROUTINE finalize
