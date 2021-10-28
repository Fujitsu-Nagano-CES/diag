MODULE diag_fft
!-------------------------------------------------------------------------------
!
!    FFT module using fftw3
!                                                   (S. Maeyama, 12 Dec. 2012)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  include "fftw3.f"

  private

  integer(kind=8), save :: plan_backward_xy, plan_forward_xy, plan_backward_x,  &
                            plan_time2freq, plan_time2freq_r2c

  public   fft_pre, fft_backward_xy, fft_forward_xy, fft_backward_x1,  &
           fft_time2freq_init, fft_time2freq,         &
           fft_time2freq_r2c_init, fft_time2freq_r2c


 CONTAINS


SUBROUTINE fft_pre
!-------------------------------------------------------------------------------
!
!    Initialization of FFT
!                                                   (S. Maeyama, 12 Dec. 2012)
!
!-------------------------------------------------------------------------------
  complex(kind=DP), dimension(0:nxw,0:2*nyw-1) :: wwkk
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: wwxy

  complex(kind=DP), dimension(0:2*nxw-1) :: wwkx, wwxx

!  integer :: ierr, omp_get_max_threads

!%%% setting for fftw3 %%%
    call dfftw_plan_dft_c2r_2d( plan_backward_xy,  &
                                2*nxw, 2*nyw,      &
                                !int(2*nxw, kind=4), int(2*nyw, kind=4),      &
                                wwkk, wwxy,        &
                                !FFTW_ESTIMATE )
                                FFTW_MEASURE )
    call dfftw_plan_dft_r2c_2d( plan_forward_xy,  &
                                2*nxw, 2*nyw,     &
                                !int(2*nxw, kind=4), int(2*nyw, kind=4),     &
                                wwxy, wwkk,       &
                                !FFTW_ESTIMATE )
                                FFTW_MEASURE )
!%%% setting for fftw3_omp %%%
!    call dfftw_init_threads( ierr )
!    call dfftw_plan_with_nthreads( omp_get_max_threads() )
!    call dfftw_plan_dft_c2r_2d( plan_backward_xy,  &
!                                2*nxw, 2*nyw,      &
!                                wwkk, wwxy,        &
!                                FFTW_ESTIMATE )
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    call dfftw_plan_with_nthreads( 1 )
    call dfftw_plan_dft_1d( plan_backward_x,   &
                            2*nxw,             &
                            !int(2*nxw, kind=4),             &
                            wwkx, wwxx,        &
                            FFTW_BACKWARD,     &
                            FFTW_ESTIMATE )

END SUBROUTINE fft_pre


SUBROUTINE fft_backward_xy ( gww, wwxy )
!-------------------------------------------------------------------------------
!
!    Execution of backward 2D FFT: (kx,ky) -> (x,y)
!                                                   (S. Maeyama, 12 Dec. 2012)
!
!-------------------------------------------------------------------------------
  complex(kind=DP), intent(in), &
    dimension(-nx:nx,0:global_ny)              :: gww
  real(kind=DP), intent(out), &
    dimension(0:2*nxw-1,0:2*nyw-1)             :: wwxy

  complex(kind=DP), dimension(0:nxw,0:2*nyw-1) :: wwkk
  integer :: mx, my

    wwkk(:,:) = ( 0._DP, 0._DP )
    do my = 0, global_ny
      do mx = 0, nx
        wwkk(mx,my) = gww(mx,my)
      end do
    end do
    do my = 1, global_ny
      do mx = -nx, -1
        wwkk(-mx,2*nyw-my) = conjg( gww(mx,my) )
      end do
    end do
    mx = 0
      do my = 1, global_ny
        wwkk(mx,2*nyw-my) = conjg( wwkk(mx,my) )
      end do

    call dfftw_execute_dft_c2r( plan_backward_xy, wwkk, wwxy )


END SUBROUTINE fft_backward_xy


SUBROUTINE fft_forward_xy ( wwxy, ww )
!--------------------------------------
!  Execution of FFT
  
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1), intent(in)    :: wwxy
  complex(kind=DP), dimension(-nx:nx,0:global_ny), intent(out) :: ww

  complex(kind=DP), dimension(0:nxw,0:2*nyw-1) :: wwkk
  real(kind=DP) :: ceff_norm
  integer :: mx, my
  
    call dfftw_execute_dft_r2c( plan_forward_xy, wwxy, wwkk )

    ceff_norm = 1._DP / real(2*nxw * 2*nyw, kind=DP)

    do my = 0, global_ny
      do mx = 0, nx
        ww(mx,my) = ceff_norm * wwkk(mx,my)
      end do
    end do
    do my = 1, global_ny
      do mx = -nx, -1
        ww(mx,my) = ceff_norm * conjg( wwkk(-mx,2*nyw-my) )
      end do
    end do
    my = 0
      do mx = 1, nx
        ww(-mx,-my) = conjg( ww(mx,my) )
      end do


END SUBROUTINE fft_forward_xy


!SUBROUTINE fft_backward_x ( gww, wwxxky )
!!-------------------------------------------------------------------------------
!!
!!    Execution of backward 1D FFT: (kx,ky) -> (x,ky)
!!                                                   (S. Maeyama, 12 Dec. 2012)
!!
!!-------------------------------------------------------------------------------
!  complex(kind=DP), intent(in), &
!    dimension(-nx:nx,0:global_ny)              :: gww
!  complex(kind=DP), intent(out), &
!    dimension(0:2*nxw-1,0:global_ny)           :: wwxxky
!
!  complex(kind=DP), dimension(0:2*nxw-1) :: wwkx
!  integer :: mx, my
!
!    do my = 0, global_ny
!
!      do mx = 0, nx
!        wwkx(mx) = gww(mx,my)
!      end do
!      do mx = -nx, -1
!        wwkx(2*nxw+mx) = gww(mx,my)
!      end do
!      do mx = nx+1, 2*nxw-nx-1
!        wwkx(mx) = ( 0._DP, 0._DP )
!      end do
!
!      call dfftw_execute_dft( plan_backward_x, wwkx, wwxxky(:,my) )
!
!    end do
!
!END SUBROUTINE fft_backward_x


SUBROUTINE fft_backward_x1 ( gww, wwxx )
!-------------------------------------------------------------------------------
!
!    Execution of backward 1D FFT: (kx) -> (x)
!                                                   (S. Maeyama, 12 Dec. 2012)
!
!-------------------------------------------------------------------------------
  complex(kind=DP), intent(in), &
    dimension(-nx:nx)              :: gww
  complex(kind=DP), intent(out), &
    dimension(0:2*nxw-1)           :: wwxx

  complex(kind=DP), dimension(0:2*nxw-1) :: wwkx
  integer :: mx

    do mx = 0, nx
      wwkx(mx) = gww(mx)
    end do
    do mx = -nx, -1
      wwkx(2*nxw+mx) = gww(mx)
    end do
    do mx = nx+1, 2*nxw-nx-1
      wwkx(mx) = ( 0._DP, 0._DP )
    end do

    call dfftw_execute_dft( plan_backward_x, wwkx, wwxx )

END SUBROUTINE fft_backward_x1


SUBROUTINE fft_time2freq_init ( nsample )
!-------------------------------------------------------------------------------
!
!    Execution of 1D FFT: sample(time) -> spectrum(freq)
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: nsample
  complex(kind=DP), dimension(0:nsample-1) :: sample, spectrum

!    call dfftw_plan_with_nthreads( 1 )
    call dfftw_plan_dft_1d( plan_time2freq,    &
                            nsample,           &
                            sample, spectrum,  &
                            FFTW_FORWARD,      &
                            FFTW_ESTIMATE )

END SUBROUTINE fft_time2freq_init


SUBROUTINE fft_time2freq ( nsample, sample, spectrum )
!-------------------------------------------------------------------------------
!
!    Execution of 1D FFT: sample(time) -> spectrum(freq)
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: nsample
  complex(kind=DP), intent(in), &
    dimension(0:nsample-1)           :: sample
  complex(kind=DP), intent(out), &
    dimension(0:nsample-1)           :: spectrum

    call dfftw_execute_dft( plan_time2freq, sample, spectrum )
    spectrum(:) = spectrum(:) / real(nsample, kind=DP)

END SUBROUTINE fft_time2freq


SUBROUTINE fft_time2freq_r2c_init ( nsample )
!-------------------------------------------------------------------------------
!
!    Execution of 1D FFT: Real[sample(time)] -> Complex[spectrum(freq)]
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: nsample
  real(kind=DP), dimension(0:nsample-1) :: sample
  complex(kind=DP), dimension(0:nsample/2) :: spectrum

!    call dfftw_plan_with_nthreads( 1 )
    call dfftw_plan_dft_r2c_1d( plan_time2freq_r2c,    &
                                int(nsample,kind=4),               &
                                sample, spectrum,      &
                                FFTW_ESTIMATE )

END SUBROUTINE fft_time2freq_r2c_init


SUBROUTINE fft_time2freq_r2c ( nsample, sample, spectrum )
!-------------------------------------------------------------------------------
!
!    Execution of 1D FFT: Real[sample(time)] -> Complex[spectrum(freq)]
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: nsample
  real(kind=DP), intent(in), &
    dimension(0:nsample-1)           :: sample
  complex(kind=DP), intent(out), &
    dimension(0:nsample/2)           :: spectrum

    call dfftw_execute_dft_r2c( plan_time2freq_r2c, sample, spectrum )
    spectrum(:) = spectrum(:) / real(nsample, kind=DP)

END SUBROUTINE fft_time2freq_r2c


END MODULE diag_fft
