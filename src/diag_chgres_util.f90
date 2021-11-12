!-------------------------------------------------------------------------------
!
!    diag_interp: interpolation complex array
!                                                   (FUJITSU LTD, November 2021)
!
!-------------------------------------------------------------------------------
MODULE diag_interp
  use diag_header, only: DP
  implicit none
  private

  real(kind=DP), parameter :: zero = 0.0_DP, one = 1.0_DP

  !----------------------------------------------------------------------------
  ! interpolator class for complex 5d(x, y, z, v, m) data
  !----------------------------------------------------------------------------
  type, public :: interp_5d
     complex(kind=DP), dimension(:,:,:,:,:), pointer :: f
     real(kind=DP), dimension(:), allocatable :: x, y, z, v, m
     integer :: ilox = 1, iloy = 1, iloz = 1, ilov = 1, ilom = 1
     logical :: initialized = .false.
   contains
     procedure, public :: initialize => initialize_interp_5d
     procedure, public :: interpolate => interpolate_interp_5d
     procedure, public :: finalize => finalize_interp_5d
  end type interp_5d


CONTAINS

  !-------------------------------------------------------------------------
  ! returns the indices in `xl` that bound `x`, to use for interpolation.
  ! and set mflag as blow
  !   if            x < xl(1)   then ileft=1,   iright=2,    mflag=-1
  !   if   xl(i) <= x < xl(i+1) then ileft=i,   iright=i+1,  mflag=0
  !   if   xl(n) == x           then ileft=n-1, iright=n,    mflag=0
  !   if    xl(n) < x           then ileft=n-1, iright=n,    mflag=1
  !-------------------------------------------------------------------------
  pure subroutine dintrv(xl, x, ilo, ileft, iright, mflag)
    implicit none
    real(kind=DP), dimension(:), intent(in) :: xl
    real(kind=DP), intent(in) :: x
    integer, intent(inout) :: ilo
    integer, intent(out) :: ileft, iright, mflag
    integer :: ihi, istep, imid, n

    n = size(xl)
    if ( n == 1 ) then
       return
    end if

    ihi = ilo + 1
    if ( ihi >= n ) then
       if ( x >= xl(n) ) then
          if ( x == xl(n) ) then
             mflag = 0
          else
             mflag = 1
          end if
          ileft = n - 1
          iright = n
          return
       end if
       if ( n <= 1 ) then
          mflag = -1
          ileft = 1
          iright = 2
          return
       end if
       ilo = n - 1
       ihi = n
    endif

    if ( x >= xl(ihi) ) then
       istep = 1
       do
          ilo = ihi
          ihi = ilo + istep
          if ( ihi >= n ) then
             if ( x >= xl(n) ) then
                if ( x == xl(n) ) then
                   mflag = 0
                else
                   mflag = 1
                end if
                ileft = n-1
                iright = n
                return
             end if
             ihi = n
          else if ( x >= xl(ihi) ) then
             istep = istep*2
             cycle
          endif
          exit
       end do
    else
       if ( x >= xl(ilo) ) then
          mflag = 0
          ileft = ilo
          iright = ilo + 1
          return
       end if
       istep = 1
       do
          ihi = ilo
          ilo = ihi - istep
          if ( ilo <= 1 ) then
             ilo = 1
             if ( x < xl(1) ) then
                mflag = -1
                ileft = 1
                iright = 2
                return
             end if
          elseif ( x < xl(ilo) ) then
             istep = istep*2
             cycle
          endif
          exit
       end do
    endif

    do
       imid = (ilo + ihi) / 2
       if ( imid == ilo ) then
          mflag = 0
          ileft = ilo
          iright = ilo + 1
          return
       end if
       if ( x < xl(imid) ) then
          ihi = imid
       else
          ilo = imid
       endif
    end do
  end subroutine dintrv
  
  !-------------------------------------------------------------------------
  ! constructor for interp_5d class
  !-------------------------------------------------------------------------
  subroutine initialize_interp_5d(me_, nx, ny, nz, nv, nm)
    implicit none
    class(interp_5d), intent(inout) :: me_
    integer, intent(in) :: nx, ny, nz, nv, nm
    integer :: i

    call me_%finalize()

    if ( nx < 2 .or. ny < 2 .or. nz < 2 .or. nv < 2 .or. nm < 2 ) then
       return
    end if

    allocate(me_%x(nx))
    allocate(me_%y(ny))
    allocate(me_%z(nz))
    allocate(me_%v(nv))
    allocate(me_%m(nm))

    do i = 1, nx
       !me_%x(i) = real(i-1)
       me_%x(i) = real(-(nx-1)/2 + (i-1))
    end do
    do i = 1, ny
       me_%y(i) = real(i-1)
    end do
    do i = 1, nz
       me_%z(i) = real(i-1)
    end do
    do i = 1, nv
       me_%v(i) = real(i-1)
    end do
    do i = 1, nm
       me_%m(i) = real(i-1)
    end do

    me_%initialized = .true.
  end subroutine initialize_interp_5d

  !-------------------------------------------------------------------------
  ! interpolation complex 5d data
  !-------------------------------------------------------------------------
  subroutine interpolate_interp_5d(me_, x, y, z, v, m, f, istat)
    implicit none
    class(interp_5d), intent(inout) :: me_
    real(kind=DP), intent(in) :: x, y, z, v, m
    complex(kind=DP), intent(out) :: f
    integer, intent(out), optional :: istat

    integer, dimension(2) :: ix, iy, iz, iv, im
    real(kind=DP) :: p1, p2, p3, p4, p5
    real(kind=DP) :: q1, q2, q3, q4, q5
    integer :: mflag
    complex(kind=DP) :: &
         fx1111, fx2111, fx1211, fx2211, fx1121, fx2121, fx1221, fx2221, &
         fxy111, fxy211, fxy121, fxy221, fxyz11, fxyz21, fxyzv1, fx1112, &
         fx2112, fx1212, fx2212, fx1122, fx2122, fx1222, fx2222, fxy112, &
         fxy212, fxy122, fxy222, fxyz12, fxyz22, fxyzv2
    
    if ( me_%initialized .eqv. .false. .or. .not. associated(me_%f)) then
       f = zero
       if ( present(istat) ) istat = -1
       return
    end if
    
    call dintrv(me_%x, x, me_%ilox, ix(1), ix(2), mflag)
    call dintrv(me_%y, y, me_%iloy, iy(1), iy(2), mflag)
    call dintrv(me_%z, z, me_%iloz, iz(1), iz(2), mflag)
    call dintrv(me_%v, v, me_%ilov, iv(1), iv(2), mflag)
    call dintrv(me_%m, m, me_%ilom, im(1), im(2), mflag)

    q1 = (x - me_%x(ix(1))) / (me_%x(ix(2)) - me_%x(ix(1)))
    q2 = (y - me_%y(iy(1))) / (me_%y(iy(2)) - me_%y(iy(1)))
    q3 = (z - me_%z(iz(1))) / (me_%z(iz(2)) - me_%z(iz(1)))
    q4 = (v - me_%v(iv(1))) / (me_%v(iv(2)) - me_%v(iv(1)))
    q5 = (m - me_%m(im(1))) / (me_%m(im(2)) - me_%m(im(1)))
    p1 = one - q1
    p2 = one - q2
    p3 = one - q3
    p4 = one - q4
    p5 = one - q5

    fx1111 = p1*me_%f(ix(1),iy(1),iz(1),iv(1),im(1)) &
         +   q1*me_%f(ix(2),iy(1),iz(1),iv(1),im(1))
    fx2111 = p1*me_%f(ix(1),iy(2),iz(1),iv(1),im(1)) &
         +   q1*me_%f(ix(2),iy(2),iz(1),iv(1),im(1))
    fx1211 = p1*me_%f(ix(1),iy(1),iz(2),iv(1),im(1)) &
         +   q1*me_%f(ix(2),iy(1),iz(2),iv(1),im(1))
    fx2211 = p1*me_%f(ix(1),iy(2),iz(2),iv(1),im(1)) &
         +   q1*me_%f(ix(2),iy(2),iz(2),iv(1),im(1))
    fx1121 = p1*me_%f(ix(1),iy(1),iz(1),iv(2),im(1)) &
         +   q1*me_%f(ix(2),iy(1),iz(1),iv(2),im(1))
    fx2121 = p1*me_%f(ix(1),iy(2),iz(1),iv(2),im(1)) &
         +   q1*me_%f(ix(2),iy(2),iz(1),iv(2),im(1))
    fx1221 = p1*me_%f(ix(1),iy(1),iz(2),iv(2),im(1)) &
         +   q1*me_%f(ix(2),iy(1),iz(2),iv(2),im(1))
    fx2221 = p1*me_%f(ix(1),iy(2),iz(2),iv(2),im(1)) &
         +   q1*me_%f(ix(2),iy(2),iz(2),iv(2),im(1))
    fx1112 = p1*me_%f(ix(1),iy(1),iz(1),iv(1),im(2)) &
         +   q1*me_%f(ix(2),iy(1),iz(1),iv(1),im(2))
    fx2112 = p1*me_%f(ix(1),iy(2),iz(1),iv(1),im(2)) &
         +   q1*me_%f(ix(2),iy(2),iz(1),iv(1),im(2))
    fx1212 = p1*me_%f(ix(1),iy(1),iz(2),iv(1),im(2)) &
         +   q1*me_%f(ix(2),iy(1),iz(2),iv(1),im(2))
    fx2212 = p1*me_%f(ix(1),iy(2),iz(2),iv(1),im(2)) &
         +   q1*me_%f(ix(2),iy(2),iz(2),iv(1),im(2))
    fx1122 = p1*me_%f(ix(1),iy(1),iz(1),iv(2),im(2)) &
         +   q1*me_%f(ix(2),iy(1),iz(1),iv(2),im(2))
    fx2122 = p1*me_%f(ix(1),iy(2),iz(1),iv(2),im(2)) &
         +   q1*me_%f(ix(2),iy(2),iz(1),iv(2),im(2))
    fx1222 = p1*me_%f(ix(1),iy(1),iz(2),iv(2),im(2)) &
         +   q1*me_%f(ix(2),iy(1),iz(2),iv(2),im(2))
    fx2222 = p1*me_%f(ix(1),iy(2),iz(2),iv(2),im(2)) &
         +   q1*me_%f(ix(2),iy(2),iz(2),iv(2),im(2))

    fxy111 = p2*fx1111 + q2*fx2111
    fxy211 = p2*fx1211 + q2*fx2211
    fxy121 = p2*fx1121 + q2*fx2121
    fxy221 = p2*fx1221 + q2*fx2221
    fxy112 = p2*fx1112 + q2*fx2112
    fxy212 = p2*fx1212 + q2*fx2212
    fxy122 = p2*fx1122 + q2*fx2122
    fxy222 = p2*fx1222 + q2*fx2222

    fxyz11 = p3*fxy111 + q3*fxy211
    fxyz21 = p3*fxy121 + q3*fxy221
    fxyz12 = p3*fxy112 + q3*fxy212
    fxyz22 = p3*fxy122 + q3*fxy222

    fxyzv1 = p4*fxyz11 + q4*fxyz21
    fxyzv2 = p4*fxyz12 + q4*fxyz22

    f = p5*fxyzv1 + q5*fxyzv2
    
    if ( present(istat) ) istat = 0
    return
  end subroutine interpolate_interp_5d

  !-------------------------------------------------------------------------
  ! destructor for interp_5d class
  !-------------------------------------------------------------------------
  subroutine finalize_interp_5d(me_)
    implicit none
    class(interp_5d), intent(inout) :: me_

    if ( associated(me_%f) ) nullify(me_%f)
    if ( allocated(me_%x) ) deallocate(me_%x)
    if ( allocated(me_%y) ) deallocate(me_%y)
    if ( allocated(me_%z) ) deallocate(me_%z)
    if ( allocated(me_%v) ) deallocate(me_%v)
    if ( allocated(me_%m) ) deallocate(me_%m)
    me_%ilox = 1
    me_%iloy = 1
    me_%iloz = 1
    me_%ilov = 1
    me_%ilom = 1
    me_%initialized = .false.
  end subroutine finalize_interp_5d

end MODULE diag_interp


!-------------------------------------------------------------------------------
!
!    diag_cache_cnt: cache for cnt(x,y,z) array
!                                                   (FUJITSU LTD, November 2021)
!
!-------------------------------------------------------------------------------
MODULE diag_cache_cnt
  use diag_header
  use diag_rb, only : rb_cnt_ivimisloop
  implicit none

  !----------------------------------------------------------------------------
  ! cache buffer class for cnt(x, y, z) block indexed by (v, m)
  !----------------------------------------------------------------------------
  type, public :: cntbuff
     complex :: idx = (-1.0, -1.0) ! (v, m*i)
     complex(kind=DP) :: pd(2*nx +1, global_ny +1, 2*global_nz)
     integer :: age = 0
   contains
     procedure, public :: set => set_idx_pd
  end type cntbuff

  public initialize_cache, get_blk, finalize_cache
  integer :: ips = 0, stp = 0
  complex :: ei = (0.0, 1.0)
  integer :: num_blk = 0
  type(cntbuff), dimension(:), allocatable :: cnt_lst

CONTAINS

  subroutine set_idx_pd(this, v, m, s, t)
    class(cntbuff), intent(inout) :: this
    integer, intent(in) :: v, m, s, t
    this%idx = real(v) + real(m)*ei
    if ( v < 0 .or. m < 0 ) then
       this%age = 0
       return
    end if
    call rb_cnt_ivimisloop(v, m, s, t, this%pd)
    this%age = 1
    return
  end subroutine set_idx_pd

  !-------------------------------------------------------------------------
  ! initialize chache list
  !-------------------------------------------------------------------------
  SUBROUTINE initialize_cache( is, istp, nb, istatus )
    integer, intent(in) :: is, istp
    integer, optional, intent(in) :: nb
    integer, optional, intent(out) :: istatus

    call finalize_cache
    
    ! adjust number of blks
    if ( present(nb) ) then
       num_blk = nb
    else
       num_blk = global_nv *2
    end if
    if ( num_blk < 2 ) num_blk = 2

    ! store ips and stp
    ips = is
    stp = istp

    ! allocate buffer
    if ( present(istatus) ) then
       allocate( cnt_lst(num_blk), stat=istatus )
    else
       allocate( cnt_lst(num_blk) )
    end if
    return
  end SUBROUTINE initialize_cache

  !-------------------------------------------------------------------------
  ! get block from chache list
  !-------------------------------------------------------------------------
  SUBROUTINE get_blk( v, m, p )
    integer, intent(in) :: v, m
    complex(kind=DP), dimension(:,:,:), intent(out) :: p
    complex :: idx
    integer :: i, j, vacancy
    logical :: found

    if ( v < 1 .or. v > 2*global_nv .or. m < 0 .or. m > global_nm ) then
       p = (0.0, 0.0)
       return
    end if

    idx = real(v) + real(m)*ei
    vacancy = 0; found = .false.
    do i = 1, num_blk
       if ( cnt_lst(i)%idx == idx ) then ! cache hit
          p = cnt_lst(i)%pd
          cnt_lst(i)%age = 1
          found = .true.
       else if ( cnt_lst(i)%age > 0 ) then
          cnt_lst(i)%age = cnt_lst(i)%age + 1
       else if ( vacancy == 0 ) then
          vacancy = i
       end if
    end do

    ! found in cache
    if ( found .eqv. .true. ) then
       return
    end if

    ! not found in cache
    if ( vacancy > 0 ) then ! vacancy exists
       call cnt_lst(vacancy)%set(v, m, ips, stp)
       p = cnt_lst(vacancy)%pd
       return
    end if

    ! select the oldest one
    vacancy = 1
    do i = 2, num_blk
       if ( cnt_lst(vacancy)%age < cnt_lst(i)%age ) then
          vacancy = i
       end if
    end do
    call cnt_lst(vacancy)%set(v, m, ips, stp)
    p = cnt_lst(vacancy)%pd
    return
  end SUBROUTINE get_blk

  !-------------------------------------------------------------------------
  ! finalize chache list
  !-------------------------------------------------------------------------
  SUBROUTINE finalize_cache
    integer :: i
    if ( .not. allocated(cnt_lst) ) return
    deallocate( cnt_lst )
    num_blk = 0
    return
  end SUBROUTINE finalize_cache

end MODULE diag_cache_cnt
