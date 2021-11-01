subroutine read_fort_cnt(path, nx, ny, nz, nv, nm, time, cnt)
  implicit none
  integer, parameter :: DP = selected_real_kind(14)
  character(*), intent(in) :: path
  integer, intent(in) :: nx, ny, nz, nv, nm
  real(kind=DP), intent(out) :: time
  complex(kind=DP), dimension(-nx:nx, 0:ny, -nz:nz-1, 1:2*nv, 0:nm), &
       intent(out) :: cnt

  open(unit=1, file=trim(path), status="old", action="read", form="unformatted")
  read(unit=1) time, cnt
  close(unit=1)
  return
end subroutine read_fort_cnt
