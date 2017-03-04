program test_QUADS
  implicit none
  integer :: element_id
  integer :: itype 
  integer :: nint 
  integer :: iout 

  real    :: thick
  real    :: ym
  real    :: pr 
  real    :: xx(2,4)
  real    :: s(8,8)


  thick = 1.0 
  ym    = 2.83e7
  pr    = 0.3
  xx    = 0.0 
  s     = 0.0 

  element_id = 1
  itype      = 0
  nint       = 1
  iout       = 5

  call QUADS (element_id, itype, nint, thick, ym, pr, xx, s, iout)

end program