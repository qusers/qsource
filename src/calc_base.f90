!       (C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden
!       calc_base.f90
!       by John Marelius
!       shared data for trajectory calculation modules

module calc_base
use topo
implicit none
        
real(8), allocatable :: xin(:) ! All functions share this vector which contains a frame (J)
! real(8), allocatable :: xin2d(:,:) ! Schlitters formula needs all coordinates at the same time

end module calc_base
