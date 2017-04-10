!------------------------------------------------------------------------------!
!  Q version 5.7                                                               !
!  Code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      !
!  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          !
!  Petra Wennerstrom, Kajsa Ljunjberg, John Marelius, Martin Nervall,          !
!  Johan Sund, Ake Sandgren, Alexandre Barrozo, Masoud Kazemi, Paul Bauer,     !
!  Miha Purg, Irek Szeler                                                      !
!  latest update: March 29, 2017                                               !
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!  (C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden
!  calc_base.f90
!  by John Marelius
!  shared data for trajectory calculation modules
!------------------------------------------------------------------------------!
module calc_base
use topo
implicit none

real(8), allocatable :: xin(:) ! All functions share this vector which contains a frame (J)
! real(8), allocatable :: xin2d(:,:) ! Schlitters formula needs all coordinates at the same time

end module calc_base
