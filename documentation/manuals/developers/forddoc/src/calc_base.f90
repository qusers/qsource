!------------------------------------------------------------------------------!
!  Q v.5.7 makefile                                                            !
!  Code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      !
!  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          !
!  Kajsa Ljunjberg, John Marelius, Martin Nervall                              !
!  Maintainers: Beat Amrein, Alexandre Barrozo, Paul Bauer, Mauricio Esguerra, !
!  Irek Szeler                                                                 !
!  latest update: july 13, 2015                                                !
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
