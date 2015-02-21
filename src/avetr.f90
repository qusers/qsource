!------------------------------------------------------------------------------!
!  Q V.5.7                                                                     !
!  code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      !
!  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          !
!  Kajsa Ljunjberg, John Marelius, Martin Nervall                              !
!  Maintainers: Beat Amrein, Alexandre Barrozo, Paul Bauer,                    !
!  Mauricio Esguerra, Irek Szeler, Masoud Karemi                               !
!  latest update: February 9, 2015                                             !
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!  (c) 2004 Uppsala Molekylmekaniska HB, Uppsala, Sweden                       !
!  avetr.f90                                                                   !
!  by Martin Nervall                                                           !
!  date: March 2004                                                            !
!  average co-ordinates from qdyn trajectory files and write pdb-structure     !
!  tested to reproduce average structures from vmd                             !
!------------------------------------------------------------------------------!
module avetr
use prep
implicit none

integer, parameter             :: ave_pdb = 11
integer(4), private            :: ncoords, n_sets = 0
real(4), allocatable, private  :: x_in(:), x_sum(:), x2_sum(:)
real(8), private               :: rmsd
contains
!todo: *choose which frames, add more trajectories, divide x_sum every 100 steps

!******************************************************
!main subroutine
!******************************************************
subroutine avetr_calc
  integer :: i, allocation_status
  character(len=1) :: ans
  logical :: fin
  n_sets = 0
  call trajectory
  ncoords = trj_get_ncoords()
  allocate(x_in(ncoords), x_sum(ncoords), x2_sum(ncoords), &
                                stat=allocation_status)
  if (allocation_status .ne. 0) then
    write(*,*) 'out of memory!'
    return
  end if
  do while(trj_read_masked(x_in))  !add from first file
    call add_coordinates
  end do

  !add from multiple files
  fin = .false.
  do while(.not. fin)
    call get_string_arg(ans, '-----> add more frames? (y or n): ')
    if (ans .eq. 'y') then
          call trajectory
          do while(trj_read_masked(x_in))  !add from additional files
                call add_coordinates
          end do
        else 
          fin = .true.
        end if
  end do


  call average
  call write_average
  deallocate(x_in, x_sum, x2_sum, stat=allocation_status)
end subroutine avetr_calc

!******************************************************
!sum the coordinates and the sqared coordinates
!******************************************************
subroutine add_coordinates
        x_sum = x_sum + x_in
        x2_sum = x2_sum + x_in**2
        n_sets = n_sets +1
end subroutine add_coordinates

!******************************************************
!make average and rmsd
!******************************************************
subroutine average
        x_sum = x_sum / n_sets
        x2_sum = x2_sum / n_sets
        rmsd = sqrt(sum(x2_sum - x_sum**2)/ncoords)
end subroutine average

!******************************************************
!write average coords to pdb file.
!variables used from prep: mask
!variables used from topo: xtop
!******************************************************
subroutine write_average
        !assign masked coordinates to right atom in topology
        call mask_put(mask, xtop, x_sum)
    call writepdb
        write(*,'(a,f6.3,a)') 'root mean square co-ordinate deviation ', rmsd, ' a'
        x_sum = 0
        x2_sum = 0
end subroutine write_average

end module avetr
