!	(c) 2000 uppsala molekylmekaniska hb, uppsala, sweden
!	trj.f90
!	by john marelius
!	q trajectory data, access and dcd format i/o

module trj

use atom_mask
use misc
implicit none

	character(*), private, parameter	:: module_version = '5.01'
	character(*), private, parameter	:: module_date = '2003-06-02'

	type, private :: rec1
		sequence
		character(len=4)	::	trj_type
		integer(4)				::	n_frames, n_steps_before, interval, n_steps
		integer(4)				::	interval_velcheck, unused6, unused7, n_degf
		integer(4)				::	n_fixed, unused10, const_p, unused12
		integer(4)				::	unused13, unused14, unused15, unused16
		integer(4)				::	unused17, unused18, unused19, charmm_version
	end type rec1

	type(rec1), private					::	r1

	character(len=80), private	::	topology
	integer, private						::	mask_rows
	integer, parameter					::	max_mask_rows = 20
	character(len=80)						::	mask_row(max_mask_rows)

	integer, private						::	lun = 0
	integer, private						::	current_frame = 0
		
	!mask 	
	integer, private							::	ncoords
	real(4), private, allocatable	::	xmasked(:)
	type(mask_type), private			::	mask
contains


subroutine trj_startup
	mask_rows = 0
	call mask_startup
end subroutine trj_startup


subroutine trj_shutdown
	call trj_close
	call mask_shutdown
end subroutine trj_shutdown


subroutine trj_initialize(frames, steps_before, interval, steps,&
	degf, topfile)
!arguments
	integer, intent(in)			::	frames, steps_before, interval, steps
	integer, intent(in)			::	degf
	character(*), intent(in)	::	topfile

	!clear file number
	lun = 0
	
	!set defaults
	r1%trj_type = 'cord' !coordinate trajectory (not velocity)
	r1%interval_velcheck = 0
	r1%const_p = 0
	r1%charmm_version = 0
	r1%n_fixed = 0 !required by vmd program

	!set values
	r1%n_frames = frames
	r1%n_steps_before = steps_before
	r1%interval = interval
	r1%n_steps = steps
	r1%n_degf = degf
	topology = topfile

	call mask_initialize(mask)

end subroutine trj_initialize


integer function trj_add(line)
	!arguments
	character(*)				::	line
	
	if(lun /= 0) then
		!file open - can't add now
		write(*,910) 
910		format('>>>>> error: trajectory file open, cannot add atoms to mask')
		trj_add = 0
	elseif(mask_rows < max_mask_rows) then
		mask_rows = mask_rows + 1
		mask_row(mask_rows) = line
		trj_add = mask_add(mask, line)
	else
		write(*,900) max_mask_rows
900		format('>>>>> error: too many trajectory atom mask rows (max',i3,').')
		trj_add = 0
	end if
end function trj_add

logical function trj_store_mask(line)
	!arguments
	character(*)				::	line
	
	if(mask_rows < max_mask_rows) then
		mask_rows = mask_rows + 1
		mask_row(mask_rows) = line
		trj_store_mask = .true.
	else
		write(*,900) max_mask_rows
900		format('>>>>> error: too many trajectory atom mask rows (max',i3,').')
		trj_store_mask = .false.
	end if
end function trj_store_mask


integer function trj_commit_mask()
	!locals
	integer						::	row

	if(lun /= 0) then
		!file open - can't add now
		write(*,910) 
910		format('>>>>> error: trajectory file open, cannot add atoms to mask')
	else
		do row = 1, mask_rows
			trj_commit_mask = mask_add(mask, mask_row(row))
		end do
	end if
	trj_commit_mask = mask%included
end function trj_commit_mask


integer function trj_count()
	trj_count = mask%included
end function trj_count

logical function trj_create(filename, append)
!arguments
	character(*)				::	filename
	logical, optional			::	append
!locals
	integer						::	filestat
	character(len=10)			::	writemode
	character(len=80)			::	titlerow

	!allocated masked coordinate array
	ncoords = 3*mask%included
	allocate(xmasked(ncoords))

	titlerow = 'q dcd trajectory version ' // module_version

	writemode = 'rewind'
	if(present(append)) then
		if(append) writemode='append'
	end if

	lun = freefile()
	open(unit=lun, file=filename, status='unknown', form='unformatted', action='write',&
		position=writemode, err=900)
	if(writemode /= 'append') then
		!write the first 3 dcd format records
		write(lun, err=920) r1
		write(lun, err=920) 2+mask_rows, titlerow, topology, mask_row(1:mask_rows)
		write(lun, err=920) int(mask%included, 4)
	end if
	trj_create = .true.
	return
	
	!error handling
	!open error
900	write(*,910) trim(filename)
910	format('>>>>> error: failed to open trajectory ',a)
	trj_create = .false.
	return
	
	!write header error
920	write(*,930) 
930	format('>>>>> error: failed to write trajectory header.')
	trj_create = .false.
	return

end function trj_create


!******************************************************
!write real(4) coords to trajectory. writes only 
! atoms in current mask.
!******************************************************
logical function trj_write(x)
!arguments
	real(8) 					::	x(:)
	
	!extract coordinates
	call mask_get(mask, x, xmasked)
	!write x record
	write(lun, err=900) xmasked(1:ncoords:3)
	!write y record
	write(lun, err=900) xmasked(2:ncoords:3)
	!write z record
	write(lun, err=900) xmasked(3:ncoords:3)

	trj_write = .true.
	return

	!error handling
900	trj_write = .false.
	
end function trj_write

!******************************************************
!read a frame from tjacectory file and returns coordinates
! associated with atomnumbers in topology.
!******************************************************
logical function trj_read(x)
!arguments
	real(8) 					::	x(:)
	!read x record to temp variable xmasked
	read(lun, err=900, end=900) xmasked(1:ncoords:3)
	!read y record
	read(lun, err=900, end=900) xmasked(2:ncoords:3)
	!read z record
	read(lun, err=900, end=900) xmasked(3:ncoords:3)


	!assign masked coordinates to right atom in topology
	call mask_put(mask, x, xmasked)

	trj_read = .true.

	current_frame = current_frame + 1
	return

	!error handling
900	trj_read = .false.
	
end function trj_read

!******************************************************
!read a frame from tjacectory file containing masked atoms,
! returns only masked coordinates. !real(4)!
!******************************************************
logical function trj_read_masked(x)
!arguments
	real(4) 					::	x(:)
	!read x record to temp variable xmasked
	read(lun, err=900, end=900) x(1:ncoords:3)
	!read y record
	read(lun, err=900, end=900) x(2:ncoords:3)
	!read z record
	read(lun, err=900, end=900) x(3:ncoords:3)

	trj_read_masked = .true.

	current_frame = current_frame + 1
	return

	!error handling
900	trj_read_masked = .false.
	
end function trj_read_masked

!******************************************************

logical function trj_seek(frame)
!arguments
	integer						::	frame
!locals
	integer						::	i

	!optionally fast-forward to selected frame
	if(frame > r1%n_frames .or. frame < 1) then
		write(*,962) frame
962		format('>>>>> error: frame',i6,' does not exist.')
		close(lun)
		trj_seek = .false.
	elseif(frame > current_frame) then
		do i=current_frame+1, frame
			!read x y and z records
			read(lun, err=980)
			read(lun, err=980)
			read(lun, err=980)
		end do
		current_frame = frame
		trj_seek = .true.
	elseif(frame < current_frame) then
		do i=current_frame, frame+1, - 1
			!read x y and z records
			backspace(lun, err=980)
			backspace(lun, err=980)
			backspace(lun, err=980)
		end do
		current_frame = frame
		trj_seek = .true.
	else
		!already correctly positiones
		trj_seek = .true.
	end if
	return

980	write(*,990) i, frame
990	format('>>>>> error: read error at frame',i6,' while looking for frame',i6)
	trj_seek = .false.
	return

end function trj_seek

subroutine trj_close
!locals
	logical						::	is_open
	inquire(lun, opened=is_open)
	if(is_open) close(lun)
	lun = 0
	if(allocated(xmasked)) deallocate(xmasked)
	call mask_finalize(mask)
	ncoords = 0
	current_frame = 0
end subroutine trj_close

!***********************************************

logical function trj_open(filename)
!arguments
	character(*)				::	filename

!locals
	character(len=80)			::	titlerow
	integer						::	atoms
	integer						::	i

	!open file
	lun=freefile()
	open(unit=lun, file=filename, status='old', form='unformatted', action='read',&
		err=900)
	!read header
	read(lun, err=920) r1
	if(r1%trj_type /= 'cord') then
		write(*,940)
940		format('>>>>> error: not a dcd coordinate trajectory.')
		close(lun)
		trj_open = .false.
		return
	end if
!	write(*,456) r1%n_frames
456 format(i10)
	read(lun, err=920) mask_rows, titlerow, topology, mask_row(1:mask_rows-2)
	if(titlerow(1:16) /= 'q dcd trajectory') then
		write(*,941)
941		format('>>>>> error: not a q dcd trajectory.')
		close(lun)
		trj_open = .false.
		return
	end if
		
	mask_rows = mask_rows - 2
	read(lun, err=920) atoms
	
	!re-create mask
	call mask_initialize(mask)
	if(allocated(xmasked)) deallocate(xmasked)

	trj_open = .true.
	do i = 1, mask_rows
		if(mask_add(mask, mask_row(i)) == 0) then
			write(*,950) trim(mask_row(i))
		end if
	end do
950	format('>>> warning: no atoms in mask ',a)
	if(mask%included /= atoms) then
		write(*,960) atoms, mask%included
		write(*,*) 'ignored'
!		close(lun)
!		trj_open = .false.
!		return
	end if
960	format('>>>>> error: different number of atoms in mask and trajectory',/,&
			i6, ' atoms in trajectory,',i6,' atoms in mask.')

	!allocated masked coordinate array
	ncoords = 3*atoms
	allocate(xmasked(ncoords))

	current_frame = 1

	return
	
	!error handling
	!open error
900	write(*,910) trim(filename)
910	format('>>>>> error: failed to open trajectory ',a)
	trj_open = .false.
	return
	
	!read header error
920	write(*,930) 
930	format('>>>>> error: failed to read trajectory header.')
	trj_open = .false.
	close(lun)
	return

end function trj_open


integer function trj_intersection(m)
!arguments
	type(mask_type)				::	m
	trj_intersection = count(m%mask .and. mask%mask)
end function trj_intersection


subroutine trj_clone_mask(m)
!arguments
	type(mask_type)				::	m
	if(size(m%mask) == size(mask%mask)) then
		m = mask
	end if
end subroutine trj_clone_mask	

integer function trj_get_ncoords()
 trj_get_ncoords = ncoords
end function trj_get_ncoords	

end module trj
