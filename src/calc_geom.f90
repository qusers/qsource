! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! calc_geom.f90
! by John Marelius
! geometry analysis functions
!TODO: precision not fixed

module CALC_GEOM
	use CALC_BASE
	implicit none

!constants
	integer, parameter			::	MAX_MEASUREMENTS = 99

!module variables
	real(8)						:: pi,deg2rad,rad2deg
	integer, private			::	Nmeas = 0
	type GEOM_TYPE
		integer					::	i, j, k, l, cod, qcod(max_states)
		logical					::	is_fep
		integer					::	kind
	end type GEOM_TYPE
	type(GEOM_TYPE), private	::	geom(MAX_MEASUREMENTS)
contains

subroutine geom_initialize
	pi = 4 * atan(1.)
    deg2rad = pi/180.
    rad2deg = 1./deg2rad
end subroutine geom_initialize

subroutine geom_finalize
end subroutine geom_finalize

integer function dist_add(desc)
	!arguments
	character(*)				::	desc

	!locals
	character(len=80)			::	line
	integer						::	readstat
	integer						::	b

	if(Nmeas == MAX_MEASUREMENTS) then
		write(*,10) MAX_MEASUREMENTS
		return
	end if
10	format('Sorry, the maximum number of geometry calculations is ',i2)
	!add a new dist measurement
	Nmeas = Nmeas + 1
	write(*,'(a)', advance='no') 'dist: Enter the two atom numbers: '
	read(*,'(a)', iostat=readstat) line
	read(line, *, iostat=readstat) geom(Nmeas)%i, geom(Nmeas)%j
	if(readstat /= 0 .or. geom(Nmeas)%i < 1 .or. geom(Nmeas)%i > nat_pro .or. &
		geom(Nmeas)%j < 1 .or. geom(Nmeas)%j > nat_pro) then
		write(*, 900) 
900		format('>>>>> ERROR: Invalid atom number(s).')
		Nmeas = Nmeas - 1
		dist_add = 0
		return
	end if

	geom(Nmeas)%cod = 0 !clear force field code
	geom(Nmeas)%qcod(:) = 0
	geom(Nmeas)%is_fep = .false.
	!search topology for a bond between these atoms
	do b = 1, nbonds
		if((bnd(b)%i == geom(Nmeas)%i .and. bnd(b)%j == geom(Nmeas)%j) .or. &
			(bnd(b)%j == geom(Nmeas)%i .and. bnd(b)%i == geom(Nmeas)%j)) then
			geom(Nmeas)%cod = bnd(b)%cod
			exit
		end if
	end do
	!search FEP strategy for bond between the atoms
	if (use_fep) then
	!only if topology, else boom
	if ((geom(Nmeas)%i.le.nat_solute).and.(geom(Nmeas)%j.le.nat_solute)) then
		do b = 1, nqbond
			if((qbnd(b)%i==geom(Nmeas)%i .and. qbnd(b)%j==geom(Nmeas)%j) .or. &
				(qbnd(b)%j==geom(Nmeas)%i .and. qbnd(b)%i==geom(Nmeas)%j)) then
			do states = 1, nstates
			geom(Nmeas)%qcod(states) = qbnd(b)%cod(states)
			end do
			geom(Nmeas)%cod = 0
			geom(Nmeas)%is_fep = .true.
			exit
			end if 
		end do
	end if
	end if
	dist_add = Nmeas
	if(geom(Nmeas)%cod > 0) then
		write(desc, 19) geom(Nmeas)%i, geom(Nmeas)%j
19		format('distance, bond energy between atoms',i5,' and',i5)
	elseif(geom(Nmeas)%is_fep) then
		write(desc, 219) geom(Nmeas)%i, geom(Nmeas)%j
219             format('distance, qbond energy between atoms',i5,' and',i5)
	else
		write(desc, 20) geom(Nmeas)%i, geom(Nmeas)%j
20		format('distance between atoms',i5,' and',i5)
	end if

end function dist_add


integer function angle_add(desc)
	!arguments
	character(*)				::	desc

	!locals
	character(len=80)			::	line
	integer						::	readstat
	integer						::	i, j, k
	integer						::	b

	if(Nmeas == MAX_MEASUREMENTS) then
		write(*,10) MAX_MEASUREMENTS
		return
	end if
10	format('Sorry, the maximum number of geometry calculations is ',i2)
	!add a new dist measurement
	Nmeas = Nmeas + 1
	write(*,'(a)', advance='no') 'angle: Enter the three atom numbers: '
	read(*,'(a)', iostat=readstat) line
	read(line, *, iostat=readstat) i, j, k
	if(readstat /= 0 .or. i < 1 .or. i > nat_pro .or. &
		j < 1 .or. j > nat_pro .or. &
		k < 1 .or. k > nat_pro .or. &
		i == j .or. i == k .or. j == k) then
		write(*, 900) 
900		format('>>>>> ERROR: Invalid atom number(s).')
		Nmeas = Nmeas - 1
		angle_add = 0
		return
	end if

	geom(Nmeas)%i = i
	geom(Nmeas)%j = j
	geom(Nmeas)%k = k
	geom(Nmeas)%cod = 0 !clear force field code
	geom(Nmeas)%qcod(:) = 0
	geom(Nmeas)%is_fep = .false.
	!search topology for this angle
	do b = 1, nangles
		if(ang(b)%j == j .and. &
			((ang(b)%i == i .and. ang(b)%k == k) .or. &
			 (ang(b)%i == k .and. ang(b)%k == i))) then
			geom(Nmeas)%cod = ang(b)%cod
			exit
		end if
	end do
        !search FEP strategy for angle between the atoms
        if (use_fep) then
        !only if topology, else boom
        if ((geom(Nmeas)%i.le.nat_solute).and.(geom(Nmeas)%k.le.nat_solute)) then
                do b = 1, nqangle
                        if((qang(b)%i==geom(Nmeas)%i .and. qang(b)%k==geom(Nmeas)%k) .or. &
                                (qang(b)%i==geom(Nmeas)%k .and. qang(b)%k==geom(Nmeas)%i)) then
                        do states = 1, nstates
                        geom(Nmeas)%qcod(states) = qang(b)%cod(states)
                        end do
                        geom(Nmeas)%cod = 0
                        geom(Nmeas)%is_fep = .true.
                        exit
                        end if
                end do
        end if
        end if
	angle_add = Nmeas
	if(geom(Nmeas)%cod > 0) then
		write(desc, 19) i, j, k
19		format('angle, angle energy between atoms',i5,',',i5,',',i5)
        elseif(geom(Nmeas)%is_fep) then
                write(desc, 319) i, j, k
319             format('angle, qangle energy between atoms',i5,',',i5,',',i5)
	else
		write(desc, 20) i, j, k
20		format('angle between atoms',i5,',',i5,',',i5)
	end if

end function angle_add


integer function torsion_add(desc)
	!arguments
	character(*)				::	desc

	!locals
	character(len=80)			::	line
	integer						::	readstat
	integer						::	i, j, k, l
	integer						::	b

	if(Nmeas == MAX_MEASUREMENTS) then
		write(*,10) MAX_MEASUREMENTS
		return
	end if
10	format('Sorry, the maximum number of geometry calculations is ',i2)
	!add a new dist measurement
	Nmeas = Nmeas + 1
	write(*,'(a)', advance='no') 'torsion: Enter the four atom numbers: '
	read(*,'(a)', iostat=readstat) line
	read(line, *, iostat=readstat) i, j, k, l
	if(readstat /= 0 .or. i < 1 .or. i > nat_pro .or. &
		j < 1 .or. j > nat_pro .or. &
		k < 1 .or. k > nat_pro .or. &
		l < 1 .or. l > nat_pro .or. &
		i == j .or. i == k .or. i == l .or. j == k .or. &
		j == l .or. k == l) then
		write(*, 900) 
900		format('>>>>> ERROR: Invalid atom number(s).')
		Nmeas = Nmeas - 1
		torsion_add = 0
		return
	end if

	geom(Nmeas)%i = i
	geom(Nmeas)%j = j
	geom(Nmeas)%k = k
	geom(Nmeas)%l = l
	geom(Nmeas)%cod = 0 !clear force field code
	geom(Nmeas)%qcod(:)=0
	geom(Nmeas)%is_fep = .false.
	!search topology for this angle
	do b = 1, ntors
		if((tor(b)%i == i .and. tor(b)%j == j .and. &
			tor(b)%k == k .and. tor(b)%l == l) .or. &
		   (tor(b)%i == l .and. tor(b)%j == k .and. &
			tor(b)%k == j .and. tor(b)%l == i)) then
			geom(Nmeas)%cod = tor(b)%cod
			exit
		end if
	end do
        !search FEP strategy for torsion between the atoms
        if (use_fep) then
        !only if topology, else boom
        if ((geom(Nmeas)%i.le.nat_solute).and.(geom(Nmeas)%j.le.nat_solute) .and. &
		(geom(Nmeas)%k.le.nat_solute).and.(geom(Nmeas)%l.le.nat_solute)) then
                do b = 1, nqtor
                if((qtor(b)%i == i .and. qtor(b)%j == j .and. &
                        qtor(b)%k == k .and. qtor(b)%l == l) .or. &
                   (qtor(b)%i == l .and. qtor(b)%j == k .and. &
                        qtor(b)%k == j .and. qtor(b)%l == i)) then
                        do states = 1, nstates
                        geom(Nmeas)%qcod(states) = qtor(b)%cod(states)
                        end do
                        geom(Nmeas)%cod = 0
                        geom(Nmeas)%is_fep = .true.
                        exit
                        end if
                end do
        end if
        end if
	write(*,*)i, j, k, l
	torsion_add = Nmeas
	if(geom(Nmeas)%cod > 0) then
		write(desc, 19) i, j, k, l
19		format('torsion, torsion energy between atoms',i5,',',i5,',',i5,',',i5)
        elseif(geom(Nmeas)%is_fep) then
                write(*, 419) i, j, k, l
419             format('torsion, qtorsion energy between atoms',i5,',',i5,',',i5,',',i5)
	else
		write(desc, 20) i, j, k, l
20		format('torsion between atoms',i5,',',i5,',',i5,',',i5)
	end if

end function torsion_add


subroutine dist_calc(i)
	!arguments
	integer, intent(in)			::	i

	!locals
	real(8)						::	rji(3)
	real(8)						::	r2, r, V , fexp
	V = 0
	if(i < 1 .or. i > Nmeas) return
	
	rji(:) = xin(3*geom(i)%i-2:3*geom(i)%i) - xin(3*geom(i)%j-2:3*geom(i)%j)
	r2 = dot_product(rji, rji)
	r = sqrt(r2)
	write(*,100, advance='no') r
100	format(f10.2)
	if(geom(i)%cod > 0) then !calc energy
		V = 0.5 * bondlib(geom(i)%cod)%fk * (r - bondlib(geom(i)%cod)%bnd0)**2
		write(*,110, advance='no') V
        end if
	if (geom(i)%is_fep) then !calc q energy
		do states = 1 , nstates
			if(geom(i)%qcod(states)>0) then
			        fexp = exp(-qbondlib(geom(i)%qcod(states))%amz*(r-qbondlib(geom(i)%qcod(states))%r0))
			        V = V + (lamda(states) &
				      * qbondlib(geom(i)%qcod(states))%Dmz*(fexp*fexp-2.*fexp) &
				      + 0.5*qbondlib(geom(i)%qcod(states))%fk*(r-qbondlib(geom(i)%qcod(states))%r0)**2)
			end if
		end do
		write(*,110, advance='no') V
	end if
110	format(f8.2)
	
end subroutine dist_calc


subroutine angle_calc(i)
	!arguments
	integer, intent(in)			::	i

	!locals
	real(8)						::	rji(3), rjk(3), scp
	real(8)						::	a, V
	V = 0
	if(i < 1 .or. i > Nmeas) return
	
	rji(:) = xin(3*geom(i)%i-2:3*geom(i)%i) - xin(3*geom(i)%j-2:3*geom(i)%j)
	rjk(:) = xin(3*geom(i)%k-2:3*geom(i)%k) - xin(3*geom(i)%j-2:3*geom(i)%j)
	scp = dot_product(rji, rjk) &
		/ sqrt(dot_product(rji, rji)*dot_product(rjk, rjk))
	IF(scp>1.) scp = 1.
	IF(scp< -1.) scp = -1.
	a = acos(scp)
	write(*,100, advance='no') a*rad2deg
100	format(f11.2)

	if(geom(i)%cod > 0) then !calc energy
		V = 0.5 * anglib(geom(i)%cod)%fk &
			*(a-anglib(geom(i)%cod)%ang0*deg2rad)**2
		write(*,110, advance='no') V
        end if
	if (geom(i)%is_fep) then !calc q energy
		do states = 1 , nstates
                        if(geom(i)%qcod(states)>0) then
                                V = V + ( lamda(states) &
                                      *   (0.5*qanglib(geom(i)%qcod(states))%fk &
                                      *   (a-qanglib(geom(i)%qcod(states))%ang0)**2) )
! no conversion factor here because qangles are already converted in
! qatom_load_fep
                        end if
                end do
		write(*,100, advance='no') V
	end if
110	format(f8.2)
	
end subroutine angle_calc


subroutine torsion_calc(i)
	!arguments
	integer, intent(in)			::	i

	!locals
	real(8)						::	rji(3), rjk(3), rkl(3), rnj(3), rnk(3), scp
	real(8)						::	phi, sgn, arg, V
	integer						::	ic
	V = 0
	if(i < 1 .or. i > Nmeas) return
	
	rji(:) = xin(3*geom(i)%i-2:3*geom(i)%i) - xin(3*geom(i)%j-2:3*geom(i)%j)
	rjk(:) = xin(3*geom(i)%k-2:3*geom(i)%k) - xin(3*geom(i)%j-2:3*geom(i)%j)
	rkl(:) = xin(3*geom(i)%l-2:3*geom(i)%l) - xin(3*geom(i)%k-2:3*geom(i)%k)
	!cross products
	rnj(1) = rji(2) * rjk(3) - rji(3) * rjk(2)
	rnj(2) = rji(3) * rjk(1) - rji(1) * rjk(3)
	rnj(3) = rji(1) * rjk(2) - rji(2) * rjk(1)
	rnk(1) = - rjk(2) * rkl(3) + rjk(3) * rkl(2)
	rnk(2) = - rjk(3) * rkl(1) + rjk(1) * rkl(3)
	rnk(3) = - rjk(1) * rkl(2) + rjk(2) * rkl(1)
	scp = dot_product(rnj, rnk) &
		/ sqrt(dot_product(rnj, rnj)*dot_product(rnk, rnk))
	IF(scp>1.0) scp = 1.0
	IF(scp< - 1.0) scp = - 1.0
	phi = acos(scp)
	sgn = rjk(1) *(rnj(2) * rnk(3) - rnj(3) * rnk(2) ) + rjk(2)&
	*(rnj(3) * rnk(1) - rnj(1) * rnk(3) ) + rjk(3) *(rnj(1)   &
	* rnk(2) - rnj(2) * rnk(1) )
	IF(sgn<0) phi = - phi
	write(*,100, advance='no') phi*180/pi
100	format(f10.2)

	ic = geom(i)%cod
	if(ic > 0) then !calc energy
		arg = torlib(ic)%rmult * phi - torlib(ic)%deltor * pi / 180.
		V = torlib(ic)%fk *(1.0 + cos(arg) ) / real(torlib(ic)%paths )
		write(*,110, advance='no') V
        end if
        if (geom(i)%is_fep) then !calc q energy
                do states = 1 , nstates
                        if(geom(i)%qcod(states)>0) then
			        arg = qtorlib(ic)%rmult*phi-qtorlib(ic)%deltor
                                V = V + ( lamda(states) &
                                      * qtorlib(ic)%fk*(1.0+cos(arg)) )
                        end if
                end do
                write(*,100, advance='no') V
	end if
110	format(f8.2)
	
end subroutine torsion_calc


subroutine dist_heading(i)
	integer						::	i
	
	write(*,'(a)', advance='no') 'dist(Ang)'

	if((geom(i)%cod > 0).or.(any(geom(i)%qcod(:).gt.0))) then
		write(*,'(a8)', advance='no') 'V_bond'
	end if
end subroutine dist_heading


subroutine angle_heading(i)
	integer						::	i
	
	write(*,'(a)', advance='no') 'angle(Deg)'

	if((geom(i)%cod > 0).or.(any(geom(i)%qcod(:).gt.0))) then
		write(*,'(a8)', advance='no') 'V_angle'
	end if
end subroutine angle_heading


subroutine torsion_heading(i)
	integer						::	i
	
	write(*,'(a)', advance='no') 'tors(Deg)'

	if((geom(i)%cod > 0).or.(any(geom(i)%qcod(:).gt.0))) then
		write(*,'(a8)', advance='no') 'V_tors'
	end if
end subroutine torsion_heading


end module CALC_GEOM
