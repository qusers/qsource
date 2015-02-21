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
!  (c) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden                       !
!  qfep.f90                                                                    !
!  by Johan Aqvist, Karin Kolmodin & John Marelius                             !
!  qfep free energy analysis program for fep, evb & umbrella sampling          !
!------------------------------------------------------------------------------!
program qfep
use nrgy        
use parse

implicit none
        character(*), parameter                 ::      module_version = '5.01'
        character(*), parameter                 ::      module_date = '2003-06-03'

        integer,parameter ::mxpts=20000,mxbin=100,mxstates=4
        character(80)      ::filnam, line
        integer           ::i,j,ifile,ipt,istate,ibin,nfiles,nstates,err, &
                              nskip,nbins,nmin,idum,noffd,nnoffd,offel

        type(offdiag_save), dimension(mxstates) :: offd

        real(8) ::rt,gapmin,gapmax,sum,dv,gaprange, &
                        xint,dvg,veff1,veff2,dga,dgb,dgg,alpha_b,scale_hij, &
                        veff,min                        
        real(8),dimension(mxbin)              ::sumg,sumg2,avdvg,avc1,avc2,avr
                                                                                  
        real(8),dimension(mxbin,4)            ::binsum
        integer,dimension(mxbin)              ::nbinpts,ptsum

        type(q_energies), dimension(mxstates)   :: eq   
        type(q_energies), dimension(mxstates)   :: aveq
        real(8),dimension(mxstates)                             :: dvv,dgv,alfa,coeff

        real(8),dimension(3)                  ::u,y
        real(8),allocatable              ::dgf(:),dgr(:),dgfsum(:),dgrsum(:),dg(:)
        real(8),dimension(mxstates,mxstates)  ::a,mu,eta,rxy0
        real(8),allocatable                   ::hij(:,:),d(:),e(:)
        type fep_data_type
                integer                                 ::      npts
                real(8)                                 ::      lambda(mxstates)
                real(8), pointer                ::      v(:,:), r(:,:) !indices are state, point
                real(8), pointer                ::      vg(:), gap(:), c1(:), c2(:) !index is point
        end type fep_data_type
        type(fep_data_type), allocatable        ::      fep(:) !index is file
        type(fep_data_type)                             ::      feptmp !temporary storage for one file

        real(8),dimension(mxstates,mxbin)     ::avdvv,sumv,sumv2 

        integer                                                         ::      f

        !header
        write(*,100) module_version,  module_date
        write(*,*)
100     format('# qfep',t30,'version ',a,t50,'(modified on ',a,')')
        
        !------------------------------------------ 
        ! input of parameters
        !------------------------------------------
        call prompt ('--> number of energy files: ')
        read (*,*) nfiles
        write (*,1) nfiles
1       format('# number of files                 =',i6)
        call prompt ('--> no. of states, no. of predefined off-diag elements: ')
        read (*,*) nstates, noffd
        write (*,2) nstates, noffd
2       format('# number of states                 =',i6,/, &
                   '# number of off-diagonal elements =',i6)


        !size of secular determinant is allocated

        allocate(hij(nstates,nstates),d(nstates),e(nstates),stat=err)
        if(err /= 0) then
                write(*,*) 'error: out of memory when allocation hij array.'
                stop 'qfep5 terminated abnormally: out of memory.'
        end if

        ! continue to read input
        call prompt ('--> give kt & no, of pts to skip: ')
        read (*,*) rt,nskip
        write (*,3) rt,nskip
3       format('# kt                              =',f6.3,/, &
                   '# number of data points to skip   =',i6)

        call prompt ('--> give number of gap-bins: ')
        read (*,*) nbins
        write (*,4) nbins
4       format('# number of gap-bins              =',i6)

        call prompt ('--> give minimum # pts/bin: ')
        read (*,*) nmin
        write (*,5) nmin
5       format('# minimum number of points per bin=',i6)

        do istate=2,nstates
                write(line,7) istate
                call prompt(line)
                read (*,*) alfa(istate)
                write (*,6) istate,alfa(istate)
        end do
6       format('# alpha for state ',i2,'              =',f6.2)
7       format('--> give alpha for state ',i2,':')

        scale_hij=0.0
        if (noffd /=0) then
                call prompt ('--> hij scaling:')
                read (*,*) scale_hij
                write (*,8) scale_hij
8       format('# scale factor for hij            =',f6.2)
        else
                call prompt ('--> number of off-diagonal elements:')
                read (*,*) nnoffd
                write (*,9) nnoffd
                if(nnoffd >0) write(*,11) 
                do offel=1,nnoffd
                        call prompt ('--> i, j, a_ij, mu_ij, eta_ij, r_xy0: ')
                        read (*,*) i, j, a(i,j), mu(i,j), eta(i,j), rxy0(i,j)
                        write(*,12) i, j, a(i,j), mu(i,j), eta(i,j), rxy0(i,j)
                end do
9       format('# number of off-diagonal elements =',i6)
11      format('#   i   j   a(i,j)  mu(i,j) eta(i,j) rxy0(i,j)')
12      format('#',2i4,4f9.2) 
        end if

        call prompt ('--> linear combination of states defining reaction coord: ')
        read (*,*) (coeff(i),i=1,nstates)
        write(*,13) coeff(1:nstates)            
13      format('# linear combination co-efficients=',8f6.2)

        !allocate large arrays
        allocate(fep(nfiles), feptmp%v(nstates,mxpts), feptmp%r(nstates,mxpts), &
                feptmp%vg(mxpts), feptmp%gap(mxpts), &
                feptmp%c1(mxpts), feptmp%c2(mxpts), &
                dgf(0:nfiles+1),dgr(0:nfiles+1), &
                dgfsum(0:nfiles+1),dgrsum(0:nfiles+1),dg(0:nfiles+1), &
                stat=err)
        if(err /= 0) then
                stop 'qfep5 terminated abnormally: out of memory when allocating arrays.'
        end if

!---------------------------------
! energy files are opened and read
!---------------------------------
        binsum(:,:)=0.
        eq(:)%total=0.

        gapmin=999.
        gapmax=-999.
        f = freefile()
        feptmp%c1(:) = 0.
        feptmp%c2(:) = 0.

        write(*,*)
        write(*,*)
        write(*,15)
        write(*,16)
15      format('# part 0: average energies for all states in all files')
16      format('# file             state   pts   lambda    eqtot   eqbond',&
        '  eqang   eqtor   eqimp    eqel   eqvdw  eel_qq  evdw_qq eel_qp  evdw_qp eel_qw evdw_qw eqrstr')

        do ifile=1,nfiles
                write(line,14) ifile
14              format('--> name of file number',i4,':')
                call prompt(line)
                read(*,*) filnam
                write (*,*) ''
                if(openit(f,filnam,'old','unformatted','read') /= 0) then
                        stop 'qfep5 terminated abnormally: failed to open energy file.'
                end if

                !read 1st record to get lambdas
                idum = get_ene(f, eq(:), offd, nstates,nnoffd)
                if(idum /= 0) then 
                        !we have a problem
                        write(*,'(a,a)') 'error: unexpected end-of-file in first record of ', &
                                trim(filnam)
                        if(idum > 0) then
                        !the problem is in energies
                                write(*,'(a,i1)') 'while reading energies of state ',idum
                                write(*,'(a,i1,a)') 'maybe there are only ',idum-1, ' states?'
                                stop 'qfep5 terminated abnormally: failed to read energies.'
                        else
                        !idum < 0 indicates problems with offdiags
                                write(*,'(a)') 'while reading off-diagonal elements.'
                                write(*,'(a,i1,a)') 'maybe there are less than ',nnoffd, &
                                        ' off-diagonal elements?'
                                stop 'qfep5 terminated abnormally: failed to read off-diagonals.'
                        end if
                end if
                                                                
                fep(ifile)%lambda(:) = eq(:)%lambda

        rewind (f)

                ipt = 0
                aveq(:) = aveq(:) * 0. !set all fields to zero using multiplication operator
                feptmp%gap(:) = 0.
                do while(get_ene(f, eq(:), offd, nstates, nnoffd) == 0) !keep reading till eof
                        ipt = ipt + 1
                        if(ipt > nskip) then
                                aveq(:) = aveq(:) + eq(:) !use qenergies + operator 
                        end if

!-------------------------------------------
! correct h_ii with alfa, and modify h_ij...
!-------------------------------------------

            alfa(1)=0.
                        eq(:)%total=eq(:)%total+alfa(:)
                        feptmp%v(1:nstates, ipt) = eq(1:nstates)%total
                        if (nnoffd .ne. 0) then
                                feptmp%r(:,ipt) = offd(:)%rkl
            end if


                        do i=1, noffd
                                hij(offd(i)%i, offd(i)%j) = offd(i)%hij
                        end do

                        if ( scale_hij .gt. 0.0 ) then
                                hij(1,2) = scale_hij*hij(1,2)
                        else
                                do i=1,nstates
                                        do j=1,nstates
                                                if (i==j) then
                                                        hij(i,j)=eq(i)%total
                                                else
                                                        if (a(i,j)==0.0) then
                                                                hij(i,j) = a(j,i)*exp(-mu(j,i)*(offd(1)%rkl-rxy0(j,i)))* &
                                                                exp(-eta(j,i)*(offd(1)%rkl-rxy0(j,i))**2)                                                                               
                                                        else 
                                                                hij(i,j) = a(i,j)*exp(-mu(i,j)*(offd(1)%rkl-rxy0(i,j)))* &
                                                                exp(-eta(i,j)*(offd(1)%rkl-rxy0(i,j))**2)                                                                        
                                                        end if 
                                                end if
                                        end do
                                end do
            end if

!-----------------------------------------------------------
! ground state energy is calculated from secular determinant
!-----------------------------------------------------------

                        if (nstates==2) then
                                feptmp%vg(ipt)=0.5*(eq(1)%total+eq(2)%total)-  &
                                        0.5*sqrt( (eq(1)%total-eq(2)%total)**2 + 4.*hij(1,2)**2 )
                                if(nnoffd > 0) then
                                        feptmp%c1(ipt)=1./(1.+((feptmp%vg(ipt)-eq(1)%total)/hij(1,2))**2)
                                        feptmp%c2(ipt)=1-feptmp%c1(ipt)
                                end if
                        else 
                                call tred2(hij,nstates,nstates,d,e)
                                call tqli(d,e,nstates,nstates,hij)
                                feptmp%vg(ipt)=minval(d)
                        end if 

                        do istate=1,nstates
                                feptmp%gap(ipt)=feptmp%gap(ipt)+feptmp%v(istate,ipt)*coeff(istate)
                        end do


              if(ipt .gt. nskip) then
                 if(feptmp%gap(ipt) .lt. gapmin) gapmin=feptmp%gap(ipt)
                 if(feptmp%gap(ipt) .gt. gapmax) gapmax=feptmp%gap(ipt)
              end if
         end do  !(ipt)
                 close(f)
                fep(ifile)%npts = ipt

                !copy feptmp to d
                allocate(fep(ifile)%v(nstates, fep(ifile)%npts), &
                        fep(ifile)%vg(fep(ifile)%npts), &
                        fep(ifile)%gap(fep(ifile)%npts), &
                        fep(ifile)%c1(fep(ifile)%npts), &
                        fep(ifile)%c2(fep(ifile)%npts), &
                        stat=err)
                if(err /= 0) then
                        stop 'qfep5 terminated abnormally: out of memory when allocating arrays.'
                end if
                fep(ifile)%v(:,:) = feptmp%v(1:nstates, 1:fep(ifile)%npts)
                fep(ifile)%vg(:) = feptmp%vg(1:fep(ifile)%npts)
                fep(ifile)%gap(:) = feptmp%gap(1:fep(ifile)%npts)
                fep(ifile)%c1(:) = feptmp%c1(1:fep(ifile)%npts)
                fep(ifile)%c2(:) = feptmp%c2(1:fep(ifile)%npts)
                if(nnoffd > 0) then 
                        allocate(fep(ifile)%r(nstates, fep(ifile)%npts), stat=err)
                        if(err /= 0) then
                                stop 'qfep5 terminated abnormally: out of memory when allocating arrays.'
                        end if
                        fep(ifile)%r(:,:) = feptmp%r(1:nstates, 1:fep(ifile)%npts)
                end if

        
                ! average energies for each file is calculated

                if(ipt <= nskip) then !skipped all points in the file
                        write(*,900) trim(filnam)
900                     format('>>>>> error: number of data sets in ',a,&
                                        ' is less than number of points to skip!')
                        stop 'qfep5 terminated abnormally: too few data points in file.'
                else
                        aveq(:) = aveq(:) * (1./(fep(ifile)%npts-nskip)) !use new * operator
                end if
                do istate=1,nstates
                write(*,17) trim(filnam), istate, fep(ifile)%npts-nskip, fep(ifile)%lambda(istate), &
                        aveq(istate)%total,aveq(istate)%q%bond,aveq(istate)%q%angle,&
                        aveq(istate)%q%torsion, &
                        aveq(istate)%q%improper,aveq(istate)%qx%el,aveq(istate)%qx%vdw, &
                        aveq(istate)%qq%el,aveq(istate)%qq%vdw,aveq(istate)%qp%el,&
                        aveq(istate)%qp%vdw, &
                        aveq(istate)%qw%el,aveq(istate)%qw%vdw,aveq(istate)%restraint

                end do
17      format(a,t23,i2,1x,i6,f9.6,14f8.2)
        end do  !ifile

        if(nfiles > 1) then !the following is meaningless for a single file
                dgf=0.
                dgfsum=0.
                sum=0.
                veff1=0.
                veff2=0.
                dv=0.
                do ifile=1,nfiles-1
                        !check that number of points > number to skip
                        if(nskip >= fep(ifile)%npts) then
                                write(*,999) ifile, nskip, fep(ifile)%npts
999                             format('file',i5,' contains only',i5,' points. can''t skip',i5)
                        end if
                        do ipt=nskip+1,fep(ifile)%npts
                                do istate=1,nstates
                                        veff1=veff1+fep(ifile)%lambda(istate)*fep(ifile)%v(istate,ipt)
                                        veff2=veff2+fep(ifile+1)%lambda(istate)*fep(ifile)%v(istate,ipt)
                                end do
                                dv=veff2-veff1
                                veff1=0.
                                veff2=0.
                                sum=sum+exp(-dv/rt)
                        end do
                        sum=sum/real(fep(ifile)%npts-nskip)
                        dgf(ifile)=-rt*dlog(sum)
                        dgfsum(ifile+1)=dgfsum(ifile)+dgf(ifile)
                        sum=0.
                end do
                dgrsum=0.
                dgr=0.
                sum=0.
                veff1=0.
                veff2=0.
                dv=0.
                do ifile=nfiles,2,-1
                        do ipt=nskip+1,fep(ifile)%npts
                                do istate=1,nstates
                                veff1=veff1+fep(ifile)%lambda(istate)*fep(ifile)%v(istate,ipt)
                                veff2=veff2+fep(ifile-1)%lambda(istate)*fep(ifile)%v(istate,ipt)
                        end do
                                dv=veff2-veff1
                                veff1=0.
                                veff2=0.
                                sum=sum+exp(-dv/rt)
                        end do
                        sum=sum/real(fep(ifile)%npts-nskip)
                        dgr(ifile)=-rt*dlog(sum)
                        dgrsum(ifile-1)=dgrsum(ifile)+dgr(ifile)
                        sum=0.
                end do
                write(*,*) 
                write(*,*) 
                write(*,21)
                write(*,22)
21              format('# part 1: free energy perturbation summary:')
22              format('# lambda(1)      dgf sum(dgf)      dgr sum(dgr)     <dg>')
        dg(1)=0.0
        do ifile=2,nfiles
            dg(ifile)=dg(ifile-1)+0.5*(dgf(ifile-1)-dgr(ifile))
        end do
        do ifile=1,nfiles
           write (*,23) &
               fep(ifile)%lambda(1),dgf(ifile-1),dgfsum(ifile), &
                             dgr(ifile+1),dgrsum(ifile),dg(ifile)
        end do
23              format(2x,f9.6,5f9.3)

                write (*,*)
                write (*,'(a,f9.2)') '# min energy-gap is: ',gapmin
                write (*,'(a,f9.2)') '# max energy-gap is: ',gapmax

                !-----------------------------------
                ! reaction free energy is calculated
                !-----------------------------------

        write (*,*)
            write (*,*)
                write(*,24)
                write(*,25)
24              format('# part 2: reaction free energy summary:')
25              format('# lambda(1)  bin energy gap      dga     dgb     dgg    # pts    c1**2    c2**2')
26              format(2x,f9.6,i5,2x,4f9.2,2x,i5,2f9.3)
                gaprange=gapmax-gapmin      !range of reaction coordinate
                xint=gaprange/real(nbins)   !divide r.c. into bins
                do ifile=1,nfiles         
                        avdvv=0.
                        avdvg=0.
                        sumv=0.
                        sumg=0.
                    avc1=0.
            avc2=0.
                        avr=0.
                        dvv=0.
                        dvg=0.
                        nbinpts=0

                        do ipt=nskip+1,fep(ifile)%npts
                                ibin=int((fep(ifile)%gap(ipt)-gapmin)/xint)+1
                                veff=0.
                                do istate=1,nstates
                                        veff=veff+fep(ifile)%lambda(istate)*fep(ifile)%v(istate,ipt)
                                end do  !states
                                dvv(1:nstates)=fep(ifile)%v(1:nstates,ipt)-veff
                                dvg=fep(ifile)%vg(ipt)-veff
                                avdvv(:,ibin)=avdvv(:,ibin)+dvv(:)
                                avdvg(ibin)=avdvg(ibin)+dvg
                        avc1(ibin)=avc1(ibin)+fep(ifile)%c1(ipt)
                avc2(ibin)=avc2(ibin)+fep(ifile)%c2(ipt)
                                !only gives first r_xy distance
                                if(nnoffd > 0)  avr(ibin)=avr(ibin)+fep(ifile)%r(1,ipt)          
                                nbinpts(ibin)=nbinpts(ibin)+1
                        end do          !ipt
                        do ibin=1,nbins
                                if ( nbinpts(ibin) .ne. 0 ) then
                    avc1(ibin)=avc1(ibin)/real(nbinpts(ibin))
                    avc2(ibin)=avc2(ibin)/real(nbinpts(ibin))
                                        avr(ibin)=avr(ibin)/real(nbinpts(ibin))                         
                                        avdvv(:,ibin)=avdvv(:,ibin)/nbinpts(ibin)
                                        avdvg(ibin)=avdvg(ibin)/nbinpts(ibin)

                                end if 
                        end do !ibin
                        do ipt=nskip+1,fep(ifile)%npts
                                ibin=int((fep(ifile)%gap(ipt)-gapmin)/xint)+1
                                veff=0.
                                do istate=1,nstates
                                        veff=veff+fep(ifile)%lambda(istate)*fep(ifile)%v(istate,ipt)
                                end do  !istate

                                do istate=1,nstates
                                        dvv(istate)=fep(ifile)%v(istate,ipt)-veff-avdvv(istate,ibin)
                                end do
                                dvg=fep(ifile)%vg(ipt)-veff-avdvg(ibin)
                                sumv(:,ibin)=sumv(:,ibin)+exp(-dvv(:)/rt)
                                sumg(ibin)=sumg(ibin)+exp(-dvg/rt)
                        end do   !ipt

                        do ibin=1,nbins
                                if (nbinpts(ibin).ge.nmin) then
                                    binsum(ibin,2)=binsum(ibin,2)+avc1(ibin)*nbinpts(ibin)
                                        binsum(ibin,3)=binsum(ibin,3)+avc2(ibin)*nbinpts(ibin)
                                        binsum(ibin,4)=binsum(ibin,4)+avr(ibin)*nbinpts(ibin)  !bin-averaged r_xy
                                        sumv(:,ibin)=sumv(:,ibin)/real(nbinpts(ibin)) 
                                        sumg(ibin)=sumg(ibin)/real(nbinpts(ibin))

                                ptsum(ibin)=ptsum(ibin)+nbinpts(ibin)

                                        do istate=1,nstates
                                                sumv2(istate,ibin)=-rt*dlog(sumv(istate,ibin))+avdvv(istate,ibin)
                                        end do
                                        sumg2(ibin)=-rt*dlog(sumg(ibin))+avdvg(ibin) 
                                        ! these are the diabatic free energy curves
                                        dgv(:)=dg(ifile)+sumv2(:,ibin) 
                                        ! this is the reaction free energy
                                        dgg=dg(ifile)+sumg2(ibin)

                                    binsum(ibin,1)=binsum(ibin,1)+dgg*int(nbinpts(ibin))

                                        write (*,26) fep(ifile)%lambda(1),ibin, &
                                                gapmin+real(ibin)*xint-xint/2., dgv(1),dgv(2),dgg, &
                                                int(nbinpts(ibin)),avc1(ibin),avc2(ibin) 
                                end if
                        end do  !ibin
                end do      !ifile
write(*,*)
write(*,*)
write(*,27)
write(*,28)

27              format('# part 3: bin-averaged summary:')
28              format('# bin  energy gap  <dgg> <dgg norm> pts  <c1**2> <c2**2> <r_xy>')

        do ibin=1,nbins
        if (ptsum(ibin).ge.nmin) then
        binsum(ibin,1)=binsum(ibin,1)/real(ptsum(ibin)) ! bin-averaged reaction free energy
        binsum(ibin,2)=binsum(ibin,2)/real(ptsum(ibin)) ! bin-averaged c1**2
        binsum(ibin,3)=binsum(ibin,3)/real(ptsum(ibin)) ! bin-averaged c2**2
        binsum(ibin,4)=binsum(ibin,4)/real(ptsum(ibin)) ! bin-averaged r_xy
        end if
        end do
        min=minval(binsum(:,1))
        do ibin=1,nbins
                if (ptsum(ibin).ge.nmin) then
 29             format(i4,1x,3f9.2,2x,i5,3f8.3,4f8.2)
                write(*,29) ibin,gapmin+real(ibin)*xint-xint/2.,binsum(ibin,1),  &
                binsum(ibin,1)-min,int(ptsum(ibin)),binsum(ibin,2),binsum(ibin,3),binsum(ibin,4)
        end if
        end do !ibin
        end if !nfiles >1

        !clean up
        do ifile=1,nfiles
                deallocate(fep(ifile)%v, fep(ifile)%vg, fep(ifile)%gap, &
                fep(ifile)%c1, fep(ifile)%c2)
                if(nnoffd > 0) deallocate(fep(ifile)%r)
        end do
        deallocate(fep)

        deallocate(hij,d,e,stat=err)
!.......................................................................

contains

!------------------------------

subroutine prompt (outtxt)
        character(*) outtxt
#if defined (__osf__)
        !prompt to stderr using unit 5=stdin on osf/1=dec unix
        integer, parameter                      ::      f=5
        !write (f,'($,a)') outtxt
#elseif defined (_win32)
        !open the err file on win32
        integer, save :: f
        if(f==0) then
                f=17
                open(unit=f,file='err')
        end if
        !write (f1,'($,a)') outtxt
#else
        !otherwise prompt to stdout
        integer, parameter                      ::      f=6
        !write (f2,'($,a)') outtxt
#endif
write (f,'(a,$)') outtxt             
end subroutine prompt

!------------------------------

subroutine tred2(a,n,np,d,e)
!------------------------------------------------------------
! this subroutine reduces a symmetric matrix to tridiagonal
! form. the tridiagonal matrix can further be diagonalized by
! the subroutine tqli.
! these subroutines were copied by karin kolmodin 20 nov. 1997
! from http://rsc.anu.au/hws/courses/mathmeth/node70.html
! and rewritten in f90.
!------------------------------------------------------------ 
          real(8),dimension(:)      ::d,e
          real(8),dimension(:,:)    ::a
      integer                   ::np,i,j,k,l,n,err
          real(8)                   ::scale,f,g,h,hh
           
      if(n.gt.1)then
        do i=n,2,-1  
           l=i-1
           h=0.
           scale=0.
          if(l.gt.1)then
            do k=1,l
              scale=scale+abs(a(i,k))
            end do
            if(scale.eq.0.)then
              e(i)=a(i,l)
            else
              do k=1,l
                a(i,k)=a(i,k)/scale
                h=h+a(i,k)**2
              end do
              f=a(i,l)
              g=-sign(sqrt(h),f)
              e(i)=scale*g
              h=h-f*g
              a(i,l)=f-g
              f=0.
              do j=1,l
                a(j,i)=a(i,j)/h
                g=0.
                do k=1,j
                  g=g+a(j,k)*a(i,k)
                end do
                if(l.gt.j)then
                  do k=j+1,l
                    g=g+a(k,j)*a(i,k)
                  end do
                end if
                e(j)=g/h
                f=f+e(j)*a(i,j)
              end do
              hh=f/(h+h)
              do j=1,l
                f=a(i,j)
                g=e(j)-hh*f
                e(j)=g
                do k=1,j
                  a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
                end do
              end do
            end if
          else
            e(i)=a(i,l)
          end if
          d(i)=h
        end do
      end if
      d(1)=0.
      e(1)=0.
      do i=1,n
        l=i-1
        if(d(i).ne.0.)then
          do j=1,l
            g=0.
            do k=1,l
              g=g+a(i,k)*a(k,j)
            end do
            do k=1,l
              a(k,j)=a(k,j)-g*a(k,i)
            end do
          end do
        end if
        d(i)=a(i,i)
        a(i,i)=1.
        if(l.ge.1)then
          do j=1,l
            a(i,j)=0.
            a(j,i)=0.
          end do
        end if
      end do
end subroutine tred2 

!-----------------------------------

subroutine tqli(d,e,n,np,z)
!------------------------------------------------------------ 
! this subroutine diagonalizes a tridiagonal matrix which has 
! been prepared by the subroutine tred2.
! these subroutines were copied by karin kolmodin 20 nov. 1997
! from http://rsc.anu.au/hws/courses/mathmeth/node70.html
! and rewritten in f90
!------------------------------------------------------------  
   implicit none
      real(8),dimension(:)               ::d,e
          real(8),dimension(:,:)             ::z
          integer                            ::i,n,np,k,l,m,iter
          real(8)                            ::dd,g,r,s,c,p,f,b

        if (n.gt.1) then
        do i=2,n
          e(i-1)=e(i)
        end do
        e(n)=0.
        do l=1,n
          iter=0
1         do m=l,n-1
            dd=abs(d(m))+abs(d(m+1))
            if (abs(e(m))+dd.eq.dd) go to 2
          end do
          m=n
2         if(m.ne.l)then
            if(iter.eq.30)stop 'too many iterations'
            iter=iter+1
            g=(d(l+1)-d(l))/(2.*e(l))
            r=sqrt(g**2+1.)
            g=d(m)-d(l)+e(l)/(g+sign(r,g))
            s=1.
            c=1.
            p=0.
            do i=m-1,l,-1
              f=s*e(i)
              b=c*e(i)
              if(abs(f).ge.abs(g))then
                c=g/f
                r=sqrt(c**2+1.)
                e(i+1)=f*r
                s=1./r
                c=c*s
              else
                s=f/g
                r=sqrt(s**2+1.)
                e(i+1)=g*r
                c=1./r  
                s=s*c
              end if
              g=d(i+1)-p
              r=(d(i)-g)*s+2.*c*b
              p=s*r
              d(i+1)=g+p
              g=c*r-b
              do k=1,n
                f=z(k,i+1)
                z(k,i+1)=s*z(k,i)+c*f
                z(k,i)=c*z(k,i)-s*f
              end do
            end do
            d(l)=d(l)-p
            e(l)=g
            e(m)=0.
            go to 1
          end if
        end do
      end if
      return
end subroutine tqli

!-----------------------------------

end program qfep
