!_____________________
        MODULE configuration
!_____________________
	IMPLICIT NONE; SAVE

        INTEGER   :: spin_or_boson       ! system type
        INTEGER   :: onsite_nn 
        INTEGER   :: close_open          
	INTEGER   :: kink_max            ! number of intervals
        INTEGER   :: dim, dir            ! dimension-directions
        INTEGER   :: triangular          ! triangular lattice switch
        INTEGER, allocatable :: Ns(:)    ! system sizes
        INTEGER   :: Nsite               ! Nsite = prod N(i)         
        DOUBLE PRECISION     :: inv_temp ! inverse temperature

!Hilbert space and correlators           	
        INTEGER    ::	Nmax               ! number of tau points in G 
        INTEGER    ::   fl_max             ! max occupation number        
 	DOUBLE PRECISION, allocatable  :: gr(:,:)   	! G-function
        DOUBLE PRECISION, allocatable  :: gr0(:)        ! p=0 G-function
 	DOUBLE PRECISION, allocatable  :: grnm(:,:,:)   ! G_{nm}-function 
	DOUBLE PRECISION, allocatable  :: timegr(:), timeq(:)
   	DOUBLE PRECISION, allocatable  :: corr_nn(:,:) ! density-density corr.	
   	DOUBLE PRECISION, allocatable  :: ttt(:)
   	DOUBLE PRECISION                        :: jc
	                                ! times where G and dens. -dens.
	                                ! correlators are calculateded  
      INTEGER, PARAMETER :: NJJ=100
	DOUBLE PRECISION   :: zjj, jjcorr(0:NJJ)
	INTEGER :: kwcheck

        DOUBLE PRECISION :: dtgr,dtq,twin,twinq                                 
	DOUBLE PRECISION :: tcc_win 
	DOUBLE PRECISION :: tc_win
	                                ! windows to collect statistics 
	                                ! for G and dens.-dens.
	                      
!world-line diagram description and linking          
	INTEGER, allocatable          :: fl(:)     
	                               ! interval state/color/occ.number
	INTEGER, allocatable          :: site(:)   
	                               ! interval-site correspondence
	DOUBLE PRECISION, allocatable :: tau(:), time(:)
                                       ! duration in time 
                                       ! and absolute time
	INTEGER, allocatable          :: ml(:), mr(:), mk(:) 
                                       ! describes left/right associations
                                       ! in time and through kink
	INTEGER, allocatable          :: ms(:,:)
	                               ! describes same-time associations  
	                               ! between intervals in real space      
	INTEGER, allocatable          :: back(:) 
	                               !inverse directions
	                               
! Managing the associations - linked list
        INTEGER                 :: nmnm
	INTEGER, allocatable   	:: nlist(:)
	INTEGER, allocatable   	:: numnam(:) 
 	                   
! worm parameters
        INTEGER              :: ira, masha  ! worm end points
        INTEGER, allocatable :: Rim(:)      ! ira-masha real space distance
        LOGICAL              :: present     ! yes/no for the worm presence 
        DOUBLE PRECISION :: tim             !ira-masha im. time distance
 
! Random number generator variables  
      DOUBLE PRECISION :: UGEN(97), C7, CD, CM 
      integer          :: ij,kl,j,I97, J97 
      double precision :: x                    ! Particular random number 
     
!input parameters from the Hamiltonian and their code conversions

        INTEGER          :: Mz      !  magnetization/particle number 
        DOUBLE PRECISION :: Mzmark  !  about canonical fluctuations
        LOGICAL          :: CANONICAL ! TRUE if canonical 
        DOUBLE PRECISION :: onsite   ! on_site interaction
        DOUBLE PRECISION :: z_coupling1
                          ! nearest-neighbor int. along chains
        DOUBLE PRECISION :: z_coupling2
                          ! nearest-neighbor int. perp. to chains
        DOUBLE PRECISION :: xy_coupling1
                          ! nearest-neighbor flip-flop amplt. 
                          ! (hopping) along chains
        DOUBLE PRECISION :: xy_coupling2
                          ! nearest-neighbor flip-flop amplt. 
                          ! (hopping) perp to chains  	
        DOUBLE PRECISION :: mag_field, random_field(2),mag0
                          ! magnetic field/chemical potential 
                          ! and random mag. field/chem. pot.
        DOUBLE PRECISION :: peierls_z,peierls_xy
                          ! Peierls parameter in chains: dJz 
        DOUBLE PRECISION :: hole_conc
                          ! concentration of (static here) holes (if any)

! all of the above converted to: 
!                    new_value_time = old_value_time * t_quanta
!                    new_value_energy =old_value_energy / t_quanta
! where t_quanta is time-granulation to ensure that im. time 
! sum rules are strictly enforced 

	DOUBLE PRECISION :: t_quanta, beta, beta2, U0, Jz1,Jz2
	DOUBLE PRECISION :: source          
	DOUBLE PRECISION, allocatable :: Int_z(:,:),Int_xy(:,:)
	DOUBLE PRECISION :: s_hop1,s_hop2,H_0, W(2)
	DOUBLE PRECISION, allocatable :: Hsite(:)	 
	DOUBLE PRECISION :: dJz, dJxy
	DOUBLE PRECISION :: Esp
	DOUBLE PRECISION :: wmso
	DOUBLE PRECISION, allocatable :: Fss(:), Fs2(:)
	integer, allocatable           :: vector(:,:) ! vector for the site 
! tabulated values of potential energies for different interval states
        DOUBLE PRECISION, allocatable :: w_hop(:,:),s_h(:,:) 
! tabulated values of hopping/exchange terms for different interval states         		
        DOUBLE PRECISION, allocatable :: Pat(:)
! particle number/spin projection vs interval state correspondence  
  

!MC management: probabilities of procedures  
	        
	DOUBLE PRECISION, PARAMETER :: w_del=0.12d0, w_jump=0.07d0
        DOUBLE PRECISION, PARAMETER :: w_move=0.1d0, w_cre=1.d0       
! del appears 2 times, jump appears 8 times, move appears 2 times    
        DOUBLE PRECISION, PARAMETER :: w1=w_del, w2=w1+w_del
        DOUBLE PRECISION, PARAMETER :: w3=w2+w_jump,w4=w3+w_jump
        DOUBLE PRECISION, PARAMETER :: w5=w4+w_jump,w6=w5+w_jump 
        DOUBLE PRECISION, PARAMETER :: w7=w6+w_jump,w8=w7+w_jump 
        DOUBLE PRECISION, PARAMETER :: w9=w8+w_jump,w10=w9+w_jump
        DOUBLE PRECISION, PARAMETER :: w11=w10+w_move,w12=w11+w_move
         
! counters, estimators, etc         
         
        INTEGER  :: Mwind , Mpat      ! MAX winding number and 
                                      ! magnetization expected        
	INTEGER :: km, kmm, dfl , kinks
	DOUBLE PRECISION :: p_ther
	INTEGER          :: p_mes, p_prn, p_wr
	                  ! thermolization, measuring, printing and writing 
	                  ! in the number of updates
        INTEGER :: aa(0:30),bb(0:30) 
	DOUBLE PRECISION :: prtcls, Udiag
	DOUBLE PRECISION :: Upot, Ekin , Energy   ! energies
	DOUBLE PRECISION :: step                  ! MC update number
	DOUBLE PRECISION :: stat                  ! partition function
        DOUBLE PRECISION :: Sq, StatSq            ! static structure factor
        DOUBLE PRECISION :: Sqe, StatSqe          ! equal_time structure factor
	DOUBLE PRECISION, allocatable :: Nav(:)     ! particles
	DOUBLE PRECISION :: nsquare, stat2
        DOUBLE PRECISION, allocatable ::  ccS(:), ssS(:) !tabulated cos and sin

	                  ! on-site density square and its statistics   
	DOUBLE PRECISION, allocatable :: arr(:)       
	DOUBLE PRECISION, allocatable :: Wind(:,:) 
	                               ! winding numbers

        INTEGER, PARAMETER :: prnt_max=100000
        INTEGER :: prntout
        DOUBLE PRECISION :: amin, amax, tmin, tmax 
        DOUBLE PRECISION :: Sq_prn(prnt_max)
        DOUBLE PRECISION :: rs_prn(prnt_max)
        DOUBLE PRECISION :: E_prn(prnt_max) 
        DOUBLE PRECISION :: N_prn(prnt_max)     
        DOUBLE PRECISION :: Sqe_prn(prnt_max)
        DOUBLE PRECISION :: UD_prn(prnt_max),UpDown
        DOUBLE PRECISION :: histUD(0:100), binN(0:100)
        DOUBLE PRECISION, PARAMETER :: UDbin=0.03d0
	                               
! new/continue switches 
	INTEGER :: old_new1, old_new2  
  
! some useful constants   
        LOGICAL, PARAMETER :: plus=.TRUE., minus=.FALSE.
	DOUBLE PRECISION, PARAMETER :: pi=3.141592653589793d0
        DOUBLE PRECISION, PARAMETER :: eilera=0.5772156649d0
	DOUBLE PRECISION, PARAMETER :: katalana=0.9159655942d0
	DOUBLE PRECISION, PARAMETER :: e_number=2.718281828459d0
	DOUBLE PRECISION, PARAMETER :: pi_d2=pi/2.d0, pi_m2=pi*2.d0
	DOUBLE PRECISION, PARAMETER :: pi_m4=pi*4.d0,vol=pi_m2**3
        DOUBLE PRECISION, PARAMETER :: piSq=pi*4.d0/3.d0
	DOUBLE PRECISION, PARAMETER :: R0=1.d-12, R1=1.d0-1.d-12  
  
! variables in all updates   
        INTEGER          :: here  ! updated site 
        DOUBLE PRECISION :: pot   ! pot energy change
	DOUBLE PRECISION :: p_a   ! Acceptance ratio 
	DOUBLE PRECISION :: thop  ! hopping matrix element   
  
                             END MODULE configuration
!********************************************************
        USE configuration 
	IMPLICIT NONE; SAVE 
        integer :: i,i1,id
        DOUBLE PRECISION :: ist,TotalSteps
        REAL*4 :: result(11)
        logical :: alive
      id=5
      OPEN(1,file='parwv2.dat')
      READ(id,*)  spin_or_boson ; print*, 'spin/boson  ', spin_or_boson 
      READ(id,*)  onsite_nn     ; print*, 'onsite/nn   ', onsite_nn
      READ(id,*)  dim           ; print*, 'dim         ', dim
      READ(id,*)  triangular    ; print*, 'lattice type', triangular 
      allocate(Ns(dim))
        IF(triangular==0) then
            dir=2*dim              ! square lattice 
        else
            dir=2*dim+2            ! triangular lattice 
        ENDIF
                                    print*, 'directions', dir
         Nsite=1
         do i=1,dim 
      READ(id,*)  Ns(i)           ; print*, 'i,Ns(i)      ', i,Ns(i)
         Nsite=Nsite*Ns(i)
         end do 

      READ(id,*)  Nmax		     ; print*, 'Nmax        ', Nmax

	allocate(vector(2,Nsite)) 
      allocate(gr(0:Ns(1)-1,0:2*Nmax),gr0(0:2*Nmax))
      allocate(timegr(0:2*Nmax),timeq(0:2*Nmax))
      allocate(corr_nn(0:Ns(1)-1,0:6*Nmax))
      allocate(ttt(0:6*Nmax+1))      
      READ(id,*)  fl_max          ; print*, 'fl_max      ', fl_max
      allocate(grnm(0:Ns(1),0:fl_max,0:fl_max))              
      READ(id,*)  kink_max        ; print*, 'kink_max    ', kink_max
      allocate(fl(kink_max),site(kink_max))
      allocate(tau(kink_max),time(kink_max))  
      allocate(ml(kink_max),mr(kink_max),mk(kink_max))
      allocate(ms(dir,kink_max),back(dir),Rim(dim))
      allocate(nlist(kink_max),numnam(kink_max))  
      allocate(Fss(-3:fl_max),Fs2(-3:fl_max))
      allocate(w_hop(-3:fl_max,-3:fl_max))
      allocate(Pat(-3:fl_max),s_h(-3:fl_max,-3:fl_max))
      allocate(Int_z(Nsite,dir),Int_xy(Nsite,dir))
      allocate(Hsite(Nsite)) 	
      allocate(ccS(Nsite),ssS(Nsite))
      READ(id,*)  Mz              ; print*, 'Mz          ', Mz
      READ(id,*)  Mzmark          ; print*, 'Mzmark      ',Mzmark
      READ(id,*)  CANONICAL       
      READ(id,*)  onsite          ; print*, 'onsite      ', onsite
      READ(id,*)  z_coupling1 	 ; print*, 'Jz1         ', z_coupling1
      READ(id,*)  xy_coupling1	 ; print*, 'Jx1         ', xy_coupling1
      READ(id,*)  z_coupling2 	 ; print*, 'Jz2         ', z_coupling2
      READ(id,*)  xy_coupling2	 ; print*, 'Jx2         ', xy_coupling2      
      READ(id,*)  peierls_z       ; print*, 'Peierls_z   ', peierls_z
      READ(id,*)  peierls_xy      ; print*, 'Peierls_xy  ', peierls_xy     
      READ(id,*)  hole_conc       ; print*, 'Hole conc.  ', hole_conc
      READ(id,*)  mag_field	 ; print*, 'H           ', mag_field
      READ(id,*)  random_field(1); print*, 'W(1)   ', random_field(1)
	READ(id,*)  random_field(2); print*, 'W(2)   ', random_field(2)
      READ(id,*)  close_open      ; print*, 'close/open  ', close_open
      READ(id,*)  t_quanta        ; print*, 't_quanta    ', t_quanta 
      READ(id,*)  source          ; print*, 'source      ', source 
      READ(id,*)  Mwind           ; print*, 'Mwind       ', Mwind           
      READ(id,*)  Mpat            ; print*, 'Mpat        ', Mpat    
      allocate(Nav(-Mpat:Mpat+1))
      allocate(arr(dir))
      allocate(Wind(dir/2,-Mwind:Mwind))
      READ(id,*)  twin    	; print*, 'G-window    ', twin 
      READ(id,*)  tcc_win         ; print*, '<nn>-window ', tcc_win
      READ(id,*)  jc              ; print*, '<nn> freq.  ', jc
      READ(id,*)  p_ther  ; print*, 'therm. time in updates ',p_ther
      READ(id,*)  p_mes   ; print*, 'measure after # updates',p_mes
      READ(id,*)  p_prn   ; print*, 'print every # updates  ',p_prn
      READ(id,*)  p_wr    ; print*, 'write every # updates  ',p_wr
      READ(id,*)  old_new1
      READ(id,*)  old_new2
      READ(id,*)  TotalSteps      ; print*, 'TotalSteps  ', TotalSteps
      READ(id,*)  inv_temp		 ; print*, 'beta        ', inv_temp
      READ(id,*)  ij              ; print*, 'seed1       ', ij 
      READ(id,*)  kl              ; print*, 'seed2       ', kl

                	CLOSE(1)
                	      				    
      mag0=mag_field
      IF(spin_or_boson==1) then
      onsite=onsite/2.d0
      mag_field=mag_field+onsite
      endif

	source=source/(inv_temp*Nsite); 
	source=source/t_quanta; source=source/t_quanta
	beta=inv_temp*t_quanta; beta=AINT(beta)
	beta2=beta/2.d0
    
	U0=onsite/t_quanta	
        Jz1=z_coupling1/t_quanta
        Jz2=z_coupling2/t_quanta 
        dJz=peierls_z/t_quanta
        dJxy=peierls_xy*0.5d0/t_quanta 
	s_hop1=xy_coupling1*0.5d0/t_quanta
        s_hop2=xy_coupling2*0.5d0/t_quanta
	H_0=mag_field/t_quanta
	W=random_field/t_quanta
	
	print*, 'parameters done'      	

       

	print*, 'sum of probabilities =', w12 
!	print*, w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12; pause
!       Initializing random number generator             
!        ij = 4836;  kl = 2748; 
        call SRANMAR(ij,kl)    
        call hamiltonian;  call initial; call testing
        print*, 'START MONTE CARLO'

       open(100, file='flag',status='replace') 
       close(100)
                      
	ee: do;  step=step+1.d0; 
	i1=step/1.d9; ist=step-i1*1.d9; i1=ist
        inquire(file='flag',exist=alive) 
        IF((TotalSteps>0.001 .and. step/1.d6>TotalSteps)
     x  .or. (.not.alive)) exit

	call manager; 
        IF(step>p_ther) then;  
!        IF(p_ther.ne.0) then; call WR_cnf; 
!        old_new1=1; call initial; endif 
        p_ther=0 
        
        IF(mod(i1,p_mes)==0) call measure(i1); endif
!	IF(mod(i1,p_wr)==0) print*, i1/p_wr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        IF(mod(i1,p_prn)==0) call printing(i1)
!        IF(mod(i1,p_wr)==0)  then; 
!	                     Print*, '!!!!! WRITE !!!!!'
!                             call WR_cnf; call WR_stat
!                             PRINT*, '!!!!! DONE !!!!!'; endif


!         IF(mod(i1,1)==0) call testing

         enddo ee      	
	call printing(i1)
	call WR_cnf;call WR_stat

      
      CONTAINS          

!_________MC management_________________________________________ 

          SUBROUTINE manager         
	  INTEGER :: name, next,nextl
          DOUBLE PRECISION :: tck

    	  x=rndm()
                		
          IF(.not.present) THEN; km=1; dfl=2*x; dfl=2*dfl-1; 
          IF(ABS(dfl).ne.1)stop 'dfl'; call create_worm(dfl)
	  ELSE IF(x<w1)  THEN; km=2;   call delete_worm(ira)
	  ELSE IF(x<w2)  THEN; km=3;   call delete_worm(masha)	  
	  ELSE IF(x<w3)  THEN; km=4;   call insert_right(ira)
	  ELSE IF(x<w4)  THEN; km=5;   call insert_right(masha)
 	  ELSE IF(x<w5)  THEN; km=6;   call insert_left(ira)
	  ELSE IF(x<w6)  THEN; km=7;   call insert_left(masha)
 	  ELSE IF(x<w7)  THEN; km=8;   call delete_right(ira)
 	  ELSE IF(x<w8)  THEN; km=9;   call delete_right(masha)
 	  ELSE IF(x<w9)  THEN; km=10;  call delete_left(ira)
	  ELSE IF(x<w10) THEN; km=11;  call delete_left(masha)
	  ELSE IF(x<w11) THEN; km=12;  call move_worm(ira)
	  ELSE IF(x<w12) THEN; km=13;  call move_worm(masha)
	  ELSE; print*, 'probabilities do not sum to 1-W', w12; stop
	  END IF

!________________________________________________________________________________

         tck=prtcls/beta
         IF(CANONICAL .AND. ABS(tck-Nsite) > Mzmark) stop 'tck failed'
!_____________________________________________________________ CANONICAL !!!!!!!    
	
!	   call testing
     	   END SUBROUTINE manager
     	         	 
!______ reading the saved configuration_______________
	SUBROUTINE RD_cnf
	integer :: name,d
	
	open(1,file='cnf.dat');  
	read(1,*) nmnm,kinks,present,ira,masha,Rim,tim
 	read(1,*) prtcls,Udiag
 	read(1,*) arr
 	DO j=1, nmnm;  
	read(1,*) nlist(j); name=nlist(j); numnam(name)=j
	read(1,*)ml(name),mr(name),mk(name),(ms(d,name), d=1,dir )	
	read(1,*)fl(name),site(name),tau(name),time(name)
 	END DO 
	   DO j=nmnm+1,kink_max; read(1,*) nlist(j); 
	                   numnam(nlist(j))=j; ENDDO
 	read(1,*) ugen 
        read(1,*) c7,cd,cm 
        read(1,*) i97,j97      
        close(1)

!        call testing 
        print*, 'configuration reading done'
	END SUBROUTINE RD_cnf
	    
!______ writing the configuration in file______________
	SUBROUTINE WR_cnf
	integer :: name,d
	
	open(1,file='cnf.dat');  
	write(1,*) nmnm,kinks,present,ira,masha,Rim,tim
	write(1,*) prtcls,Udiag
	write(1,*) arr
	DO j=1, nmnm; 
	name=nlist(j); write(1,*) name
	write(1,*)ml(name),mr(name),mk(name),(ms(d,name), d=1,dir )	
	write(1,*)fl(name),site(name),tau(name),time(name)
	END DO 
	   DO j=nmnm+1,kink_max; write(1,*) nlist(j); ENDDO	
	write(1,*) ugen 
        write(1,*) c7,cd,cm 
        write(1,*) i97,j97     
	close(1)
	    
	END SUBROUTINE WR_cnf 

!_______________________________________________________ 
       SUBROUTINE WR_stat
       integer :: kp
       DOUBLE PRECISION :: xx 
         
      do kp=0,Ns(1)-1
      do j=0,2*Nmax-1
      corr_nn(kp,j)=corr_nn(kp,j)+corr_nn(kp,j+2*Nmax)
     .+corr_nn(kp,j+4*Nmax)
      corr_nn(kp,j+2*Nmax)=0.d0; corr_nn(kp,j+4*Nmax)=0.d0
      enddo
      enddo   
                  
       open(1, file='stat.dat') 
       write(1,*) step, stat, Upot,Ekin,Energy
	 write(1,*) zjj, jjcorr
       write(1,*) StatSq, Sq 
       write(1,*) StatSqe, Sqe, UpDown
       write(1,*) Nav      
       write(1,*) Wind
       write(1,*) timegr, gr, grnm, corr_nn, nsquare,stat2       
       write(1,*) prntout
       do kp=1,prntout
       write(1,*) kp,Sq_prn(kp),rs_prn(kp),UD_prn(kp)
       write(1,*) E_prn(kp),N_prn(kp),Sqe_prn(kp)
       enddo
       write(1,*) histUD
       close(1)
	    
       open(1, file='histUD.dat')
       do j=1,100
       write(1,*) (j-0.5d0)*UDbin, histUD(j)/(histUD(0)*UDbin)
       enddo
       close(1)


       open(1, file='green.dat')
       xx=1.d0*Nav(Mpat+1)/stat/Nsite
       xx=xx/gr(0,2*Nmax)
       do j=0,Ns(1)-1
       write(1,*) j, gr(j,0)*xx, gr(j,2*Nmax)*xx
       enddo         
       close(1)

       open(1, file='green1.dat')
       do j=0,2*Nmax
       write(1,*)timegr(j),gr(0,j)*xx,gr(1,j)*xx
       enddo
       close(1)


!       open(1, file='green.dat')
!       write(1,*) Nmax
!       write(1,*) timegr, gr, 1.d0*Nav(Mpat+1)/stat/Nsite
!       close(1)

       open(1, file='gr0.dat')
       do j=0,2*Nmax
       write(1,*)timegr(j),gr0(j)/(gr0(0)+1.d-16)
       enddo
       close(1)

   
       open(1, file='grnm.dat') 
       write(1,*) Ns(1),fl_max, grnm,nsquare, stat2
       write(1,*) 1.d0*Nav(Mpat+1)/stat/Nsite       
       close(1)       
       
        open(1, file='corr.dat') 
        write(1,*) Nmax
        write(1,*) timegr, corr_nn, nsquare, stat2       
        close(1)
 
       open(1, file='counters.dat')      
       DO j=1,14 
       p_a=bb(j)/(aa(j)+1.d-5) 
       write(1,*) j, aa(j),bb(j),p_a 
       ENDDO 
       close(1) 
       aa=0; bb=0
       
       xx=1.d-14;do j=-Mpat,Mpat; xx=xx+Nav(j); enddo
       open(1,file='pat.dat') 
       do j=-Mpat,Mpat
       write(1,*) j, Nav(j)/(xx+1.d-14)
       enddo     
       close(1) 
      
       open(1,file='wind.dat')
       do kp=1,dim; 
       xx=1.d-14; do j=-Mwind,Mwind; xx=xx+Wind(kp,j); enddo        
       do j=-Mwind,Mwind
       write(1,*) j, Wind(kp,j)/xx 
       enddo     
       enddo
       close(1) 
       
!       call testing 
       END SUBROUTINE WR_stat 

!______reading previous statistics ________________
       SUBROUTINE RD_stat  
       integer :: kp, i
     
       open(1, file='stat.dat') 
       read(1,*) step, stat, Upot,Ekin,Energy 
	 read(1,*) zjj, jjcorr
       read(1,*) StatSq, Sq
       read(1,*) StatSqe, Sqe,UpDown
       read(1,*) Nav      
       read(1,*) Wind 
       read(1,*) timegr,gr,grnm, corr_nn, nsquare, stat2
       read(1,*) prntout
       do kp=1,prntout
       read(1,*) i,Sq_prn(kp),rs_prn(kp),UD_prn(kp)
       read(1,*) E_prn(kp),N_prn(kp),Sqe_prn(kp)
       enddo
       read(1,*) histUD
       close(1) 
        
       print*, 'statistics reading done' 
       END SUBROUTINE RD_stat 

!_____random number generator____________________________________
      DOUBLE PRECISION FUNCTION RNDM() 
      IMPLICIT NONE 
      DOUBLE PRECISION, parameter :: r1=1.d-15, r2=1.d-15      
      DOUBLE PRECISION :: RVAL  
 
      RVAL = UGEN(I97) - UGEN(J97) 
      IF ( RVAL < 0.D0 ) RVAL = RVAL + 1.0 
      UGEN(I97) = RVAL; I97=I97-1; IF (I97 == 0) I97 =97 
      J97=J97-1; IF (J97 == 0) J97=97 
      C7=C7-CD; IF (C7 < 0.D0 ) C7=C7+CM 
      RVAL = RVAL - C7 
      IF ( RVAL .LT. 0.D0 ) RVAL = RVAL + 1.0 
 
      RNDM = max(RVAL-r1,r2) 
      END FUNCTION RNDM 

!______rndm initialization_______________________________________                
       SUBROUTINE SRANMAR(IJ,KL) 
       IMPLICIT NONE 
       INTEGER IRAND 
       INTEGER :: i7, k7, l7, j7,ii7, jj7, m7, ij, kl 
       DOUBLE PRECISION :: S7, T7 
 
!      IF ( IJ<0 .OR. IJ>31328 .OR. KL<0 .OR. KL>30081) THEN
!        PRINT '(A)',' The first seed must be between 0 and 31328'
!        PRINT '(B)',' The second seed must be between 0 and 30081' 
!      STOP;  ENDIF 
 
      I7 = MOD(IJ/177, 177) + 2; J7 = MOD(IJ, 177) + 2 
      K7 = MOD(KL/169, 178) + 1; L7 = MOD(KL, 169) 
 
      DO  ii7 = 1, 97 
      S7 = 0.D0;  T7 = 0.5D0 
         DO  jj7 = 1, 24 
            M7 = MOD(MOD(I7*J7, 179)*K7, 179) 
            I7 = J7; J7 = K7; K7 = M7; L7 = MOD(53*L7+1, 169) 
            IF (MOD(L7*M7, 64) > 32) S7 = S7 + T7 
            T7 = 0.5D0 * T7 
         ENDDO 
      UGEN(ii7) = S7 
      ENDDO 
 
      C7 = 362436.D0 / 16777216.D0 
      CD = 7654321.D0 / 16777216.D0 
      CM = 16777213.D0 /16777216.D0 
 
      I97 = 97;  J97 = 33 
      RETURN   
      END subroutine SRANMAR
       
!________problem avoiding exponent______________________ 
         DOUBLE PRECISION FUNCTION EFUN(x) 
         IMPLICIT NONE 
         DOUBLE PRECISION, INTENT(IN) :: x 
         DOUBLE PRECISION, PARAMETER :: big=230.d0 
      
      IF (x>big) THEN;        EFUN=1.d100 
      ELSE IF (-x>big) THEN;  EFUN=1.d-100 
      ELSE;                 EFUN=exp(x) ;      END IF 
    		END FUNCTION EFUN 
    		
!________summing times mod-beta_________________________
	 DOUBLE PRECISION FUNCTION s_t(pm,ta,tb) 
         IMPLICIT NONE 
         LOGICAL, INTENT(IN) :: pm 
	 DOUBLE PRECISION, INTENT(IN) :: ta,tb 
	   IF(pm) THEN; s_t=ta+tb; IF(s_t>beta) s_t=s_t-beta 
	          ELSE;  s_t=ta-tb; IF(s_t<0.) s_t=beta+s_t;  
                         IF(s_t==beta) s_t=0. 
	   END IF 
	 END FUNCTION s_t
	  
!________summing coordinates with periodic B.C._________  
	 SUBROUTINE  shift(k8,k9) 
	 IMPLICIT NONE 
	 INTEGER, INTENT(IN) :: k8,k9 

         IF(triangular==0) then     ! square lattice

         IF(k8.le.dim) then  
      rim(k8)=rim(k8)+1;                
	IF(close_open==0) then            ! trap if close_open=1
	IF(rim(k8)==Ns(k8)) rim(k8)=0     ! homogeneous
      ENDIF
           else; 
      rim(k9)=rim(k9)-1;                
      IF(close_open==0) then            ! trap if close_open=1
	IF(rim(k9)==-1) rim(k9)=Ns(k9)-1  ! homogeneous
      ENDIF
           endif	    

           IF(k8.le.dim) then  
           rim(k8)=rim(k8)+1; IF(rim(k8)==Ns(k8)) rim(k8)=0 
           else; 
           rim(k9)=rim(k9)-1; 
           IF(rim(k9)==-1) rim(k9)=Ns(k9)-1 
           endif	    
          ELSE                      ! triangular lattice
            IF(k8.le.2) then
           rim(k8)=rim(k8)+1; IF(rim(k8)==Ns(k8)) rim(k8)=0
            else if(k8==3) then
           rim(1)=rim(1)-1; IF(rim(1)==-1) rim(1)=Ns(1)-1
           rim(2)=rim(2)+1; IF(rim(2)==Ns(2)) rim(2)=0
            else if(k8<6) then
           rim(k9)=rim(k9)-1;
           IF(rim(k9)==-1) rim(k9)=Ns(k9)-1
            else
           rim(2)=rim(2)-1; IF(rim(2)==-1) rim(2)=Ns(2)-1
           rim(1)=rim(1)+1; IF(rim(1)==Ns(1)) rim(1)=0
            endif
          ENDIF

	  END SUBROUTINE shift

	   
!________ directions of kinks____________________________   
          INTEGER FUNCTION DIRECTION(name,name1) 
          IMPLICIT NONE 
          integer, intent(in) :: name, name1 
          integer :: da 
          
            DIRECTION=0
            do da=1,dir; if (ms(da,name) == name1) then
            DIRECTION=da; exit; endif; enddo
 
           IF(DIRECTION<1 .OR. DIRECTION >dir ) then
            print*, name, name1, DIRECTION; print*
            do da=1,dir; print*, da, ms(da,name); enddo
            stop 'violation of the direction array'
            ENDIF

          END FUNCTION DIRECTION 
             
!________________________________________________________   
            INTEGER FUNCTION oneoutof(numb) 
            IMPLICIT NONE 
            integer, intent(in) :: numb
            oneoutof=rndm()*numb+1.d0
            IF(oneoutof<1 .OR. oneoutof>numb) then
            print*, oneoutof,numb,'procedure =',km
            stop 'oneouto exceeds limits'
            endif
            END FUNCTION  oneoutof        

!_________ names management_________________________________ 
           SUBROUTINE GETNAME(nick)
           IMPLICIT NONE  
	INTEGER, INTENT(OUT) :: nick 
      if (nmnm<kink_max) then 
      nmnm=nmnm+1; nick=nlist(nmnm); numnam(nick)=nmnm 
      else; print*,'GETNAME: The List is over.';stop;endif 
          END SUBROUTINE GETNAME			 
 
 	   SUBROUTINE DROPNAME(nick)
 	         IMPLICIT NONE  
      INTEGER, INTENT(IN) :: nick; INTEGER :: nm,lsn 
	 
      nm=numnam(nick) 
      if(nm<nmnm)then;lsn=nlist(nmnm);nlist(nm)=lsn;numnam(lsn)=nm 
      nlist(nmnm)=nick; numnam(nick)=nmnm; nmnm=nmnm-1 
      else if (nm==nmnm) then; nmnm=nmnm-1 
      else; print*,'DROPNAME:',nick,'not in the List.';stop;endif 
          END SUBROUTINE DROPNAME 

!************************ VITAL PROCEDURES BELOW *********************
!*******************be very careful in modifying them **************** 
         
!_________change associations for the new kink inserted_______________
!_________     (all cases are dealt with separately)   _______________ 	    	     
          SUBROUTINE insert_interval(p0,p1,p2,t1,t2)
          IMPLICIT NONE  
	  INTEGER, INTENT(IN) :: p0,p1,p2 
	  INTEGER :: nnn, nr, j, j1 
	  DOUBLE PRECISION, INTENT(IN) :: t1,t2 
	  DOUBLE PRECISION :: t3,tnnn,t0 
      
      t0=time(p0) 
      IF(t2==tau(p0)) then; 
!________________________________________________________________ 
!   x----------(p0)-----------x   ---> x...(p0)..o--(p1)-------x 
!   			                         t1	       t2 
!________________________________________________________________
 
                         DO j=1,dir; j1=back(j)  !************ 
      nnn=ms(j,p0); tnnn=-s_t(minus,t0,time(nnn)) 
 
      do1: DO; tnnn=tnnn+tau(nnn); IF(tnnn>t1) EXIT do1; nnn=mr(nnn) 
	IF(tnnn==t1) ms(j1,nnn)=p1 
      END DO do1; ms(j,p1)=nnn    
     	   
	do2: DO; IF(tnnn>=t2) EXIT do2; nnn=mr(nnn) 
 	ms(j1,nnn)=p1; tnnn=tnnn+tau(nnn) 
	END DO do2;  
                           END DO                 !***********
                            
! link in time direction 
	nr=mr(p0); mr(p0)=p1; ml(p1)=p0; mr(p1)=nr;  ml(nr)=p1 
	          tau(p0)=t1; tau(p1)=t2-t1 
! absolute times 
                  time(p1)=s_t(plus,t0,t1) 
! site assignment		 
 	          site(p1)=site(p0) 
      else; 
!___________________________________________________________________ 
!   x----------(p0)-----------x   ---> x---(p0)--o..(p1)..o--(p2)--x 
!   t0								 		   t1		t2 
!___________________________________________________________________ 
      tau(p1)=t2-t1      
      IF(p0==p2) THEN	
      !  __________________________
      !  --(p0)--o..(p1)..o--(p0)-- 
      !          t1       t2 
      !____________________________
	 
	 
                        DO j=1,dir; j1=back(j) !*********** 
      nnn=ms(j,p0); tnnn=-s_t(minus,t0,time(nnn)) 
 
	    IF(t1>0.) THEN; time(p1)=s_t(plus,t0,t1) 
      do3: DO;tnnn=tnnn+tau(nnn);IF(tnnn>t1) EXIT do3;nnn=mr(nnn) 
	IF(tnnn==t1) ms(j1,nnn)=p1; ENDDO do3; ms(j,p1)=nnn     
 
		ELSE; time(p1)=s_t(minus,t0,-t1);  ! t1<0 
      do4: DO; IF(tnnn<=t1) EXIT do4;nnn=ml(nnn);tnnn=tnnn-tau(nnn) 
	END DO do4; IF(tnnn==t1) ms(j1,nnn)=p1 
	ms(j,p1)=nnn; tnnn=tnnn+tau(nnn) 
	 
         END IF    ! whether t1 > or < 0 
 
	do5: DO; IF(tnnn>t2) THEN; EXIT do5 
	ELSE IF(tnnn==t2)THEN; nnn=mr(nnn); 
	ms(j1,nnn)=p2;EXIT do5;ENDIF 
	nnn=mr(nnn); ms(j1,nnn)=p1; tnnn=tnnn+tau(nnn) 
	END DO do5;  ms(j,p2)=nnn 
 
                           END DO               !*********** 
 
 ! link in time direction 
                mr(p0)=p1; ml(p1)=p0; mr(p1)=p0; ml(p0)=p1 
	        tau(p0)=beta-tau(p1); time(p2)=s_t(plus,t0,t2) 
! site assignment		 
	        site(p1)=site(p0) 
	                                          
   
	ELSE; 
	!_______________________________ 
	! x-(p0)--o.(p1).o-(p2)-x-(pr)- 
        !         t1     t2 
        !_______________________________
 
      t3=tau(p0);      DO j=1,dir; j1=back(j)   !*********** 
      nnn=ms(j,p0); tnnn=-s_t(minus,t0,time(nnn)) 
 
      do6: DO;tnnn=tnnn+tau(nnn);IF(tnnn>t1) EXIT do6;nnn=mr(nnn) 
           IF(tnnn==t1) ms(j1,nnn)=p1; END DO do6; ms(j,p1)=nnn    
     	   
      do7: DO; IF(tnnn>t2) EXIT do7; nnn=mr(nnn) 
 	   IF(tnnn==t2) THEN; ms(j1,nnn)=p2; 
 	   tnnn=tnnn+tau(nnn);EXIT do7 ;END IF 
 	   ms(j1,nnn)=p1; tnnn=tnnn+tau(nnn) 
	   END DO do7;  ms(j,p2)=nnn 
 
      do8: DO; IF(tnnn>=t3) EXIT do8 
	   nnn=mr(nnn); ms(j1,nnn)=p2; tnnn=tnnn+tau(nnn) 
	   END DO do8 
 
 	                        END DO             !*********** 
! link in time direction 
	nr=mr(p0);  mr(p0)=p1; ml(p1)=p0; mr(p1)=p2; ml(p2)=p1 
	mr(p2)=nr;  ml(nr)=p2; tau(p0)=t1; tau(p2)=t3-t2 
! absolute times 
        time(p1)=s_t(plus,t0,t1); time(p2)=s_t(plus,t0,t2) 
! site assignment		 
 	site(p1)=site(p0); site(p2)=site(p0) 
 
                              END IF 
       endif 
!                             ****** 
	    END SUBROUTINE insert_interval 
	     	 
!______________inverse to the previous procedure____________________	
!______________reconnecting layers when removing an interval________ 
               SUBROUTINE remove_interval(p0,p1,p2,tf) 
               IMPLICIT NONE 
               INTEGER, INTENT(IN) :: p0; 
               INTEGER :: nnn,pr,p1,p2, j, j1 
               DOUBLE PRECISION, INTENT(IN) :: tf; 
               DOUBLE PRECISION :: tnn 
 
!___________________________________________________________________ 
!  --(p0)--x...(p1)..x--(p2)--        -->  ---(p0)---------x--(p2) 
!   		     tf				                      tf 
! 
!  --(p0)--x...(p1)..x--(p2)--x--(pr)  --> ---(p0)---------x--(pr) 
!      					  tf				           tf 
!___________________________________________________________________	 
 
                        DO j=1,dir; j1=back(j)     !********** 
         nnn=ms(j,p1); tnn=-s_t(minus,time(p1),time(nnn)) 
	 IF(tnn==0.) ms(j1,nnn)=p0   
      do1:  DO; tnn=tnn+tau(nnn); IF(tnn>=tf) EXIT do1 
            nnn=mr(nnn); ms(j1,nnn)=p0; END DO do1 
                              END DO               !********** 
 
	IF(p0==p2) THEN; mr(p0)=p0; ml(p0)=p0 
	                tau(p0)=beta;  call DROPNAME(p1) 
	                mk(p0)=0 
	           ELSE; pr=mr(p2); mr(p0)=pr; ml(pr)=p0 
			     tau(p0)=tau(p0)+tf 
	                 call DROPNAME(p1); 
                         IF(p1.ne.p2) call DROPNAME(p2) 
         END IF 
         RETURN 
         END SUBROUTINE remove_interval 	

!________ change associations for the kink miving in time __________ 
          SUBROUTINE move_associations(nl,n0,tk) 
	  IMPLICIT NONE 
	  INTEGER, INTENT(IN) :: n0,nl; 
	  INTEGER :: nnn, j, j1 
	  DOUBLE PRECISION, INTENT(IN) :: tk	 
	  DOUBLE PRECISION :: tnnn,ta 
!___________________________________________________________________ 
!        ----(nl)-----x.....(n0)....  ---> ---(nl)-x.........(n0)... 
!			  0				               tk 
!___________________________________________________________________
        IF(tk==0) RETURN	     
	ta=time(n0); DO j=1,dir; j1=back(j); nnn=ms(j,n0) 
	                tnnn=-s_t(minus,ta,time(nnn)) 
	IF(tk<0) THEN;                    IF(tnnn>tk) THEN 
	do1: DO; ms(j1,nnn)=n0; nnn=ml(nnn); tnnn=tnnn-tau(nnn) 
	IF(tnnn<=tk) EXIT do1; END DO do1;          END IF 
  	IF(tnnn==tk) ms(j1,nnn)=n0 
 
	 ELSE; IF(tnnn==0.) ms(j1,nnn)=nl; tnnn=tnnn+tau(nnn) 
                                          IF(tnnn<tk) THEN 
 	do2: DO;nnn=mr(nnn); ms(j1,nnn)=nl;tnnn=tnnn+tau(nnn) 
	IF(tnnn>=tk) EXIT do2; END DO do2;          END IF 
	IF(tnnn==tk) nnn=mr(nnn)  
 
               END IF;                         ms(j,n0)=nnn 
	                END DO     !************************* 
	tau(nl)=tau(nl)+tk; tau(n0)=tau(n0)-tk;	 
	IF(tk>0) THEN; time(n0)=s_t(plus,ta,tk) 
	         ELSE;  time(n0)=s_t(minus,ta,-tk); END IF 
 
          END SUBROUTINE move_associations 

!************** PROCEDURES BELOW DO AN ALMOST IDENTICAL JOB AS ABOVE
!************** BY PREPARING THE POTENTIAL ENERGY CHANGE FOR THE UPDATES 
!************** (THE CONFIGURATION IS NOT UPDATED) 

!________potential change  for the new kink inserted_________________ 
         SUBROUTINE potential_insert(p0,t1,t2) 
	 IMPLICIT NONE 
	 INTEGER, INTENT(IN) :: p0 
	 INTEGER :: nnn,her , j, j1 
	 DOUBLE PRECISION, INTENT(IN) :: t1,t2 
	 DOUBLE PRECISION :: tnnn,dt,Iz(dir) 
	 
         pot=0.d0; her=site(p0); 
         DO j=1,dir; Iz(j)= Int_z(her,j); ENDDO 

!________________________________________________________________ 
!	x---------(p0)-----------x   ---> x---(p0)--o..(p1)..o 
!       t0	         	                        t1       t2 
!________________________________________________________________
 
        DO j=1,dir; j1=back(j); nnn=ms(j,p0) !*************** 
	tnnn=-s_t(minus,time(p0),time(nnn)) 
 
	IF(t1>0.) THEN; 
	do1: DO; tnnn=tnnn+tau(nnn); IF(tnnn>t1) EXIT do1 
	                nnn=mr(nnn);	END DO do1 
	ELSE; 
	do2: DO; IF(tnnn<=t1) EXIT do2; 
	nnn=ml(nnn); tnnn=tnnn-tau(nnn); 
	END DO do2; tnnn=tnnn+tau(nnn); END IF 
	pot=pot+Iz(j)*Fss(fl(nnn))*(tnnn-t1) 
 
	IF(tnnn<t2) THEN; 
	do3: DO; nnn=mr(nnn); dt=tau(nnn); tnnn=tnnn+dt 
	pot=pot+Iz(j)*Fss(fl(nnn))*dt; IF(tnnn>=t2) EXIT do3; 
	END DO do3; END IF 
	pot=pot+Iz(j)*Fss(fl(nnn))*(t2-tnnn); 
	END DO;                               !************** 
 
	    END SUBROUTINE potential_insert 
	    
!________________________________________________________________ 
!_______ potential change when removing an interval______________ 
         SUBROUTINE potential_remove(p1,tf) 
         IMPLICIT NONE 
         INTEGER :: nnn,p1,her, j, j1 
         DOUBLE PRECISION, INTENT(IN) :: tf; 
	 DOUBLE PRECISION :: tnnn, dt,Iz(dir) 
!________________________________________________________________ 
!        --(p0)--x...(p1)..x--(p2)--   -->  ---(p0)----x--(p2->p0) 
! 			       tf  			               tf 
!________________________________________________________________

      pot=0.d0; her=site(p1); 
      DO j=1,dir; Iz(j)= Int_z(her,j); ENDDO 
      
        DO j=1,dir; j1=back(j); nnn=ms(j,p1) !******** 
	tnnn=-s_t(minus,time(p1),time(nnn))+tau(nnn) 
        pot=pot+Iz(j)*Fss(fl(nnn))*tnnn 
                  IF(tnnn<tf) THEN; 
	do1: DO; nnn=mr(nnn); dt=tau(nnn); tnnn=tnnn+dt 
	pot=pot+Iz(j)*Fss(fl(nnn))*dt; IF(tnnn>=tf) EXIT do1; 
	END DO do1;          END IF 
 
      pot=pot+Iz(j)*Fss(fl(nnn))*(tf-tnnn)     
      END DO                                  !******* 
      RETURN     
                  END SUBROUTINE potential_remove 
      
!________________________________________________________________ 
!________ potential change for the kink miving in time___________ 
          SUBROUTINE potential_move(nl,n0,tk)
          IMPLICIT NONE  
	  INTEGER, INTENT(IN) :: n0,nl; 
	  INTEGER :: nnn,her, j, j1 
	  DOUBLE PRECISION, INTENT(IN) :: tk 
	  DOUBLE PRECISION :: tnnn,dt,Iz(dir) 
!________________________________________________________________ 
!        ----(nl)-----x.....(n0)....  ---> ---(nl)-x.........(n0)... 
!			  0				               tk 
!________________________________________________________________
 
      pot=0.d0; her=site(nl); 
      DO j=1,dir; Iz(j)= Int_z(her,j); ENDDO 
      
      DO j=1,dir; j1=back(j); nnn=ms(j,n0)          !************** 
	 tnnn=-s_t(minus,time(n0),time(nnn)) 
 
	    IF(tk<0) THEN; pot=pot+Iz(j)*Fss(fl(nnn))*tnnn 
	            IF(tnnn>tk) THEN; 
	do1: DO; nnn=ml(nnn); dt=tau(nnn); tnnn=tnnn-dt 
	pot=pot-Iz(j)*Fss(fl(nnn))*dt; IF(tnnn<=tk) EXIT do1; 
	END DO do1;            END IF 
  	 
            ELSE; tnnn=tnnn+tau(nnn); 
         pot=pot+Iz(j)*Fss(fl(nnn))*tnnn 
                    IF(tnnn<tk) THEN; 
         do2: DO; nnn=mr(nnn); dt=tau(nnn); tnnn=tnnn+dt 
	 pot=pot+Iz(j)*Fss(fl(nnn))*dt; IF(tnnn>=tk) EXIT do2; 
	 END DO do2;          END IF 
            END IF;   
       pot=pot+Iz(j)*Fss(fl(nnn))*(tk-tnnn); END DO  !************** 
	 
         END SUBROUTINE potential_move
          
          
!********* UPDATES USING ALL OF THE ABOVE SUBROUTINES *******          
          
!_________                            _______________________
            SUBROUTINE create_worm(df) 
            IMPLICIT NONE 
            INTEGER, INTENT(IN) :: df  
            INTEGER :: n0, ngive, n1, n2, fl0,fl1,here 
	    DOUBLE PRECISION :: t0,t1,t2,tau3,tau0,pp,tv,tck 
 
!	  select interval on which to create a WORM   
      ngive=oneoutof(nmnm); 
      n0=nlist(ngive); fl0=fl(n0);  IF(fl0<0) RETURN                  
                                        !**** worm is not created on a hole 
      fl1=fl0+df 
      IF(fl1<0 .or. fl1>fl_max) RETURN  !**** obey the Hilbert space allowed 
 
!  select times for the WORM
!_________________________________________
!    t0        t1        t2	t0+tau0 
!     x---------o.........o-------------x 
!        tau1      tau2      tau3 
!_________________________________________
        tau0=tau(n0); t1=tau0*rndm(); tau3=rndm()*tau0 
	t1=AINT(t1);  tau3=AINT(tau3) 
	IF( t1==tau3 .or. t1==0. .or. tau3==tau0 ) RETURN 
	IF( t1==tau0 .or. tau3==0.) RETURN         !*** to avoid zero-time intervals 
	  
	  IF(tau3<t1) THEN; 
	  IF(ml(n0)==n0) then; t1=t1-beta; t2=tau3 
	                 else; t2=t1; t1=tau3;      ENDIF 
	              ELSE;    t2=tau3 
	  END IF 
 	  t0=time(n0); tv=t2-t1                     !*** interval parameters prepared
 	
!________________________________________________________________________________
        tck=prtcls+(Pat(fl1)-Pat(fl0))*tv   
        tck=tck/beta
        IF(CANONICAL .AND. ABS(tck-Nsite) > Mzmark) RETURN 
!_____________________________________________________________ CANONICAL !!!!!!!   

!DECISION 	 
        p_a=nmnm*wmso*w_hop(fl0,fl1)**2; 
        IF(n0==ml(n0)) p_a=p_a*2.d0 
        p_a=p_a*tau0*tau0;         
        here=site(n0);    
        
        IF(onsite_nn==1) then; pp=0.d0; else
        call potential_insert(n0,t1,t2) 
	pp=pot*(Fss(fl1)-Fss(fl0)) 
        ENDIF                    
    
        pp=pp+U0*(Fs2(fl1)-Fs2(fl0))*tv
        pp=pp-Hsite(here)*(Fss(fl1)-Fss(fl0))*tv 
        p_a=p_a*EFUN(-pp)  

        aa(km)=aa(km)+1	                           !*** counter of attempts 
    	IF(p_a<1.d0 .and. rndm()>p_a) RETURN 
    	bb(km)=bb(km)+1                            !*** counter of acc. attempts
    	 
!UPDATE	          
	call GETNAME(n1) 
        IF(tau0==beta) THEN;         
!_________________________________        
!          -------o.......o-----
!             n0     n1      n0 
!_________________________________
! link in space-time 
        call insert_interval(n0,n1,n0,t1,t2) 
	mk(n0)=0; mk(n1)=0; fl(n1)=fl1    
! flavor 
      SELECT CASE(df) 
	CASE(1);  masha=n1;  ira=n0; present=.TRUE.; tim=tv 
	CASE(-1); masha=n0; ira=n1;  present=.TRUE.; tim=beta-tv 
      END SELECT   
	ELSE; call GETNAME(n2);
!___________________________________	            
!     x-----o.....o-----x 
!       n0     n1    n2
!___________________________________  
! link in space-time 
        call insert_interval(n0,n1,n2,t1,t2) 
	mk(n1)=0; mk(n2)=0; fl(n1)=fl1; fl(n2)=fl0  
! flavor 
      SELECT CASE(df) 
	CASE(1);  masha=n1; ira=n2; present=.TRUE.; tim=tv 
	CASE(-1); masha=n2; ira=n1; present=.TRUE.; tim=beta-tv 
      END SELECT   
	END IF; 
	prtcls=prtcls+(Pat(fl1)-Pat(fl0))*tau(n1) 
	Udiag=Udiag+pp; 

!___________________________________________________________ 
!               masha    ira 
!         ....o_______o......    appear in the configuration
!___________________________________________________________

             END SUBROUTINE create_worm 
              
!__________                            ________________________
           SUBROUTINE delete_worm(imya)         
           IMPLICIT NONE 
           INTEGER, INTENT(IN) :: imya; 
           INTEGER :: n0, n1, n2, n5, fl0,fl1 
           DOUBLE PRECISION :: tau0,t1,t2,t3,pp,tck 
 
! analysing the case and leaving if the case is wrong       
      IF(imya==ira) THEN;         n1=ira; n2=masha; 
      ELSE  IF(imya==masha) THEN; n1=masha; n2=ira; 
      END IF      
      IF(n2.ne.mr(n1)) RETURN       
      
      t1=0; t2=tau(n1); n0=ml(n1); fl0=fl(n0); fl1=fl(n1) 
      here=site(n1)
!________________________________________________________________________________
        tck=prtcls+(Pat(fl0)-Pat(fl1))*t2     
        tck=tck/beta
        IF(CANONICAL .AND. ABS(tck-Nsite) > Mzmark) RETURN 
!_____________________________________________________________ CANONICAL !!!!!!!    


       IF(onsite_nn==1) then; pp=0.d0; else
	call potential_remove(n1,t2) 
	pp=pot*(Fss(fl1)-Fss(fl0)) 
       ENDIF
	pp=pp+U0*(Fs2(fl1)-Fs2(fl0))*t2
	pp=pp-Hsite(here)*(Fss(fl1)-Fss(fl0))*t2	
        p_a=EFUN(-pp)*wmso*w_hop(fl0,fl1)**2 


!DECISION          
        IF(n2==n0) THEN; 
!________________________________
!         n0    n1    n0 
!       -----o......o-----
!________________________________ 
			   
        tau0=beta; p_a=p_a*2.d0; p_a=p_a*(nmnm-1)*tau0*tau0
        aa(km)=aa(km)+1                     !*** number of attempts  
	p_a=1.d0/p_a; IF(p_a<1. .and. rndm()>p_a) RETURN 
        bb(km)=bb(km)+1                     !*** number of acc. attempts 
	call remove_interval(n0,n1,n2,t2); 
 
        ELSE; 
!________________________________
!           n0    n1    n2 
!        x-----o......o-----x 
!________________________________
        t3=t2+tau(n2); tau0=t3+tau(n0); p_a=p_a*(nmnm-2)*tau0*tau0 
        aa(km)=aa(km)+1                     !*** number of attempts 
	p_a=1.d0/p_a; IF(p_a<1. .and. rndm()>p_a) RETURN 
        bb(km)=bb(km)+1                     !*** number of accepted attempts 
        call remove_interval(n0,n1,n2,t3) 
        END IF 

!UPDATE 
        present=.FALSE.; masha=0; ira=0; tim=0.d0 
	prtcls=prtcls+(Pat(fl0)-Pat(fl1))*tau(n1) 
	Udiag=Udiag-pp 
 
        END SUBROUTINE delete_worm   !**************** ira-masha are gone        
           
!_________jump for ira or reconnection for masha __________________
!_________(insert kink to the left of the worm end)________________
          SUBROUTINE insert_left(imya)          
	  IMPLICIT NONE 
	  INTEGER, INTENT(INOUT) :: imya 
	  INTEGER :: n0,nn,nl,n1,n2,nnr,fl0,fll,flnn,fln1,kbb,k  
	  DOUBLE PRECISION :: tl,tau0,tnn,tk,t1,t2,pp,ppe 
!__________________________________________________________________________ 
!  tnn								  tk 
!..x...(nn)............x(nnr)     ..x.(nn).x-(n1)-o.(n2).x(nnr)	 
!			                -->          |		          |     ^  
! x--(nl)------o..(n0)..x(nr)	     x--(nl)-x....(n0)......x(nr)	  V k1  | k
! tl          t=0 
! ______________________________________new kink___________________________
   
        n0=imya;  
        nl=ml(n0);tl=-tau(nl);
        k=oneoutof(dir); nn=ms(k,n0); flnn=fl(nn); 
        IF(flnn==-3) RETURN;                 !*** no worm on hole
        fl0=fl(n0); fll=fl(nl); dfl=fll-fl0; fln1=flnn+dfl 
         
        IF(fln1<0.or.fln1>fl_max) RETURN     !*** obey Hilbert space            
        nnr=mr(nn); tnn=-s_t(minus,time(n0),time(nn))
 	IF(nnr==nn) THEN; tau0=tl; ELSE; tau0=max(tl,tnn); END IF         
 
!give the time position for the new kink 
       tk=rndm()*tau0; tk=AINT(tk); IF(tk==0 .or. tk==tau0) RETURN 
       t1=tk-tnn; t2=-tnn 
       
!DECISION 
        here=site(n0); SELECT CASE(dfl); 
        CASE(1);  thop=s_h(fll,flnn)*Int_xy(here,k); 
        CASE(-1); thop=s_h(flnn,fll)*Int_xy(here,k); 
                       END SELECT       
        p_a=-thop*dir*tau0*w_hop(fln1,flnn)/w_hop(fll,fl0);
        IF(onsite_nn==1) then; pp=0.d0; else  
         call potential_insert(nn,t1,t2); 
         pp=pot*(Fss(fln1)-Fss(flnn)) 
   	 call potential_move(nl,n0,tk);  
   	 pp=pp+pot*(Fss(fll)-Fss(fl0)) 
   	 ppe=Int_z(here,k)*tk*(Fss(fl0)-Fss(fll))*(Fss(fln1)-Fss(flnn)) 
   	 pp=pp-ppe       ! correct for double counting of the ppe term
        ENDIF 

        pp=pp+U0*ABS(tk)*(Fs2(fln1)-Fs2(flnn)+Fs2(fl0)-Fs2(fll)) 
        pp=pp-Hsite(here)*(Fss(fl0)-Fss(fll))*ABS(tk) 
        pp=pp-Hsite(site(nn))*(Fss(fln1)-Fss(flnn))*ABS(tk) 
        
        aa(km)=aa(km)+1  
	p_a=p_a*EFUN(-pp); IF(p_a<1. .and. rndm()>p_a) RETURN 
        bb(km)=bb(km)+1  

!UPDATE 
        call move_associations(nl,n0,tk) 
	call GETNAME(n1) 
	IF(nn==nnr) THEN; n2=nn; ELSE; call GETNAME(n2); END IF 
	call insert_interval(nn,n1,n2,t1,t2) 
      
	fl(n1)=fln1; fl(n2)=flnn; 
	mk(n0)=n1; mk(n1)=n0; mk(n2)=0; imya=n2; 
       
	kbb=back(k) 
	SELECT CASE(dfl); 
	CASE(1);  arr(k)=arr(k)+1;     call shift(k,kbb) 
        CASE(-1); arr(kbb)=arr(kbb)+1; call shift(kbb,k) 
        END SELECT 	 
        kinks=kinks+1; Udiag=Udiag+pp 
 
        IF(ms(k,n0)/=n1  .OR. ms(kbb,n1)/=n0) then
        print*, 'insert left links violated!'
        stop
        ENDIF

	END SUBROUTINE insert_left 
	      
!_________anti-jump of ira or anti-reconnection for masha_________ 
!_________(delete kink to the left of the worm end)_______________ 
          SUBROUTINE delete_left(imya)   
	  IMPLICIT NONE 
	  INTEGER, INTENT(INOUT) :: imya 
	  INTEGER :: n0,nn,nl,n1,n2,fl0,fll,flnn,fln1,fln2,kbb,k 
	  DOUBLE PRECISION :: tau0,t1,tk,t2,pp,ppe 
!_____________________________________________________________________ 
!  t            0       t1                                imya 
!..x....(nn)............x-(nnr)      ..x.(nn).x--(n1)--o..(n2).x-(nnr) 
!	          	        <--           |			 
! x--(nl)-------o..(n0)..x-(nr)	      x-(nl)--x.......(n0)......x-(nr) 
!               t1 
!_____________________________________________________________________
 
        n2=imya; n1=ml(n2); IF(mk(n1)==0) RETURN;   !*** nothing to delete
        nn=ml(n1); fln2=fl(n2); flnn=fl(nn) 
        IF(flnn/=fln2) RETURN 	
        fln1=fl(n1); t1=tau(n1); n0=mk(n1) 
        IF(tau(n0)<=t1) RETURN;                     !*** wrong case 
        fl0=fl(n0); dfl=fln1-fln2 
        nl=ml(n0); fll=fl(nl); tau0=t1+min(tau(nn),tau(nl)) 
        k=direction(n0,n1)

!DECISION 
        here=site(n0);  SELECT CASE(dfl);  
        CASE(1);  thop=s_h(fll,flnn)*Int_xy(here,k); 
        CASE(-1); thop=s_h(flnn,fll)*Int_xy(here,k) 
                         END SELECT 
        p_a=thop*dir*tau0*w_hop(fln1,fln2)/w_hop(fll,fl0)
        IF(onsite_nn==1) then; pp=0.d0; else 
        call potential_remove(n1,t1);  
        pp=pot*(Fss(fln1)-Fss(flnn)) 
        call potential_move(nl,n0,t1);    
        pp=pp-pot*(Fss(fll)-Fss(fl0)) 
        ppe=Int_z(here,k)*t1*(Fss(fll)-Fss(fl0))*(Fss(flnn)-Fss(fln1)) 
        pp=pp-ppe 
        ENDIF
 
        pp=pp+U0*t1*(Fs2(fln1)-Fs2(flnn)+Fs2(fl0)-Fs2(fll)) 
        pp=pp-Hsite(here)*(Fss(fl0)-Fss(fll))*t1 
        pp=pp-Hsite(site(nn))*(Fss(fln1)-Fss(flnn))*t1 
         
	p_a=p_a*EFUN(-pp) 
        aa(km)=aa(km)+1  
	p_a=1./p_a;	IF(p_a<1..and.rndm()>p_a) RETURN 
        bb(km)=bb(km)+1  
        
!UPDATE  
        IF(nn==n2) THEN;     call remove_interval(nn,n1,n2,t1) 
	ELSE; t2=t1+tau(n2); call remove_interval(nn,n1,n2,t2) 
	END IF 
	                     call move_associations(nl,n0,t1); 
	mk(n0)=0; imya=n0  
	kbb=back(k); SELECT CASE(dfl); 
	CASE(1);  arr(k)=arr(k)-1;     call shift(kbb,k) 
        CASE(-1); arr(kbb)=arr(kbb)-1; call shift(k,kbb) 
                     END SELECT 
         		 
  	kinks=kinks-1;  Udiag=Udiag-pp 
 
	    END SUBROUTINE delete_left 
 
!________ jump for ira or reconnection for masha____________________
!________ (insert kink to the left of the worm)_____________________
          SUBROUTINE insert_right(imya) 
	  IMPLICIT NONE 
	  INTEGER, INTENT(INOUT) :: imya 
	  INTEGER :: n0,nn,nl,n1,n2,nnr,fl0,fll,flnn,fln1,kbb,k 
	  DOUBLE PRECISION :: tr,tau0,tnn,tk,t1,t2,pp,ppe 
!____________________________________________________________________________ 
!   tnn								      0      tk 
!...x....(nn).............x(nnr)    ..x.(nn)..o-(n1)-x..(n2)..x(nnr)	 
!			     	--> 	             |		         ^ k  | 
! x..(nl).......o--(n0)--x(nr)       x...(nl)........x--(n0)-x(nr)	 |    v k1 
! tl            0        tr
!_____________________________________________________________________________
         
        n0=imya; nl=ml(n0); tr=tau(n0);
        k=oneoutof(dir); nn=ms(k,n0); flnn=fl(nn); IF(flnn==-3) RETURN
        fl0=fl(n0); fll=fl(nl); dfl=fll-fl0; fln1=flnn-dfl  
	IF(fln1<0.or.fln1>fl_max) RETURN 
         
 	nnr=mr(nn); tnn=-s_t(minus,time(n0),time(nn)) 
 	IF(nnr==nn) THEN; tau0=tr; 
 	ELSE; tau0=min(tr,tnn+tau(nn));END IF  	
 
! give the time position for the new kink 
      tk=rndm()*tau0; tk=AINT(tk); IF(tk==0 .or. tk==tau0) RETURN 
      t1=-tnn; t2=s_t(plus,t1,tk); IF(t1>t2) t1=t1-beta 
	 
!DECISION 
        here=site(n0);  SELECT CASE(dfl); 
        CASE(1); thop=s_h(fll,fln1)*Int_xy(here,k) ; 
        CASE(-1);thop=s_h(fln1,fll)*Int_xy(here,k) ; 
                        END SELECT
                         
        p_a=thop*dir*tau0*w_hop(fln1,flnn)/w_hop(fll,fl0)
        IF(onsite_nn==1) then; pp=0.d0; else 
         call potential_insert(nn,t1,t2); 
         pp=pot*(Fss(fln1)-Fss(flnn)) 
         call potential_move(nl,n0,tk);   pp=pp+pot*(Fss(fll)-Fss(fl0)) 
         ppe=Int_z(here,k)*tk*(Fss(fll)-Fss(fl0))*(Fss(fln1)-Fss(flnn)) 
	  pp=pp+ppe 
        ENDIF 

        pp=pp+U0*tk*(Fs2(fll)-Fs2(fl0)+Fs2(fln1)-Fs2(flnn)) 
        pp=pp-Hsite(here)*(Fss(fll)-Fss(fl0))*tk 
        pp=pp-Hsite(site(nn))*(Fss(fln1)-Fss(flnn))*tk 
          
        aa(km)=aa(km)+1   
	p_a=p_a*EFUN(-pp); IF(p_a<1. .and. rndm()>p_a) RETURN 
        bb(km)=bb(km)+1  

!UPDATE    
        call move_associations(nl,n0,tk) 
	call GETNAME(n1); 
	IF(nn==nnr) THEN; n2=nn; ELSE; call GETNAME(n2); END IF 
	call insert_interval(nn,n1,n2,t1,t2) 
      
	fl(n1)=fln1; fl(n2)=flnn; 
	mk(n0)=n2; mk(n2)=n0; mk(n1)=0; imya=n1; 
	kbb=back(k); SELECT CASE(dfl); 
	CASE(1);  arr(k)=arr(k)+1;     call shift(k,kbb) 
        CASE(-1); arr(kbb)=arr(kbb)+1; call shift(kbb,k) 
                    END SELECT 	 
  	Udiag=Udiag+pp; kinks=kinks+1 

        IF(ms(k,n0)/=n2  .OR. ms(kbb,n2)/=n0) then
        print*, 'insert right links violated!'
        stop
        ENDIF         
 
 	    END SUBROUTINE insert_right 
  
!_________anti-jump for masha or anti-reconnection for ira_____________
!_________(delete kink to the right of the worm end)___________________   
          SUBROUTINE delete_right(imya) 
	  IMPLICIT NONE 
	  INTEGER, INTENT(INOUT) :: imya 
 	  INTEGER :: n0,nn,nl,n1,n2,fl0,fll,flnn,fln1,fln2,kbb,k 
	  DOUBLE PRECISION :: tau0,t1,t2,pp,ppe 
!_____________________________________________________________________ 
!  tnn				                         0 imya t1 
!...x....(nn)...........x(nnr)	       ..x.(nn)..o-(n1).x.(n2)..x(nnr) 
!			        <--- 	                |	 
! x..(nl).......o--(n0)--x(nr)         x...(nl).........x--(n0)--x(nr) 
!               0 
!_____________________________________________________________________
	n1=imya; nn=ml(n1); n2=mr(n1); IF(mk(n2)==0) RETURN 
	fln2=fl(n2); flnn=fl(nn); fln1=fl(n1); 
	IF(flnn/=fln2) RETURN
	fln1=fl(n1); t1=tau(n1); n0=mk(n2); nl=ml(n0) 
	IF(t1>=tau(nl)) RETURN                          !*** wrong case 
	fl0=fl(n0); fll=fl(nl);  dfl=flnn-fln1 	 
	tau0=t1+min(tau(n2),tau(n0)); k=direction(n0,n2) 
	 
!DECISION 
        here=site(n0); SELECT CASE(dfl); 
        CASE(1); thop=s_h(fll,fln1)*Int_xy(here,k) ; 
        CASE(-1);thop=s_h(fln1,fll)*Int_xy(here,k) ; 
                       END SELECT 
        p_a=thop*dir*tau0*w_hop(fln1,fln2)/w_hop(fll,fl0)
        IF(onsite_nn==1) then; pp=0.d0; else 
        call potential_remove(n1,t1); 
        pp=pot*(Fss(fln1)-Fss(flnn)) 
   	 call potential_move(nl,n0,-t1);  pp=pp+pot*(Fss(fl0)-Fss(fll)) 
   	 ppe=Int_z(here,k)*t1*(Fss(fl0)-Fss(fll))*(Fss(flnn)-Fss(fln1)) 
   	 pp=pp-ppe 
        ENDIF

        pp=pp+U0*t1*(Fs2(fll)-Fs2(fl0)+Fs2(fln1)-Fs2(flnn)) 
        pp=pp-Hsite(here)*(Fss(fll)-Fss(fl0))*t1 
        pp=pp-Hsite(site(nn))*(Fss(fln1)-Fss(flnn))*t1 
           
	p_a=p_a*EFUN(-pp) 
        aa(km)=aa(km)+1  
     	p_a=1./p_a;	IF(p_a<1..and.rndm()>p_a) RETURN 
        bb(km)=bb(km)+1  

!UPDATE   
        IF(nn==n2) THEN;     call remove_interval(nn,n1,n2,t1) 
	ELSE; t2=t1+tau(n2); call remove_interval(nn,n1,n2,t2) 
	END IF 
	                     call move_associations(nl,n0,-t1)
	mk(n0)=0; imya=n0 	 
	kbb=back(k);  SELECT CASE(dfl); 
	CASE(1);  arr(k)=arr(k)-1;       call shift(kbb,k) 
        CASE(-1); arr(kbb)=arr(kbb)-1;   call shift(k,kbb) 
                      END SELECT 		 
    	Udiag=Udiag-pp;  kinks=kinks-1 			 
    	    			 
	    END SUBROUTINE delete_right 
 
 
!______   move the worm end in time __________________________________ 
          SUBROUTINE move_worm(n0)
          IMPLICIT NONE 
	  INTEGER, INTENT(IN) :: n0; INTEGER :: nl 
	  DOUBLE PRECISION :: tk,t1,t2,pp,tv,tck 
!____________________________________________________________ 
!----(nl)-----x.....(n0)....  --->  ---(nl)-x.........(n0)... 
!	      0					            tk 
!____________________________________________________________	 
     	nl=ml(n0); t1=tau(nl); t2=t1+tau(n0); 
     	tk=rndm()*t2; tk=AINT(tk) 
	IF(tk==0 .or. tk==t2 .or. tk==t1) RETURN; tk=-t1+tk

!_________________________________________________________
       tck=prtcls+tk*(Pat(fl(nl))-Pat(fl(n0)))   
        tck=tck/beta
        IF(CANONICAL .AND. ABS(tck-Nsite) > Mzmark) RETURN 
!________________________________________________________  ! CANONICAL !!!!!!!!!!
        

!DECISION 
      IF(onsite_nn==1) then; pot=0.d0; else
      call potential_move(nl,n0,tk);
      ENDIF 

      here=site(n0) 
      pp=(-Hsite(here)*tk+pot)*(Fss(fl(nl))-Fss(fl(n0))) 
      pp=pp+U0*(Fs2(fl(nl))-Fs2(fl(n0)))*tk
      p_a=EFUN(-pp) 
        aa(km)=aa(km)+1  
      IF(p_a<1. .and. rndm()>p_a) RETURN 
        bb(km)=bb(km)+1  

!UPDATE  
       call move_associations(nl,n0,tk) 
       prtcls=prtcls+tk*(Pat(fl(nl))-Pat(fl(n0))) 
       Udiag=Udiag+pp; 
       IF(n0==ira) then; tim=tim+tk;        
       else; tim=tim-tk; endif
       IF(tim < 0)  tim=tim+beta  
       IF(tim > beta) tim=tim-beta  		   
	 
		    END SUBROUTINE move_worm 

!*************************TESTING PROGRAM ***************************
! links, correspondence, beta-loops, times, etc. ********************
	SUBROUTINE testing
        IMPLICIT NONE
	INTEGER :: name, next, na, nnn,flne, j2,her,i,k,j1,kw
	DOUBLE PRECISION :: t0,tl,ts,tnnn,tf,dt,add,ad0
	DOUBLE PRECISION :: diag,prt

         print*, 'in testing; km=', km
!       go to 234           !***** if you want to skip some testing

!testing energies
        diag=0.d0; ad0=0.d0
                      DO j2=1,nmnm; next=nlist(j2); tf=tau(next)
                                      flne=fl(next); her=site(next)
        add=0.d0;   
        IF(onsite_nn==2) then
        DO j=1,dir  
	nnn=ms(j,next); tnnn=-s_t(minus,time(next),time(nnn))+tau(nnn)
	add=add+Int_z(her,j)*Fss(fl(nnn))*tnnn
	IF(tnnn<tf) THEN; a: DO; nnn=mr(nnn); dt=tau(nnn); tnnn=tnnn+dt 
	add=add+Int_z(her,j)*Fss(fl(nnn))*dt; 
	IF(tnnn>=tf) EXIT a; END DO a; END IF
        add=add+Int_z(her,j)*Fss(fl(nnn))*(tf-tnnn)
                      END DO;
        ENDIF              
        prt=add*Fss(flne)
        !IF(ABS(prt) > 10.d0 ) then
        !print*, prt
        !print*, Fss(flne)
        !print*, add
        !print*, tf/t_quanta
        !print*, '___________________'
        !pause
        !ENDIF
        diag=diag+add*Fss(flne);  
        ad0=ad0+U0*tf*Fs2(flne)-Hsite(her)*tf*Fss(flne)
                            END DO
 
	diag=diag/2.d0 ;
	IF(ABS(Udiag-diag-ad0)>1.d-4) then
	print*, 'Udiag-diag-ad0 .ne. 0'
        print*,  Udiag-diag-ad0
!        stop
        Udiag=diag+ad0;                !****** pause 
                                       !**** uncomment if serious 
        endif
 
!testing particle number 
        prt=0.	
	DO j=1,nmnm
	name=nlist(j); next=mk(name); prt=prt+Pat(fl(name))*tau(name)
	IF(next/=0) THEN;  do i=1,dir+1
	IF(i==dir+1) THEN
	print*, 'kink connections failed',name,next; stop; END IF
	IF(next==ms(i,name)) EXIT; end do;  END IF; END DO 
	IF(ABS(prt-prtcls)>1.d-1) then
	print*, 'prt-prtcls .ne. 0', km;
!        stop 
	prtcls=prt;                     !****** pause    !**** uncomment if serious  
          endif
         
!        go to 234         !***** if you want to skip some testing

!links				 
	DO j=1,nmnm; name=nlist(j); t0=time(name)
	do k=1,dir; na=ms(k,name); tl=-s_t(minus,t0,time(na))
	ts=tl+tau(na)
	   IF(ts<=0.) THEN
	   print*, 'names are not in order - times', name, k, na
	   print*, name,t0
	   print*, na,time(na),tau(na)
	   print*, 'nmnm=', nmnm
	   print*, 'procedure ',km
     	   stop '****** serious violation of time ordering******' 
     	   END IF

	IF(numnam(name)>nmnm .or. numnam(na)>nmnm) THEN
	print*, 'names are not in order ????', ira,masha,name,na
	print*, nmnm, numnam(name), numnam(na), k; 
	stop '****** serious violation of name ordering******' 
	ENDIF
	end do; END DO

	DO j=1,nmnm
        name=nlist(j);         t0=tau(name);
        IF(t0<=0) then; print*, 'tau()<=0 ', km; 
        stop '**serious violation of positive interval time length**' 
        ENDIF 
        next=name
	DO
	next=mr(next); IF(next==name) EXIT; t0=t0+tau(next); END DO
	IF(t0/=beta) then; print*, 'beta violation', t0,beta,name,km;
        print*, name, tau(name), mr(name), tau(mr(name)) 
	stop '**serious violation of the beta loop sum rule **' 
	ENDIF
	
!layers
        DO j1=1,dim
	next=name; DO k=1,Ns(j1); next=ms(j1,next); END DO
	IF(site(next)/=site(name)) THEN; print*, 'layers violation'; 
	stop '**serious violation of periodic BC**'
	END IF
	END DO

        END DO
   
!worms    
   
       IF(.not.present) then 
       IF(ira.ne.0 .or. masha.ne.0 ) then
       print*, 'ira-masha .ne 0'
       print*, ira,masha,km; 
       stop '**serious violation: do we have worms or not?**' 
       endif        
       DO j=1,dim; k=back(j); 
       IF(triangular==0) then    ! square lattice
          kw=(arr(j)-arr(k))
       ELSE                      ! triangular lattice
          if(j==1) then
          kw=(arr(j)-arr(k)+arr(6)-arr(3))
          else
          kw=(arr(j)-arr(k)+arr(3)-arr(6))
          endif
       ENDIF
       i=kw/Ns(j); i=i*Ns(j);
       IF(i-kw.ne.0) then
       print*, 'winding numbers not integer' 
       print*, j,k,Ns(j),kw, arr(j),arr(k);
       print*, 'procedure ',km,kmm 
       print*, ' ira-masha are now at '
       print*,  Rim, tim/beta 
       stop  '**serious violation of integer winding numbers**'
       endif; 
       ENDDO 
       else; 
        IF(ira==0 .OR. masha==0) then;
        print*, 'present, but ira or masha =0'
        print*, present,ira,masha,km;
        stop '**serious violation: ira-masha names not in order'
        endif
       endif 
        
! 234    continue 
       print*, 'out of testing'
	    END SUBROUTINE testing
	    
!______________________________________________________________________	    

!***********VERY IMPORTANT !!! - HAMILTONIAN IS CONSTRUCTED HERE********
!___________                      ______________________________________ 
            SUBROUTINE hamiltonian
            IMPLICIT NONE
            integer :: ir,ix,iy,iz,kf,kb,kcf,kcb,irr 
	    integer :: i,j1,n12,n1,nn(dim)
	    double precision :: rwall, rx, ry, radius
            
! H = sum_{ij} U(i-j) Fss(n(i)) Fss(n(j))
        Fss=0.d0; Fs2=0.d0
        DO j=0,fl_max; Fss(j)=j*1.d0; Fs2(j)=Fss(j)**2; ENDDO
! here Fss(n) = n --> density-density coupling 
           Fss(-3)=fl_max/2.d0;  Fs2(-3)= Fss(-3)**2 
!for the hole

          IF(spin_or_boson==1) then 
            s_h=0.d0; DO j=1,fl_max; DO j1=0,fl_max-1
            s_h(j,j1)=sqrt(j*(j1+1.d0)) 
            ENDDO; ENDDO       
            !***********tab. matr. el. for bosons b1^+ b |n=j,n1=j1> 
          ELSE 
            s_h=0.d0; DO j=1,fl_max; DO j1=0,fl_max-1
            s_h(j,j1)=sqrt(j*(fl_max-j+1.d0))
            s_h(j,j1)=s_h(j,j1)*sqrt((j1+1.d0)*(fl_max-j1)) 
            ENDDO; ENDDO         
            !********tab. matr. el.  for spins S1^+ S^-|S_z=j,S1_z=j1> 
          ENDIF
             
         wmso=source*w_del*2.d0; w_hop=0.d0

         IF(spin_or_boson==1) then 
           do j=0,fl_max; DO j1=0,fl_max         
           IF(j>j1) then; w_hop(j,j1)=sqrt(j*1.d0); 
           else IF(j1>j) then; w_hop(j,j1)=sqrt(j1*1.d0); endif       
           enddo; enddo 
           !********tab. matr. el. for the worm (bose)
           
           do j=0,fl_max; Pat(j)=j*1.d0; enddo;  Pat(-3)=0.d0   
           !********particle numbers for bosons    
          
          ELSE
           do j=0,fl_max; DO j1=0,fl_max  
           IF(j>j1) then; w_hop(j,j1)=sqrt(j*(fl_max-j+1.d0)); 
           else if(j1>j) then; w_hop(j,j1)=sqrt((j+1.d0)*(fl_max-j)); 
           endif               
           enddo; enddo 
           !********tab. matr. el. for the worm (spin) 
  
           do j=0,fl_max; Pat(j)=j*1.d0-fl_max*0.5d0; enddo; 	
           !********magnetization for spins              
          
          ENDIF 
          Pat(-3)=0.d0

      IF(triangular==0) then
! simple cubic lattice            
!directions 
      DO i=1,dim; back(i)=i+dim; back(i+dim)=i; ENDDO 
 
! lattice connections 
      vector=0 
      IF(dim>1) n12=Ns(1)*Ns(2) 
                      IF(dim==1) then;  n1=0 
      DO ix=1,Ns(1); n1=n1+1 
	vector(1,n1)=ix-1   
	nn(1)=n1+1; IF(ix==Ns(1)) nn(1)=n1-Ns(1)+1;   
	DO kf=1,dim;kb=back(kf);ms(kf,n1)=nn(kf);ms(kb,nn(kf))=n1;ENDDO 
      ENDDO 
	           ELSE IF(dim==2) then;  n1=0      
      DO iy=1,Ns(2); DO ix=1,Ns(1); n1=n1+1  
	vector(1,n1)=ix-1; vector(2,n1)=iy-1;
	nn(1)=n1+1;     IF(ix==Ns(1))   nn(1)=n1-Ns(1)+1 
	nn(2)=n1+Ns(1); IF(iy==Ns(2))   nn(2)=n1-n12+Ns(1) 
	DO kf=1,dim;kb=back(kf); ms(kf,n1)=nn(kf); ms(kb,nn(kf))=n1;ENDDO  
	ENDDO; ENDDO 
                 ELSE IF(dim==3) then;  n1=0      
      DO iz=1,Ns(3); DO iy=1,Ns(2); DO ix=1,Ns(1); n1=n1+1;
	vector(1,n1)=ix-1; vector(2,n1)=iy-1;   vector(3,n1)=iz-1;
	nn(1)=n1+1;     IF(ix==Ns(1))   nn(1)=n1-Ns(1)+1 
	nn(2)=n1+Ns(1); IF(iy==Ns(2))   nn(2)=n1-n12+Ns(1) 
	nn(3)=n1+n12;   IF(iz==Ns(3))   nn(3)=n1-Nsite+n12 
      DO kf=1,dim;kb=back(kf); ms(kf,n1)=nn(kf); ms(kb,nn(kf))=n1;ENDDO  
	ENDDO; ENDDO; ENDDO  
                 ELSE; stop ' I did not plan for four dimensions' 
                 ENDIF    

       ELSE
!triangular lattice in dim=2 only
       IF(dim.ne.2) stop 'triangular lattice is for dim=2 only'
!directions
       do i=1,3; back(i)=i+3; back(i+3)=i; enddo
       n12=Ns(1)*Ns(2); IF(n12.ne.Nsite) stop 'n12 is not Nsite'
!lattice connections
      kf=0
      do iy=1,Ns(2);   do ix=1,Ns(1); kf=kf+1      
	kb=kf+1; IF(ix==Ns(1)) kb=kb-Ns(1)
      ms(1,kf)=kb; ms(back(1),kb)=kf
	kb=kf+Ns(1); IF(iy==Ns(2)) kb=kb-n12
      ms(2,kf)=kb; ms(back(2),kb)=kf 
      enddo; enddo

      kf=0
      do iy=1,Ns(2);   do ix=1,Ns(1); kf=kf+1      
	kb=ms(4,ms(2,kf)) 
	ms(3,kf)=kb; ms(back(3),kb)=kf 
      enddo; enddo

       ENDIF


! Prepare the interaction matrix between sites 
!  zz, xy  coupling 
      kcf=1;  kcb=back(kcf)                  !*** chain direction 
      do ir=1,Nsite; do kf=1,dir
          IF(kf.ne.kcf .AND. kf.ne.kcb) then !*** perp to chains
          Int_z(ir,kf)=Jz2; Int_xy(ir,kf)=s_hop2 
          ELSE                               !*** in chain direction
          Int_z(ir,kf)=Jz1; Int_xy(ir,kf)=s_hop1
          endif
      enddo; enddo    
          
      do ir=1,Nsite;  IF(mod(ir,2)==0) then  !*** Peierls coupling
      irr=ms(kcb,ir)
      Int_z(ir,kcb)=Int_z(ir,kcb)+dJz; 
      Int_z(irr,kcf)=Int_z(ir,kcb)
      Int_xy(ir,kcb)=Int_xy(ir,kcb)+dJxy;
      Int_xy(irr,kcf)=Int_xy(ir,kcb)      
      endif; enddo    
                   
      Hsite=H_0; Esp=0.d0      
      IF(spin_or_boson==2) then 
        do ir=1,Nsite; do kf=1,dir 
        Hsite(ir)=Hsite(ir)+Int_z(ir,kf)*fl_max*0.5d0
        Esp=Esp+Int_z(ir,kf)*fl_max**2*0.125d0 
        enddo; enddo 
        !****** chem potentials when converting spins to particles
      ENDIF    
         

!onsite disorder (if any)
      rwall=Ns(1)/2-0.5d0; rx=0.d0; ry=0.d0
      do ir=1,Nsite
      IF(close_open==0) then
      Hsite(ir)=Hsite(ir)+W(1)*(2.d0*rndm()-1.d0) ! random potential 
      else
                 rx=vector(1,ir)-rwall
                 ry=vector(2,ir)-rwall
      
!	radius=rx*rx+ry*ry+1.d0
!!	IF(radius>1.*1.*1.) then
!      wmso(ir)=wmso(ir)*dexp(-0.d0*dlog(radius))     ! worm reweighing in the trap
!!	ELSE
!!	ENDIF

      Hsite(ir)=Hsite(ir)-W(1)*rx*rx-W(2)*ry*ry              
! parabolic potential 
      IF(rx<-rwall+1.d-1 .or. rx>rwall-1.d-1) Hsite(ir)=-1.d200
      IF(ry<-rwall+1.d-1 .or. ry>rwall-1.d-1) Hsite(ir)=-1.d200 
                                                    ! wall to kill periodic 
	                                              ! boundary conditions
      ENDIF     
      enddo 

      print*, 'Hamiltonian done' 
                        
	    END SUBROUTINE hamiltonian           
 
!_______configuration and statistics initialization_____________
!_______(no-kinks between the sites)____________________________
	SUBROUTINE initial
	IMPLICIT NONE       
	INTEGER :: Npat, Ncnt, Nholes, ix,k,iy,i
	DOUBLE PRECISION :: tplus, tminus,t0,b1 
	 
!initialize the linked list arrays
        do i=1,kink_max; nlist(i)=i;  numnam(i)=i; enddo 
	present=.FALSE.; ira=0; masha=0; Rim=0.d0; tim=0.d0
	nmnm=Nsite; kinks=0; arr=0;
 
!making the elements
	DO j=1,Nsite
	fl(j)=0; site(j)=j; tau(j)=beta; time(j)=0.d0
	mk(j)=0; mr(j)=j; ml(j)=j ;   ENDDO

!putting holes (if any)  
        Nholes=Nsite*hole_conc; print*, 'Number of holes =',Nholes 
        IF(Nholes.ne.0) then;        
          j=0; do1: DO ; ix=oneoutof(Nsite)
               IF(fl(ix).ne.-3) then; fl(ix)=-3; j=j+1; ENDIF
               IF(j==Nholes) EXIT do1; 
               ENDDO do1      
        ENDIF
                   
!putting particles 
        IF(spin_or_boson==1)then;Npat=Mz                ! for bosons 
                            else;Npat=Mz+fl_max*Nsite/2 ! for spins
        ENDIF                            
   
        IF(Npat.ne.0) then;        
               j=0; 
               do2: DO ; ix=oneoutof(Nsite)
               IF(fl(ix).ne.-3 .AND. fl(ix) < fl_max) THEN 
               fl(ix)=fl(ix)+1; j=j+1; ENDIF
               IF(j==Npat) EXIT do2; 
               ENDDO do2      
                
!               DO iy=1,Ns(1); DO ix=1,Ns(1); i=(iy-1)*Ns(1)+ix
!               j=mod(iy-1,3); j=mod(ix-1+j,3)
!               IF(j==0) fl(i)=fl(i)+1
!               ENDDO; ENDDO                  ! 1/3 filling?


        ENDIF                    
!particle counter initialized	         
	do j=1,nmnm; prtcls=prtcls+Pat(fl(j)); enddo
	prtcls=prtcls*beta
	
!statistics initialized	
	step=0.d0; stat=1.d-15; Energy=0.d0; Upot=0.d0 ; Ekin=0.d0; 
        StatSq=1.d-14; Sq=0.d0
        StatSqe=1.d-14; Sqe=0.d0; UpDown=0.d0
        prntout=0
        Sq_prn=0.d0; rs_prn=0.d0; E_prn=0.d0; N_prn=0.d0
        Sqe_prn=0.d0; UD_prn=0.d0; histUD=0.d0
        zjj=1.d-15; jjcorr=0.d0

        Nav=0.d0; Wind=0.d0
        nsquare=0.d0; stat2=1.d-15            
        aa=0.d0; bb=0.d0; km=30 
       
	  gr=0.d0; grnm=0.d0; corr_nn=0.d0
	  
	dtgr=inv_temp/(Nmax*2.d0*Nmax)          !** sqrt dense intervals 	  
        DO k=-Nmax+1,Nmax; 
        IF(k>=0) then; timegr(k)=k*dtgr*k; 
        else; timegr(2*Nmax+k)=-k*dtgr*k+inv_temp; endif 
        ENDDO;  timegr(2*Nmax)=inv_temp

!	dtgr=inv_temp/(Nmax*2.d0)               !** equal step intervals        
!       DO k=0,2*Nmax; timegr(k)=k*dtgr; ENDDO;        

!granulated time/energy
        dtq=dtgr*t_quanta
 	twinq=twin*t_quanta        
        timeq=timegr*t_quanta
        tc_win=tcc_win*t_quanta
        
!A tricky way to collect statistics for <nn> by using an extended interval
        DO k=0,2*Nmax-1; ttt(k)=timeq(k);
        ttt(k+2*Nmax)=timeq(k)+beta;
        ttt(k+4*Nmax)=timeq(k)+beta+beta;        
        enddo; ttt(6*Nmax)=beta+beta+beta; ttt(6*Nmax+1)=5*beta
 

!if you wish to continue, then read what was prepared before 
        IF(old_new1==1) then; call RD_cnf;  p_ther=0.; endif
        IF(old_new2==1) then; call RD_stat; p_ther=0.; endif


           do iy=1,Ns(2); do ix=1,Ns(1);
           i=ix+(iy-1)*Ns(1);
           b1=(ix-1)+0.5*(iy-1); b1=b1*piSq
           ccS(i)=dcos(b1);
           ssS(i)=dsin(b1);
           enddo; enddo

           binN=0.d0
           do ix=1,Nsite
           iy=1+((Nsite*1.d0)/(ix*1.d0)-1.d0)/UDbin
           IF(iy<101 .AND. iy>0) then
           binN(iy)=binN(iy)+1.d0
           binN(0)=binN(0)+1.d0
           endif
           enddo

        print*, 'initial done'
  	call testing
  	
        END SUBROUTINE initial
        
!________________________________________________________________
! VERY IMPORTANT PROCEDURE OF PERFORMING MEASUREMENTS 
!________________________________________________________________
       SUBROUTINE measure(i1)
       INTEGER, INTENT(IN) :: i1 
       integer ::  kp,kw,k, i,name,name1, nl,jgr,kat,here, katN
        integer :: fi,fm,ix,iy 
       integer, parameter ::    kat0=0
       double precision :: kpa,tleft,tright,addt,ht,ddt,pp
       double precision :: b1,b2,tnorm,cop,nm(0:fl_max)
       double precision :: ndistr(Nsite)
       double precision :: Mup, Mdown
       
       katN=2*Nmax       
       IF(.not.present) then           ! NO-WORM STATISTICS
               
       kpa=prtcls/beta; kp=NINT(kpa)   ! PARTICLE NUMBERS
       IF(CANONICAL .AND. kpa.ne.Nsite) RETURN     ! FOR CANONICAL calc
       stat=stat+1.d0                  ! PARTITION FUNCTION    
       IF(kp>Mpat.or.kp<-Mpat)         stop 'kp > Mpat; increase Npat' 
       Nav(kp)=Nav(kp)+1               ! PARTICLE NUMBERS HISTOGRAM   
       Nav(Mpat+1)=Nav(Mpat+1)+kp;     ! PARTICLE NUMBER AVERAGE 
       
	 Ekin=Ekin-kinks               ! kinetic energy 
	 Upot=Upot+Udiag               ! Potential energy 
	 Energy=Energy+Udiag-kinks     ! Total energy 
	                    
                                         ! winding numbers
       DO j=1,dim; k=back(j); 
          IF(triangular==0) then         ! square lattice
          kw=(arr(j)-arr(k));
          ELSE                           ! triangular lattice
          IF(j==1) then
          kw=(arr(j)-arr(k)+arr(6)-arr(3))
          else
          kw=(arr(j)-arr(k)+arr(3)-arr(6))
          ENDIF
          ENDIF
!diagnostics       
       i=kw/Ns(j); i=i*Ns(j)        
       IF(i-kw.ne.0) then
            print*, 'winding numbers not integer in measure' 
	    print*, j, Ns(j), kw, arr(j),arr(k)  
            print*, ' ira-masha are now at '
            print*,  Rim, tim/beta 
       stop   ; endif
       IF(j==1) kwcheck=kw 
       kw=kw/Ns(j)
       IF(ABS(kw)>Mwind) stop 'kw > Mwind; increase Mwind' 	       
       Wind(j,kw)=Wind(j,kw)+1; 
!       IF(kw.ne.0) then
!       print*, j, Ns(j), arr(j), arr(k)
!       print*, kw        
!       pause
!       endif    
       END DO   !WINDING NUMBERS HISTOGRAM

      kw=Nsite*jc
      IF(mod(i1,kw)==0) call sigma

! <n*n>
        go to 543
        IF(mod(i1,100*Nsite)==0) then          !
        stat2=stat2+1.d0; kpa=0.d0; nm=0.d0
        do i=1,nmnm; name=nlist(i); k=Pat(fl(name))
        kpa=kpa+tau(name)*k**2
        nm(k)=nm(k)+tau(name) 
        enddo; kpa=kpa/beta/Nsite        ! PARTICLE NUMBER SQUARE
        nsquare=nsquare+kpa              ! PARTICLE NUMBER SQUARE AVERAGE
        nm=nm/beta/Nsite
        do k=0,fl_max; grnm(Ns(1),k,k)=grnm(Ns(1),k,k)+nm(k); 
        enddo          ! stat. of occ. num. fluctuations    
        endif

! n-n correlation function
         go to 543 
         IF(mod(i1,Ns(1)*Ns(1)*20)==0) then;  
         StatSq=StatSq+1.d0     
         StatSqe=StatSqe+1.d0;
         call Equal_time_nn

         ndistr=0.d0
         do i=1,nmnm; name=nlist(i); k=site(name)
         ndistr(k)=ndistr(k)+ Pat(fl(name))*tau(name); enddo
         ndistr=ndistr/beta

         open(1,file='ndistrW.dat')
         do i=1,Nsite
         write(1,*) i, ndistr(i)
         enddo
         close(1)
       
        kpa=0.d0; Mup=0.d0; Mdown=0.d0 
        do i=1,Nsite; kpa=kpa+ndistr(i)*ccS(i); 
        IF(ndistr(i)>0.5d0) then; Mup=Mup+1.d0
                            else; Mdown=Mdown+1.d0; endif
        enddo
        tleft=0.d0; tright=0.d0
        do i=1,Nsite; tleft=tleft+ndistr(i)*ssS(i)
                      tright=tright+ndistr(i); enddo   
        IF(ABS(tright- prtcls/beta)>1.d-8) THEN
        print*, tright,  prtcls/beta
        stop 'integrated density .ne. N in static'
        ENDIF

        Sq=Sq+(kpa**2+tleft**2)/(Nsite**2) 
         
        IF(Mdown>0.d0) then; UpDown=UpDown+Mup/Mdown
                       else; UpDown=UpDown+Nsite; endif 
        i=1+(Mup/Mdown)/UDbin
        IF(i .LE. 0) then; 
        print*, i, mup, Mdown, UDbin; 
        stop 'bin number .LE. 0'; ENDIF
        IF(i < 101 ) then; histUD(i)=histUD(i)+1.d0
                           histUD(0)=histUD(0)+1.d0; ENDIF        

       ENDIF 
      
 543    continue

!       IF(mod(i1,Nsite*100)==0) call NN_0       !    40 is OK,  n-n correlator
      
!Green function, WORM IS PRESENT              

       else IF(mod(i1,10)==0) THEN;      !    10 is OK          
       IF(tim==0 .or. tim==beta) then;  jgr=0;
       else IF(tim<=beta2) then; jgr=sqrt(tim/dtq);
       else; jgr=sqrt((beta-tim)/dtq); jgr=2*Nmax-jgr-1; endif
       fi=fl(ira)+1; fm=fl(masha) 
!       jgr=tim/dtq                      ! for linear grid 
 
!  looking at ira 
             nl=ml(ira);tleft=tau(nl);tright=tau(ira);
             here=site(ira)   
             kat=jgr  
             ddt=timeq(kat)-tim; IF(ddt>0) ddt=ddt-beta	    
             addt=ABS(ddt)
             aa: DO; IF(addt>tleft .or. addt>twinq) exit aa                
             b2=min(addt+tright,twinq)
             b1=min(tleft-addt,twinq)
             tnorm=t_quanta/(b1+b2) 
                            
                   IF(onsite_nn==1) then; pot=0.d0; else 
                   call potential_move(nl,ira,ddt)
                   ENDIF
 
                  pp=(-Hsite(here)*ddt+pot)*(Fss(fl(nl))-Fss(fl(ira)))
                  pp=pp+U0*ddt*(Fs2(fl(nl))-Fs2(fl(ira)))
                  cop=tnorm*EFUN(-pp)                  
!!                  gr(Rim(1),kat)=gr(Rim(1),kat)+cop 
                  gr0(kat)=gr0(kat)+cop             

             IF(kat==0) then                   
!             grnm(Rim(1),fi,fm)=grnm(Rim(1),fi,fm)+cop 
             IF(Rim(1).ne.0) then
!!             gr(Rim(1),katN)=gr(Rim(1),katN)+cop          
             gr0(katN)=gr0(katN)+cop
             endif 
             ENDIF
                  		      
             kat=kat-1; IF(kat<0) kat=2*Nmax-1
             ddt=timeq(kat)-tim; IF(ddt>0) ddt=ddt-beta	             
             addt=ABS(ddt)                  
             END DO aa

             kat=jgr+1;  IF(kat>2*Nmax) kat=1
             ddt=timeq(kat)-tim; IF(ddt<0) ddt=ddt+beta		              
             bb: DO; IF(ddt>tright .or. ddt>twinq) exit bb 
             b2=min(tright-ddt,twinq)
             b1=min(tleft+ddt,twinq)
             tnorm=t_quanta/(b1+b2) 
                                      
           
                 IF(onsite_nn==1) then; pot=0.d0; else
                 call potential_move(nl,ira,ddt)
                 ENDIF 
                 
                pp=(-Hsite(here)*ddt+pot)*(Fss(fl(nl))-Fss(fl(ira)))
                pp=pp+U0*ddt*(Fs2(fl(nl))-Fs2(fl(ira)))
                 cop=tnorm*EFUN(-pp)
!!                gr(Rim(1),kat)=gr(Rim(1),kat)+cop
                gr0(kat)=gr0(kat)+cop
                
                IF(kat==2*Nmax) then 
!                grnm(Rim(1),fi,fm)=grnm(Rim(1),fi,fm)+cop
                IF(Rim(1).ne.0) then
!!                gr(Rim(1),kat0)=gr(Rim(1),kat0)+cop			  	
                gr0(kat0)=gr0(kat0)+cop
                endif
                ENDIF
             
             kat=kat+1;  IF(kat>2*Nmax) kat=1
             ddt=timeq(kat)-tim; IF(ddt<0) ddt=ddt+beta 
             END DO bb

!  looking at masha 
             nl=ml(masha);tleft=tau(nl);tright=tau(masha);
             here=site(masha) 

	     kat=jgr
             ddt=tim-timeq(kat); IF(ddt<0) ddt=ddt+beta         
             cc: DO; IF(ddt>tright .or. ddt>twinq ) exit cc
             b1=min(tright-ddt,twinq)
             b2=min(tleft+ddt,twinq)
             tnorm=t_quanta/(b1+b2);        
                 
              IF(onsite_nn==1) then; pot=0.d0; else            
              call potential_move(nl,masha,ddt)
              ENDIF
              
              pp=(-Hsite(here)*ddt+pot)*(Fss(fl(nl))-Fss(fl(masha)))
              pp=pp+U0*ddt*(Fs2(fl(nl))-Fs2(fl(masha)))
              cop=EFUN(-pp)*tnorm
!!              gr(Rim(1),kat)=gr(Rim(1),kat) +cop
              gr0(kat)=gr0(kat) +cop
              
              IF(kat==0) then
!              grnm(Rim(1),fi,fm)=grnm(Rim(1),fi,fm)+cop 
              IF(Rim(1).ne.0) then 
!!              gr(Rim(1),katN)=gr(Rim(1),katN)+cop                  	
              gr0(katN)=gr0(katN)+cop
              endif
              ENDIF
     
             kat=kat-1; IF(kat<0) kat=2*Nmax-1
             ddt=tim-timeq(kat); IF(ddt<0) ddt=ddt+beta  
             END DO  cc

             kat=jgr+1; IF(kat>2*Nmax) kat=1
             ddt=tim-timeq(kat); IF(ddt>0) ddt=ddt-beta
	     addt=ABS(ddt)           
             dd: DO; IF(addt>tleft .or. addt>twinq) exit dd
             b1=min(tright+addt,twinq)
             b2=min(tleft-addt,twinq)
             tnorm=t_quanta/(b1+b2)
                       
              IF(onsite_nn==1) then; pot=0.d0; else     
              call potential_move(nl,masha,ddt)
              ENDIF
              
              pp=(-Hsite(here)*ddt+pot)*(Fss(fl(nl))-Fss(fl(masha)))
              pp=pp+U0*ddt*(Fs2(fl(nl))-Fs2(fl(masha)))
              cop=EFUN(-pp)*tnorm
!!              gr(Rim(1),kat)=gr(Rim(1),kat) +cop
              gr0(kat)=gr0(kat) +cop              

             IF(kat==2*Nmax) then
!             grnm(Rim(1),fi,fm)=grnm(Rim(1),fi,fm)+cop 
             IF(Rim(1).ne.0 ) then
!!             gr(Rim(1),kat0)=gr(Rim(1),kat0) +cop 	
             gr0(kat0)=gr0(kat0) +cop
             endif
             ENDIF
  
          kat=kat+1; IF(kat>2*Nmax) kat=1
          ddt=tim-timeq(kat); IF(ddt>0) ddt=ddt-beta
	  addt=ABS(ddt)
             END DO dd 
             
	    endif
          
       END SUBROUTINE measure
!***************************************
      subroutine sigma
      implicit none
      integer :: name0, name0k, name0l
	integer :: jo, jname, k0
      double precision  :: t, factor
	double precision  :: addr, addi, addry, addiy
      
	DO jo=0,Njj                               ! frequency loop
	addr=0.d0
	addi=0.d0
	addry=0.d0
	addiy=0.d0

      DO jname=1,nmnm                           ! kink loop 
	name0=nlist(jname)                        ! name
      name0l=ml(name0)                          ! left neighbor
	IF(name0l.ne.name0) then        ! we have a kink        
!	name0k=mk(name0)                          ! kinjk association 
      t=2.d0*pi*jo*time(name0)/beta             ! phase factor

!      k0=direction(name0,name0k)      
!	factor=1.d0     
!	IF(Pat(fl(name0))>Pat(fl(name0l))) factor=-factor

!	IF(k0==1)            then; 
      addr=addr+dcos(t)
      addi=addi+dsin(t) 
!	else IF(k0==2)       then;  
!      addry=addry+factor*0.5d0*dcos(t)
!      addiy=addiy+factor*0.5d0*dsin(t)      
!	else IF(k0==3)       then; factor=-factor  
!      addr=addr+factor*0.5d0*dcos(t)
!      addi=addi+factor*0.5d0*dsin(t)
!      else;                      factor=-factor 
!      addry=addry+factor*0.5d0*dcos(t)
!      addiy=addiy+factor*0.5d0*dsin(t)      
!      ENDIF
            
      endif  ! we have a kink
      enddo  ! kink loop 
!	IF(jo==0) then
!	IF(ABS(kwcheck-addr)>1.d-6) then
!      print*, kwcheck, addr
!	pause
!	ENDIF
!	ENDIF


	jjcorr(jo)=jjcorr(jo)+addr*addr+addi*addi
!      jjcorr(jo)=jjcorr(jo)+addry*addry+addiy*addiy
	enddo  ! frequency loop   
      zjj=zjj+4.d0

      end subroutine sigma
!***********************************************


      subroutine NN_0
      implicit none
      integer :: name0, name1, name, l, jout 
      real :: t, t0, t_max
      
      DO jout=1,jc; name0=oneoutof(nmnm); name0=nlist(name0)

      t_max=tau(name0) + 2.d0*tc_win ; t0=time(name0)            
      t=tc_win ;   name1=name0
      
      do; name1=ml(name1);  t=t-tau(name1); if (t < 0.d0 ) exit
      end do      
      t=t+tau(name1)
      
      DO l=0, Ns(1)-1             
       name=name1 
       do; call cor_nn(name0,name,l); if (t > t_max) exit
       name=mr(name); t=t+tau(name)
       end do
              
       name1=ms(1,name1)
       t=time(name1)-t0; if(t>0) t=t-beta
       t=t+tc_win        
                            
       do; t=t+tau(name1); if (t > 0) exit; name1=mr(name1)
       end do
              
      END DO
            
      enddo      
      end subroutine NN_0
!_________________________________________________________________
     
      subroutine Equal_time_nn
      INTEGER :: i,k,loop,name,n2,n11(Ns(1)),cou
      DOUBLE PRECISION :: t2, b1,b2,b3,ndistr(Nsite)
      DOUBLE PRECISION :: t00, t0, tslice


      cou=1
      i=oneoutof(nmnm); name=nlist(i); 
      t00=time(name); t0=0.d0; tslice=rndm()*beta;

      e: DO; IF(t0+tau(name) .GT.  tslice ) EXIT e;
             t0=t0+tau(name); name=mr(name); enddo e 
      k=site(name); ndistr(k)=fl(name); n11(cou)=name  ! first site
      tslice=s_t(plus,t00,tslice)

      n2=name
      DO loop=1,Ns(1)-1; 
      n2=ms(1,n2); t2=-s_t(minus,tslice,time(n2))
                   a : DO; IF(t2+tau(n2)>0) exit a; 
                   t2=t2+tau(n2); n2=mr(n2);  ENDDO a
       k=site(n2); ndistr(k)=fl(n2); cou=cou+1; n11(cou)=n2
       ENDDO

       DO cou=1,Ns(1); n2=n11(cou)
       DO loop=1,Ns(2)-1;
       n2=ms(2,n2); t2=-s_t(minus,tslice,time(n2))
                   b : DO; IF(t2+tau(n2)>0) exit b;
                   t2=t2+tau(n2); n2=mr(n2); ENDDO b
       k=site(n2); ndistr(k)=fl(n2)
       ENDDO
       ENDDO

       open(1,file='ndistr.dat')
       do i=1,Nsite
       write(1,*) i,ndistr(i)
       enddo
       close(1)


       b1=0.d0; b2=0.d0; b3=0.d0
       do i=1,Nsite; b1=b1+ndistr(i)*ccS(i); 
                     b2=b2+ndistr(i)*ssS(i)
                     b3=b3+ndistr(i)
                                      enddo
       IF(ABS(b3-prtcls/beta)>1.d-8) THEN
       print*, b3,prtcls/beta 
       stop 'b3 .ne. N'
       ENDIF 
       Sqe=Sqe+(b1*b1+b2*b2)/(Nsite*1.d0*Nsite)

      end subroutine Equal_time_nn

!_________________________________________________________________
      subroutine COR_NN(name1,name2,kp)
      implicit none
      integer :: name1, name2, j, kp
      double precision :: p, delta, t1, t2, t, ta, tb, tc, d
     
      ! . . . . . . . . . . . . . . . . . . |____________| . . . . . 
      !                                     |    name2 
      !                    |<--- delta ---->|
      !.. . . . |__________|. . . . . . . . .  . . . . . . . . . . .
      !            name1              
                                              
      p=nmnm; p=Pat(fl(name1))*p*Pat(fl(name2))

      t1=tau(name1); t2=tau(name2)
 
      delta=s_t(minus,time(name2),time(mr(name1)))
      IF(delta<0) delta=delta+beta     
      IF(delta==0 .or. delta==beta) then; j=1; 
      else IF(delta<=beta2) then; j=sqrt(delta/dtq)+1;
      else; j=sqrt((beta-delta)/dtq); j=2*Nmax-j; endif 
     
      if (t2 > t1) then
         ta=delta+t1; tb=delta+t2; d=t1
      else
          ta=delta+t2; tb=delta+t1; d=t2
      end if
      
      tc=delta+t1+t2
      
      t=ttt(j)  
      do; if (t > ta) exit
          corr_nn(kp,j)=corr_nn(kp,j)+p*(t-delta)
          j=j+1; t=ttt(j)
      end do
      
           
      do; if (t > tb) exit
          corr_nn(kp,j)=corr_nn(kp,j)+p*d
          j=j+1; t=ttt(j)
      end do
      
       do; if (t > tc) exit
          corr_nn(kp,j)=corr_nn(kp,j)+p*(delta+t1+t2-t)
          j=j+1; t=ttt(j)
      end do
       
      end subroutine COR_NN           
!______________________________________________________________

!________PRINTING______________________________________________
         SUBROUTINE printing(i1)  
         DOUBLE PRECISION :: xx,yy,qq,xxn,yyn,Ecan
         DOUBLE PRECISION ::lambda,kappa,deltaN,zu, remn 
!         REAL*4 :: result(11)
	 INTEGER, INTENT(IN) :: i1
	 INTEGER :: jo   

       print*; print* 
       print*, 'MC steps (in millions)=', i1/1000000
       print*, 'stat=', stat
	 print*, 'zjj =', zjj 
       print*, '<n^2> stat       =', AINT(stat2)
       print*, 'StatSq         =', StatSq
!       print*, 'StatSqe        =', StatSqe
       print*, '  nmnm=', nmnm, 'present=', present
!       print*, 'kinks  =', kinks,'  arrows=', (arr(j), j=1,dir)
       print*, 'particles now =',prtcls/beta
!       print*, 'Udiag now     =',Udiag/inv_temp
       
	              IF(p_ther==0 . AND. stat.ne.0.d0) then
       prntout=prntout+1
       result(1)=inv_temp
       print*, 'K        =', Ekin/(stat)    ! Ekin/(stat*inv_temp)
!       print*, 'U        =', Upot/(stat*inv_temp)
!       print*, 'S corr.  =', Esp*t_quanta
!       print*, 'E        =', (Energy/(stat*inv_temp)+Esp*t_quanta)/Nsite 
       xxn=Nav(Mpat+1)/stat
       yyn=xxn+Nsite*fl_max*0.5d0 
!       print*, 'particles (now) =', prtcls/beta
!       print*, '   <M> =', xxn/(Nsite*0.5d0*fl_max)
       print*
       N_prn(prntout)=xxn
       call ERSTAT(N_prn,prntout,amax,tmax,amin,tmin)
       print 705, N_prn(prntout), (amax-amin)/2.
       remn=(amax-amin)/2.
 705   format(6x,'<N>         =',g12.5,4x,'+-',g10.3)
!      print 706, amax, tmax
!      print 706, amin, tmin
       result(8)= N_prn(prntout); result(9)=(amax-amin)/2.
       print*

       Ecan=Energy/(stat*inv_temp)+(mag_field-onsite)*xxn
       E_prn(prntout)=Ecan/Nsite
       call ERSTAT(E_prn,prntout,amax,tmax,amin,tmin)
       print 702, E_prn(prntout), (amax-amin)/2.
 702   format(6x,'E_grand+mu*<N> =',g12.5,4x,'+-',g10.3)
!      print 706, amax, tmax
!      print 706, amin, tmin
 706   format(1x, 12x, g11.4, ' at -> ', g9.2)
       result(2)= E_prn(prntout); result(3)=(amax-amin)/2.
       print*


       Sq_prn(prntout)=sqrt(Sq/StatSq) 
       call ERSTAT(Sq_prn,prntout,amax,tmax,amin,tmin)
       print 703, Sq_prn(prntout),( amax-amin)/2.
 703   format(6x,'Sq(static)     =',g12.5,4x,'+-',g10.3)
!      print 706, amax, tmax
!      print 706, amin, tmin
       result(4)= Sq_prn(prntout); result(5)=(amax-amin)/2.
       print*

       Sqe_prn(prntout)=sqrt(Sqe/StatSqe)
       call ERSTAT(Sqe_prn,prntout,amax,tmax,amin,tmin)
       print 707, Sqe_prn(prntout),( amax-amin)/2.
 707   format(6x,'Sq(same time)  =',g12.5,4x,'+-',g10.3)
!      print 706, amax, tmax
!      print 706, amin, tmin
       result(10)= Sqe_prn(prntout); result(11)=(amax-amin)/2.
       print*

       UD_prn(prntout)=UpDown/StatSq
       call ERSTAT(UD_prn,prntout,amax,tmax,amin,tmin)
       print 709, UD_prn(prntout),( amax-amin)/2.
 709   format(6x,'N(up)/N(down)  =',g12.5,4x,'+-',g10.3)
!      print 706, amax, tmax
!      print 706, amin, tmin
       print*

       IF(dim==1) then
       xx=1.d-14;do j=-Mpat,Mpat; xx=xx+Nav(j); enddo
       kappa=0.d0;do j=-Mpat,Mpat; zu=j 
       kappa=kappa+(Nav(j)*(zu-xxn)**2)/xx; 
       enddo 
       print*, 'kappa =',kappa*inv_temp/Nsite        
       xx=1.d-14;do j=-Mwind,Mwind; xx=xx+Wind(1,j);enddo   
       lambda=0.d0; do j=-Mwind,Mwind; zu=j
       lambda=lambda+(Wind(1,j)*zu**2)/xx;
       enddo
       print*, 'lambda =',lambda*Nsite/inv_temp
       xx=sqrt(kappa*lambda+1.d-15)
       print*, 'K      =',pi*xx 
              
       ELSE IF(dim==2) then
       xx=1.d-14;do j=-Mwind,Mwind; xx=xx+Wind(1,j);enddo
       yy=1.d-14;do j=-Mwind,Mwind; yy=yy+Wind(2,j);enddo      
       lambda=0.d0; do j=-Mwind,Mwind; zu=j
       lambda=lambda+(Wind(1,j)*zu**2)/xx;
       lambda=lambda+(Wind(2,j)*zu**2)/yy;
       enddo

       rs_prn(prntout)=lambda/2.d0 
       call ERSTAT(rs_prn,prntout,amax,tmax,amin,tmin)
       print 704, rs_prn(prntout),( amax-amin)/2.
 704   format(6x,'r_s/T=<W_x*W_x> =',g12.5,4x,'+-',g10.3)
!      print 706, amax, tmax
!      print 706, amin, tmin
       result(6)= rs_prn(prntout); result(7)=(amax-amin)/2.

       xx=1.d-14;do j=-Mpat,Mpat; xx=xx+Nav(j); enddo
       deltaN=0.d0;do j=-Mpat,Mpat; zu=j
       deltaN=deltaN+(Nav(j)*(zu-xxn)**2)/xx;
       enddo       
!	xx=(1.+2./Nsite)/12.
       print*, '<dN*dN>=',deltaN
       PRINT*
       print*, 'Lambda_S=',rs_prn(prntout)/inv_temp
	 print*, 'kappa=', deltaN*inv_temp/(1.d0*Nsite) 
	xx=(rs_prn(prntout)/inv_temp)
	xx=xx/(deltaN*inv_temp/(1.d0*Nsite)) 
	 print*, 'sound velocity=', sqrt(xx)
       print* 

       ELSE 
       xx=1.d-14;do j=-Mpat,Mpat; xx=xx+Nav(j); enddo
       deltaN=1.d-14;do j=-Mpat,Mpat; zu=j 
       deltaN=deltaN+(Nav(j)*(zu-xxn)**2); 
       enddo; deltaN=deltaN/xx
       print*, '<dN*dN> =',deltaN
       
       xx=1.d-14;do j=-Mwind,Mwind; xx=xx+Wind(1,j);enddo
       yy=1.d-14;do j=-Mwind,Mwind; yy=yy+Wind(2,j);enddo      
       qq=1.d-14;do j=-Mwind,Mwind; qq=qq+Wind(3,j);enddo  
       lambda=0.d0; do j=-Mwind,Mwind; zu=j
       lambda=lambda+(Wind(1,j)*zu**2)/xx;
       lambda=lambda+(Wind(2,j)*zu**2)/yy;
       lambda=lambda+(Wind(3,j)*zu**2)/qq;
       enddo; lambda=lambda/dim 

       rs_prn(prntout)=lambda
       call ERSTAT(rs_prn,prntout,amax,tmax,amin,tmin)
       print 717, 1/inv_temp, rs_prn(prntout),( amax-amin)/2.
 717   format(6x,'r_s*L/T=<W*W>/3 =',g12.7,4x,g12.5,4x,g10.3)
      print 706, amax, tmax
      print 706, amin, tmin
       result(6)= rs_prn(prntout); result(7)=(amax-amin)/2.
       print*

       open(1,file='fresult.dat')
	write(1,111) mag0, N_prn(prntout)/Nsite, remn/Nsite, 
     .rs_prn(prntout),(amax-amin)/2., deltaN
	 close(1)
 111   format(g12.5,g12.7,g10.3,1x,g12.7,g10.3,1x,g12.7,g10.3)
       ENDIF
                   ENDIF

!       print*
!       print*
!       print*, ' ira-masha are now at '
!       print*,  Rim, tim/beta

       print*
       open(1,file='result.dat')
       write(1,810) result
       close(1)       
 810   format(g10.5,5(3x,g10.5,1x,g8.2))

       open(1, file='12h.dat')
       do j=1,100
       IF(binN(j)>0.1) then
       write(1,*) (j-0.5d0)*UDbin, histUD(j)/(histUD(0)*binN(j))
       endif
       enddo
       close(1)


!       open(1, file='JJcorr.dat')
!       do j=1,100
!       xx=(j-0.5d0)*inv_temp/100.d0
!	 yy=jjcorr(0)
!	 do jo=1,Njj
!	 yy=yy+2.d0*jjcorr(jo)*dcos(xx*2.d0*pi*jo/inv_temp)
!      enddo
!       write(1,*) xx, yy/((inv_temp**2)*Nsite*zjj)
!      enddo
!      close(1)      
!      print*, 'CORRELATOR w=0', jjcorr(0)/(Nsite*zjj)

       open(1, file='KKWcorr.dat')
         write(1,*) 'MC steps (in millions)=', i1/1000000
	 do jo=0,Njj
	 xx=2.d0*pi*jo/inv_temp
	 yy=-Ekin/(Nsite*stat*inv_temp)
!	 yy= yy- jjcorr(jo)/(Nsite*zjj*inv_temp)
       yy= jjcorr(jo)/(Nsite*zjj*inv_temp)
       write(1,*) xx, yy, yy+Ekin/(Nsite*stat*inv_temp)
       enddo
       close(1)      
        
!       open(1, file='JJWTcorr.dat')
!	 do jo=1,Njj
!	 xx=2.d0*pi*jo/inv_temp
!       yy=-Ekin/(Nsite*2.d0*stat*inv_temp)
!!	 yy= yy- jjcorr(jo)/(Nsite*zjj*inv_temp)
!       yy= jjcorr(jo)/(Nsite*zjj*inv_temp)
!       write(1,*) xx*inv_temp, yy/xx
!       enddo
!       close(1)   


!         call testing

       END SUBROUTINE printing
!********************************************************
      subroutine ERSTAT(a,n,amax,tmax,amin,tmin)
!     Analizing 3/4 print-out
         
      double precision :: a, amax, tmax, amin, tmin, aa
      integer n, i
      dimension a(n)
      
      amax=-1.d200; amin=1.d200
             
      DO i=n/4+1, n;  aa=a(i)
         if (aa > amax) then
            amax=aa; tmax=i
         end if
         if (aa < amin) then
            amin=aa; tmin=i
         end if
      END DO
             
      tmax=tmax/n; tmin=tmin/n
      end subroutine ERSTAT
 
           
           
           
           
           
           
           
           
           END 
