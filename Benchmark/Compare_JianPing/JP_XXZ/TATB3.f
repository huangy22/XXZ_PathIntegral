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
        INTEGER, allocatable :: rdir(:)    ! real neigbours
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
   	INTEGER                        :: jc
	                                ! times where G and dens. -dens.
	                                ! correlators are calculateded  
        DOUBLE PRECISION :: dtgr,dtq,twin,twinq    
	DOUBLE PRECISION :: tcc_win 
	DOUBLE PRECISION :: tc_win
	                                ! windows to collect statistics 
	                                ! for G and dens.-dens.
	                      
!world-line diagram description and linking          
	INTEGER, allocatable          :: fl(:), Lat(:)    
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
	INTEGER, allocatable          :: back(:,:) 
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
         DOUBLE PRECISION :: vn
                          
                          ! (hopping) perp to chains
        DOUBLE PRECISION :: mag_field, random_field
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
	DOUBLE PRECISION, allocatable :: Int_z(:,:),Int_xy(:,:,:)
	DOUBLE PRECISION :: s_hop1,s_hop2,H_0,W
	DOUBLE PRECISION, allocatable :: Hsite(:)	 
	DOUBLE PRECISION :: dJz, dJxy
	DOUBLE PRECISION :: Esp
	DOUBLE PRECISION :: wmso
	DOUBLE PRECISION, allocatable :: Fss(:), Fs2(:)
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
!        INTEGER :: igrow,igrow1,igrow1a,igrow1b
!        INTEGER :: igrow2,igrow2a,igrow2b
!!        INTEGER :: igrow,igrow1,igrow1a,igrow1b
!        INTEGER :: igrow2,igrow2a,igrow2b
!        INTEGER :: igrow3,igrow3a,igrow3b
!        INTEGER :: igrow3,igrow3a,igrow3b
	DOUBLE PRECISION :: p_ther
	INTEGER          :: p_mes, p_prn, p_wr,start
	                  ! thermolization, measuring, printing and writing 
	                  ! in the number of updates
        INTEGER :: aa(0:30),bb(0:30) 
	DOUBLE PRECISION :: prtcls, Udiag
        DOUBLE PRECISION :: all_phi
	DOUBLE PRECISION :: Upot, Ekin , Energy, Energy2   ! energies
	DOUBLE PRECISION :: step,sum_se            ! MC update number
	DOUBLE PRECISION :: stat,stat_se           ! partition function
        DOUBLE PRECISION :: StatSqe          ! equal_time structure factor
	DOUBLE PRECISION, allocatable :: Nav(:)     ! particles
	DOUBLE PRECISION :: nsquare, stat2
        DOUBLE PRECISION, allocatable ::  ccS(:), ssS(:) !tabulated cos and sin

	                  ! on-site density square and its statistics   
	DOUBLE PRECISION, allocatable :: arr(:),arrx(:,:),arry(:,:)
        DOUBLE PRECISION, allocatable :: arrz(:,:)       
	DOUBLE PRECISION, allocatable :: Wind(:,:)
!       DOUBLE PRECISION, allocatable :: Sq(:),Sqe(:)
	                               ! winding numbers

        INTEGER, PARAMETER :: prnt_max=1000000
        INTEGER :: prntout
        DOUBLE PRECISION :: amin, amax, tmin, tmax 
        DOUBLE PRECISION :: Sq_prn(prnt_max)
        DOUBLE PRECISION :: rs_prn(prnt_max)
        DOUBLE PRECISION :: E_prn(prnt_max),Egrand_prn(prnt_max)
        DOUBLE PRECISION :: Egrand2_prn(prnt_max) 
        DOUBLE PRECISION :: N_prn(prnt_max)
        DOUBLE PRECISION :: N0_prn(prnt_max)
        DOUBLE PRECISION :: N1_prn(prnt_max)
        DOUBLE PRECISION :: N2_prn(prnt_max)
        DOUBLE PRECISION :: N3_prn(prnt_max)
        DOUBLE PRECISION :: N4_prn(prnt_max)
        DOUBLE PRECISION :: N5_prn(prnt_max)
        DOUBLE PRECISION :: N6_prn(prnt_max)
        DOUBLE PRECISION :: N7_prn(prnt_max)
        DOUBLE PRECISION :: N8_prn(prnt_max)
        DOUBLE PRECISION :: N9_prn(prnt_max)
        DOUBLE PRECISION :: N10_prn(prnt_max)
        DOUBLE PRECISION :: N11_prn(prnt_max)
        DOUBLE PRECISION :: N12_prn(prnt_max) 
        DOUBLE PRECISION :: Sqe(10000),Sqe2(10000)
        DOUBLE PRECISION :: ddc1(0:32,0:32,0:32)
        DOUBLE PRECISION :: ddc2(0:32,0:32,0:32)
        DOUBLE PRECISION :: ddc3(0:32,0:32,0:32)
        DOUBLE PRECISION :: ddc_stat(0:32,0:32,0:32)
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
        integer :: i,i1
        DOUBLE PRECISION :: ist
        REAL*4 :: result(33)
     
       dim=3
       allocate(Ns(dim))
       open(16,file='inp11')
       read(16,*) Ns(1),ij,kl
       close(16)
      
      mag_field=3.0d0
      !mag_field=0.0d0
      inv_temp=1.0d0/0.10d0
      Ns(2)=Ns(1)
      Ns(3)=Ns(1)
      xy_coupling1=0.50d0
      xy_coupling2=xy_coupling1

      OPEN(1,file='parwv')
      READ(1,*)  spin_or_boson   ; print*, 'spin/boson  ', spin_or_boson 
      READ(1,*)  onsite_nn       ; print*, 'onsite/nn   ', onsite_nn
      READ(1,*)  triangular      ; print*, 'lattice type', triangular 

      dir=6
   
      Nsite=1
      do i=1,dim 
      Nsite=Nsite*Ns(i)
      end do
      Nsite=Nsite*4

      allocate(rdir(Nsite),Lat(Nsite)) 
      READ(1,*)  Nmax; print*, 'Nmax        ', Nmax
      allocate(gr(0:Ns(1)-1,0:2*Nmax),gr0(0:2*Nmax))
      allocate(timegr(0:2*Nmax),timeq(0:2*Nmax))
      allocate(corr_nn(0:Ns(1)-1,0:6*Nmax))
      allocate(ttt(0:6*Nmax+1))      
      READ(1,*)  fl_max          ; print*, 'fl_max      ', fl_max
      allocate(grnm(0:Ns(1),0:fl_max,0:fl_max))              
      READ(1,*)  kink_max        ; print*, 'kink_max    ', kink_max
      allocate(fl(kink_max),site(kink_max))
      allocate(tau(kink_max),time(kink_max))  
      allocate(ml(kink_max),mr(kink_max),mk(kink_max))
      allocate(ms(dir,kink_max),back(dir,kink_max),Rim(dim))
      allocate(nlist(kink_max),numnam(kink_max))  
      allocate(Fss(-3:fl_max),Fs2(-3:fl_max))
      allocate(w_hop(-3:fl_max,-3:fl_max))
      allocate(Pat(-3:fl_max),s_h(-3:fl_max,-3:fl_max))
      allocate(Int_z(kink_max,dir),Int_xy(2,kink_max,dir))
      allocate(Hsite(kink_max)) 
      allocate(ccS(kink_max),ssS(kink_max))
      READ(1,*)  Mz              ; print*, 'Mz          ', Mz
      READ(1,*)  Mzmark          ; print*, 'Mzmark      ',Mzmark
      READ(1,*)  CANONICAL       
      READ(1,*)  onsite          ; print*, 'onsite      ', onsite
      READ(1,*)  vn              ; print*, 'vortex number  ', vn
      READ(1,*)  z_coupling1     ; print*, 'Jz1         ', z_coupling1
!      READ(1,*)  xy_coupling1    ; print*, 'Jx1         ', xy_coupling1
      READ(1,*)  z_coupling2     ; print*, 'Jz2         ', z_coupling2
!     READ(1,*)  xy_coupling2    ; print*, 'Jx2         ', xy_coupling2
      READ(1,*)  peierls_z       ; print*, 'Peierls_z   ', peierls_z
      READ(1,*)  peierls_xy      ; print*, 'Peierls_xy  ', peierls_xy  
      READ(1,*)  hole_conc       ; print*, 'Hole conc.  ', hole_conc
!      READ(1,*)  mag_field       ; print*, 'H           ', mag_field
      READ(1,*)  random_field    ; print*, 'W           ', random_field
      READ(1,*)  close_open      ; print*, 'close/open  ', close_open
      READ(1,*)  t_quanta        ; print*, 't_quanta    ', t_quanta 
      READ(1,*)  source          ; print*, 'source      ', source 
      READ(1,*)  Mwind           ; print*, 'Mwind       ', Mwind
      READ(1,*)  Mpat            ; print*, 'Mpat        ', Mpat    
      allocate(Nav(-Mpat:Mpat+11))
      allocate(arr(dir),arrx(dir,Nsite),arry(dir,Nsite),arrz(dir,Nsite))
      allocate(Wind(dir/2,-Mwind:Mwind))
      READ(1,*)  twin            ; print*, 'G-window    ', twin 
      READ(1,*)  tcc_win         ; print*, '<nn>-window ', tcc_win
      READ(1,*)  jc              ; print*, '<nn> freq.  ', jc
      READ(1,*)  p_ther  ; print*, 'therm. time in updates ',p_ther
      READ(1,*)  p_mes   ; print*, 'measure after # updates',p_mes
      READ(1,*)  p_prn   ; print*, 'print every # updates  ',p_prn
      READ(1,*)  p_wr    ; print*, 'write every # updates  ',p_wr   
      CLOSE(1)
             

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
        call SRANMAR(ij,kl)    
        call hamiltonian;  call initial
        print*, 'START MONTE CARLO'
                      
        ee: do;  step=step+1.d0; 
        i1=step/1.d9; ist=step-i1*1.d9; i1=ist
        do i=1,1000 
        call manager
        enddo 
        IF(step>p_ther) then; 
        p_ther=0         
        IF(mod(i1,p_mes)==0) call measure(i1); endif
        IF(mod(i1,p_prn)==0) call printing(i1)
        IF(mod(i1,p_wr)==0)  then; 
                             Print*, '!!!!! WRITE !!!!!'
                             PRINT*, '!!!!! DONE !!!!!'; endif
         enddo ee    

      
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

         tck=prtcls/beta
!         IF(CANONICAL .AND. ABS(tck-Nsite/4) > Mzmark) stop 'tck failed'

       END SUBROUTINE manager
    

!_______________________________________________________ 

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
            do da=1,rdir(site(name)); if (ms(da,name) == name1) then
            DIRECTION=da; exit; endif; enddo
 
           IF(DIRECTION<1 .OR. DIRECTION >rdir(site(name))) then
            print*, name, name1, DIRECTION; print*
            do da=1,rdir(site(name)); print*,da,ms(da,name); enddo
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
 
      DO j=1,rdir(site(p0)); j1=back(j,site(p0))  !************ 
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
	 
	 
      DO j=1,rdir(site(p0)); j1=back(j,site(p0)) !*********** 
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
 
      t3=tau(p0);      DO j=1,rdir(site(p0)); j1=back(j,site(p0))  !*********** 
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
 
                        DO j=1,rdir(site(p1)); j1=back(j,site(p1))  !********** 
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
	ta=time(n0);DO j=1,rdir(site(n0))
        j1=back(j,site(n0));nnn=ms(j,n0) 
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
         DO j=1,rdir(her); Iz(j)= Int_z(her,j); ENDDO 

!________________________________________________________________ 
!	x---------(p0)-----------x   ---> x---(p0)--o..(p1)..o 
!       t0	         	                        t1       t2 
!________________________________________________________________
 
        DO j=1,rdir(site(p0));j1=back(j,site(p0));nnn=ms(j,p0) !*************** 
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
      DO j=1,rdir(her); Iz(j)= Int_z(her,j); ENDDO 
      
        DO j=1,rdir(site(p1));j1=back(j,site(p1));nnn=ms(j,p1) !******** 
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
      DO j=1,rdir(her); Iz(j)= Int_z(her,j); ENDDO 
      
      DO j=1,rdir(site(n0));j1=back(j,site(n0));nnn=ms(j,n0) !************** 
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
        IF(CANONICAL .AND. ABS(tck-Nsite/4) > Mzmark) RETURN 
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
        IF(CANONICAL .AND. ABS(tck-Nsite/4) > Mzmark) RETURN 
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
!			       --> 	   |		          |     ^  
! x--(nl)------o..(n0)..x(nr)	   x--(nl)-x....(n0)......x(nr)	  V k1  | k
! tl          t=0 
! ______________________________________new kink___________________________
!        print *,'1em',site(imya)   
        n0=imya;  
        nl=ml(n0);tl=-tau(nl);
        k=oneoutof(rdir(site(n0))); nn=ms(k,n0); flnn=fl(nn); 
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
        CASE(1);  thop=s_h(fll,flnn)*Int_xy(1,here,k); 
        CASE(-1); thop=s_h(flnn,fll)*Int_xy(1,here,k); 
        END SELECT       

        SELECT CASE(dfl);
        CASE(1);
        if(k.eq.1.or.k.eq.2) then 
        all_phi=all_phi+Int_xy(2,here,k)
        else
        all_phi=all_phi-Int_xy(2,here,k)
        endif
        CASE(-1); 
        if(k.eq.1.or.k.eq.2) then
        all_phi=all_phi-Int_xy(2,here,k)
        else
        all_phi=all_phi+Int_xy(2,here,k)
        endif
        END SELECT
      
        p_a=-thop*rdir(site(n0))*tau0*w_hop(fln1,flnn)/w_hop(fll,fl0);
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
       
	kbb=back(k,site(n0)) 
!        print *,dfl
        IF(dfl.eq.1) then 
        arr(1)=arr(1)+arrx(k,site(n0))
        arr(2)=arr(2)+arry(k,site(n0))
        arr(3)=arr(3)+arrz(k,site(n0))
!        print *,arrx(k,site(n0)),arry(k,site(n0))
        ELSE IF(dfl.eq.-1) then 
        arr(1)=arr(1)-arrx(k,site(n0))
        arr(2)=arr(2)-arry(k,site(n0))      
        arr(3)=arr(3)-arrz(k,site(n0))  
        ENDIF
!        print *,'1',k,site(n0),site(nn),dfl,arr(1),arr(2)
!        print *,arrx(k,site(n0)),arrx(kbb,ms(k,n0))
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
!        print *,'em',site(imya) 
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
        CASE(1);  thop=s_h(fll,flnn)*Int_xy(1,here,k); 
        CASE(-1); thop=s_h(flnn,fll)*Int_xy(1,here,k) 
        END SELECT

        SELECT CASE(dfl);
        CASE(1);
        if(k.eq.1.or.k.eq.2) then
        all_phi=all_phi-Int_xy(2,here,k)
        else
        all_phi=all_phi+Int_xy(2,here,k)
        endif
        CASE(-1);
        if(k.eq.1.or.k.eq.2) then
        all_phi=all_phi+Int_xy(2,here,k)
        else
        all_phi=all_phi-Int_xy(2,here,k)
        endif
        END SELECT

 
        p_a=thop*rdir(site(n0))*tau0*w_hop(fln1,fln2)/w_hop(fll,fl0)
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
	kbb=back(k,site(n0))
          
        IF(dfl.eq.1) then
        arr(1)=arr(1)-arrx(k,site(n0))
        arr(2)=arr(2)-arry(k,site(n0))
        arr(3)=arr(3)-arrz(k,site(n0))
        ELSE IF(dfl.eq.-1) then
        arr(1)=arr(1)+arrx(k,site(n0))
        arr(2)=arr(2)+arry(k,site(n0))
        arr(3)=arr(3)+arrz(k,site(n0))
        ENDIF
!        print *,'2',k,site(nn),site(n0),dfl,arr(1),arr(2)
 
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
        k=oneoutof(rdir(site(n0)));
!        print *,n0,site(n0),rdir(site(n0)),nl,site(nl),rdir(site(nl))
        nn=ms(k,n0);flnn=fl(nn);IF(flnn==-3)RETURN
!        print *,k,nn,site(nn),rdir(site(nn))
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
        CASE(1); thop=s_h(fll,fln1)*Int_xy(1,here,k) ; 
        CASE(-1);thop=s_h(fln1,fll)*Int_xy(1,here,k) ; 
                        END SELECT
        
        SELECT CASE(dfl);
        CASE(1);
        if(k.eq.1.or.k.eq.2) then
        all_phi=all_phi+Int_xy(2,here,k)
        else
        all_phi=all_phi-Int_xy(2,here,k)
        endif
        CASE(-1);
        if(k.eq.1.or.k.eq.2) then
        all_phi=all_phi-Int_xy(2,here,k)
        else
        all_phi=all_phi+Int_xy(2,here,k)
        endif
        END SELECT
         
                 
        p_a=thop*rdir(site(n0))*tau0*w_hop(fln1,flnn)/w_hop(fll,fl0)
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
	kbb=back(k,site(n0))

        IF(dfl.eq.1) then
        arr(1)=arr(1)+arrx(k,site(n0))
        arr(2)=arr(2)+arry(k,site(n0))
        arr(3)=arr(3)+arrz(k,site(n0))
        ELSE IF(dfl.eq.-1) then
        arr(1)=arr(1)-arrx(k,site(n0))
        arr(2)=arr(2)-arry(k,site(n0))
        arr(3)=arr(3)-arrz(k,site(n0))
        ENDIF
!        print *,'3',k,site(n0),site(nn),dfl,arr(1),arr(2)
	 
  	Udiag=Udiag+pp; kinks=kinks+1 

        IF(ms(k,n0)/=n2  .OR. ms(kbb,n2)/=n0) then
        print*, 'insert right links violated!'
!        print *,site(n0),k,site(n2),kbb
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
        CASE(1); thop=s_h(fll,fln1)*Int_xy(1,here,k) ; 
        CASE(-1);thop=s_h(fln1,fll)*Int_xy(1,here,k) ; 
                       END SELECT
       
        SELECT CASE(dfl);
        CASE(1);
        if(k.eq.1.or.k.eq.2) then
        all_phi=all_phi-Int_xy(2,here,k)
        else
        all_phi=all_phi+Int_xy(2,here,k)
        endif
        CASE(-1);
        if(k.eq.1.or.k.eq.2) then
        all_phi=all_phi+Int_xy(2,here,k)
        else
        all_phi=all_phi-Int_xy(2,here,k)
        endif
        END SELECT

 
        p_a=thop*rdir(site(n0))*tau0*w_hop(fln1,fln2)/w_hop(fll,fl0)
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
	kbb=back(k,site(n0));  

        IF(dfl.eq.1) then
        arr(1)=arr(1)-arrx(k,site(n0))
        arr(2)=arr(2)-arry(k,site(n0))
        arr(3)=arr(3)-arrz(k,site(n0))
        ELSE IF(dfl.eq.-1) then
        arr(1)=arr(1)+arrx(k,site(n0))
        arr(2)=arr(2)+arry(k,site(n0))
        arr(3)=arr(3)+arrz(k,site(n0))
        ENDIF
!        print *,'4',k,site(nn),site(n0),dfl,arr(1),arr(2)

		 
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
        IF(CANONICAL .AND. ABS(tck-Nsite/4) > Mzmark) RETURN 
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
!______________________________________________________________________	    

!***********VERY IMPORTANT !!! - HAMILTONIAN IS CONSTRUCTED HERE********
!___________                      ______________________________________ 
            SUBROUTINE hamiltonian
            IMPLICIT NONE
            integer :: ir,kf,kb,kcf,kcb,irr,iz,iy,ix,j1,n12,n1
            integer :: idisp(6,27),ityp(Nsite/4) 
            integer :: ms_c(6,Nsite/4),noce(Nsite)
            integer :: ind,ixn,iyn,izn,ixm,ixp,iym,iyp,izm,izp,itp,ib
  
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


      IF(triangular==0.and.dim==3) then

      n12=Ns(1)*Ns(2) 
   
      do 519 iz=1,Ns(3)
      do 520 iy=1,Ns(2)
      do 521 ix=1,Ns(1)
        ind=(iz-1)*n12+(iy-1)*Ns(1)+ix
        izn=iz
        if(iz.gt.2 ) izn=2
        if(iz.eq.Ns(3)) izn=3
        iyn=iy
        if(iy.gt.2 ) iyn=2
        if(iy.eq.Ns(2)) iyn=3
        ixn=ix
        if(ix.gt.2 ) ixn=2
        if(ix.eq.Ns(1)) ixn=3
        ityp(ind)=(izn-1)*9+(iyn-1)*3+ixn
521   continue
520   continue
519   continue

        do 211 itp=1,27
        iz =(itp-1)/9+1
        iy =(itp-1-(iz-1)*9)/3+1
        ix =itp-(iz-1)*9-(iy-1)*3
        ixm=-1+((4 -ix)/3)*Ns(1)
        ixp= 1-(ix/3)*Ns(1)
        iym=-Ns(1)+((4 -iy)/3)*n12
        iyp= Ns(1)-(iy/3)*n12
        izm=-n12+((4 -iz)/3)*Nsite/4
        izp= n12-(iz/3)*Nsite/4
        idisp(1,itp)=ixm
        idisp(6,itp)=ixp
        idisp(2,itp)=iym
        idisp(5,itp)=iyp
        idisp(3,itp)=izm
        idisp(4,itp)=izp
211   continue


      do 152 ind=1,Ns(1)*Ns(2)*Ns(3)
      do 151 ib=1,6
      ms_c(ib,ind)=ind+idisp(ib,ityp(ind))
151   continue
152   continue

      do 153 ind=1,Nsite
      rdir(ind)=6
      Lat(ind)=mod(ind-1,4)+1
      noce(ind)=(ind-1)/4+1
153   continue


      Do n1=1,Nsite
       
      if(Lat(n1).eq.1) then

      ms(1,n1)=n1+1; back(1,n1)=1
      ms(2,n1)=n1+2; back(2,n1)=2
      ms(3,n1)=n1+3; back(3,n1)=1
      ms(4,n1)=4*(ms_c(5,noce(n1)))-2; back(4,n1)=4
      ms(5,n1)=4*(ms_c(4,noce(n1)))-1; back(5,n1)=4
      ms(6,n1)=4*(ms_c(6,noce(n1))); back(6,n1)=4

      arrx(1,n1)=0.5d0;  arry(1,n1)=1.d0; arrz(1,n1)=0.5d0
      arrx(2,n1)=0.5d0;  arry(2,n1)=0.5d0; arrz(2,n1)=1.d0
      arrx(3,n1)=1.d0;  arry(3,n1)=0.5d0; arrz(3,n1)=0.5d0
      arrx(4,n1)=-0.5d0; arry(4,n1)=-1.d0;arrz(4,n1)=-0.5d0
      arrx(5,n1)=-0.5d0;  arry(5,n1)=-0.5d0; arrz(5,n1)=-1.d0
      arrx(6,n1)=-1.d0; arry(6,n1)=-0.5d0; arrz(6,n1)=-0.5d0
 
      else if(Lat(n1).eq.2) then

      ms(1,n1)=n1-1; back(1,n1)=1
      ms(2,n1)=n1+1; back(2,n1)=1
      ms(3,n1)=n1+2; back(3,n1)=2
      ms(4,n1)=4*(ms_c(2,noce(n1)))-3; back(4,n1)=4
      ms(5,n1)=4*(ms_c(6,ms_c(2,noce(n1)))); back(5,n1)=5
      ms(6,n1)=4*(ms_c(2,ms_c(4,noce(n1))))-1; back(6,n1)=5

      arrx(1,n1)=-0.5d0;arry(1,n1)=-1.d0; arrz(1,n1)=-0.5d0
      arrx(2,n1)=0.d0;arry(2,n1)=-0.5d0;arrz(2,n1)=0.5d0
      arrx(3,n1)=0.5d0;arry(3,n1)=-0.5d0;arrz(3,n1)=0.d0
      arrx(4,n1)=0.5d0;arry(4,n1)=1.d0;  arrz(4,n1)=0.5d0
      arrx(5,n1)=-0.5d0; arry(5,n1)=0.5d0;arrz(5,n1)=0.d0
      arrx(6,n1)=0.d0;arry(6,n1)=0.5d0;arrz(6,n1)=-0.5d0


      else if(Lat(n1).eq.3) then

       ms(1,n1)=n1-1; back(1,n1)=2
       ms(2,n1)=n1-2; back(2,n1)=2
       ms(3,n1)=n1+1; back(3,n1)=3
       ms(4,n1)=4*ms_c(3,noce(n1))-3; back(4,n1)=5
       ms(5,n1)=4*ms_c(5,ms_c(3,noce(n1)))-2; back(5,n1)=6
       ms(6,n1)=4*ms_c(6,ms_c(3,noce(n1))); back(6,n1)=6

      arrx(1,n1)=0.d0;arry(1,n1)=0.5d0;arrz(1,n1)=-0.5d0
      arrx(2,n1)=-0.5d0;arry(2,n1)=-0.5d0; arrz(2,n1)=-1.0d0
      arrx(3,n1)=0.5d0;arry(3,n1)=0.0d0;arrz(3,n1)=-0.5d0
      arrx(4,n1)=0.5d0; arry(4,n1)=0.5d0;  arrz(4,n1)=1.0d0
      arrx(5,n1)=0.d0;arry(5,n1)=-0.5d0;arrz(5,n1)=0.5d0
      arrx(6,n1)=-0.5d0;arry(6,n1)=0.0d0;arrz(6,n1)=0.5d0

      else

      ms(1,n1)=n1-3; back(1,n1)=3
      ms(2,n1)=n1-2; back(2,n1)=3
      ms(3,n1)=n1-1; back(3,n1)=3
      ms(4,n1)=4*ms_c(1,noce(n1))-3; back(4,n1)=6
      ms(5,n1)=4*ms_c(5,ms_c(1,noce(n1)))-2; back(5,n1)=5
      ms(6,n1)=4*ms_c(4,ms_c(1,noce(n1)))-1; back(6,n1)=6

      arrx(1,n1)=-1.d0;  arry(1,n1)=-0.5d0;  arrz(1,n1)=-0.5d0
      arrx(2,n1)=-0.5d0;arry(2,n1)=0.5d0;arrz(2,n1)=0.0d0
      arrx(3,n1)=-0.5d0;arry(3,n1)=0.0d0; arrz(3,n1)=0.5d0
      arrx(4,n1)=1.d0;  arry(4,n1)=0.5d0;  arrz(4,n1)=0.5d0
      arrx(5,n1)=0.5d0;arry(5,n1)=-0.5d0;arrz(5,n1)=0.0d0
      arrx(6,n1)=0.5d0;arry(6,n1)=0.0d0; arrz(6,n1)=-0.5d0


      ENDIF

      ENDDO
             
      ENDIF


       do n1=1,Nsite
       do ib=1,6
!        print *,n1,ib
!       print *,'*************' 
!       print *,arrx(ib,n1)
!       print *,arry(ib,n1)
!       print *,arrz(ib,n1)
!       print *,arrx(back(ib,n1),ms(ib,n1))
!       print *,arry(back(ib,n1),ms(ib,n1))
!       print *,arrz(back(ib,n1),ms(ib,n1))

       if(abs(arrx(ib,n1)+arrx(back(ib,n1),ms(ib,n1))).gt.0.1) stop
       if(abs(arry(ib,n1)+arry(back(ib,n1),ms(ib,n1))).gt.0.1) stop
       if(abs(arrz(ib,n1)+arrz(back(ib,n1),ms(ib,n1))).gt.0.1) stop
       enddo
       enddo
      
!      stop 

        do ir=1,Nsite; do kf=1,rdir(ir)
          IF(kf.le.rdir(ir)/2) then
          Int_xy(1,ir,kf)=s_hop1
          ELSE
          Int_xy(1,ir,kf)=s_hop2
          ENDIF
          Int_xy(2,ir,kf)=0.0d0
        enddo; enddo
 
       do ir=1,Nsite; do kf=1,rdir(ir)
          IF(kf.le.rdir(ir)/2) then    
          Int_z(ir,kf)=Jz1
          ELSE                      
          Int_z(ir,kf)=Jz2
          endif
      enddo; enddo    
          
      Hsite=H_0          
!onsite disorder (if any)

      irr=Nsite/2.+1 
      do ir=1,Nsite
      IF(close_open==0) then
      Hsite(ir)=Hsite(ir)+W*(2.d0*rndm()-1.d0) ! random potential 
      else
      Hsite(ir)=Hsite(ir)+W*(ir-irr)**2        ! parabolic potential 
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
	nmnm=Nsite; kinks=0; arr=0
 
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
	step=0.d0; stat=1.d-15; Energy=0.d0; Energy2=0.d0;
        Upot=0.d0 ; Ekin=0.d0;sum_se=0.d0;stat_se=1.d-14 
        StatSqe=1.d-14; ddc1=0.d0; ddc_stat=1.d-14; UpDown=0.d0
        prntout=0; ddc2=0.d0; ddc3=0.d0
        rs_prn=0.d0; E_prn=0.d0; Egrand_prn=0.d0; Egrand2_prn=0.d0;
        N_prn=0.d0;N0_prn=0.d0;N1_prn=0.d0;N2_prn=0.d0;N3_prn=0.d0; 
        UD_prn=0.d0; histUD=0.d0

        all_phi=0.d0

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
        END SUBROUTINE initial
        
!________________________________________________________________
! VERY IMPORTANT PROCEDURE OF PERFORMING MEASUREMENTS 
!________________________________________________________________
       SUBROUTINE measure(i1)
       INTEGER, INTENT(IN) :: i1 
       integer ::  kp,kw,k, i,name,name1, nl,jgr,kat,here, katN
       integer :: fi,fm,ix,iy,iz,nn,ijk,nnn,nij,ind,k1,k2,k3,k4,k5
       integer :: id_x,id_y,id_z,istart,isx,isy,isz,j2,j3,j4,j5
       integer, parameter ::    kat0=0
       double precision :: kpa,tleft,tright,addt,ht,ddt,pp
       double precision :: b1,b2,tnorm,cop,nm(0:fl_max),sub_n_et(4)
       double precision :: se_1,se_2,se_3,se_4,tempe1,tempe2,tempe0
       double precision :: ndistr(Nsite),tempdd0,tempdd1,tempdd2
       double precision :: Mup, Mdown,sub_n(Nsite),tempdd
       double precision :: s_1,s_2,s_3,s_4,s_5,s_345_c,s_345_s,rd1,s_12       
       integer :: is_1,is_2,is_3,is_4,ind_1,ind_2,ind_3,ind_4,mmm
       double precision :: sum_is,sum_ind,aa1,aa2,aa3,bb
       katN=2*Nmax       
       IF(.not.present) then           ! NO-WORM STATISTICS
       
       kpa=prtcls/beta; kp=NINT(kpa)   ! PARTICLE NUMBERS
!       print *,kpa,Nsite/4
       IF(CANONICAL .AND. kpa.ne.Nsite/4) RETURN    ! FOR CANONICAL calc
            
       stat=stat+1.d0                  ! PARTITION FUNCTION    
       IF(kp>Mpat.or.kp<-Mpat)         stop 'kp > Mpat; increase Npat' 
       Nav(kp)=Nav(kp)+1               ! PARTICLE NUMBERS HISTOGRAM   
       Nav(Mpat+1)=Nav(Mpat+1)+kp;     ! PARTICLE NUMBER AVERAGE 
       
       sub_n=0.0d0
       do i=1,nmnm
       nn=nlist(i)
       ijk=Lat(site(nn))
       sub_n(ijk)=sub_n(ijk)+Pat(fl(nn))*tau(nn)/beta
       enddo

       s_1=sub_n(1)/(1.d0*Nsite/4.d0)
       s_2=sub_n(2)/(1.d0*Nsite/4.d0)
       s_3=sub_n(3)/(1.d0*Nsite/4.d0)
       s_4=sub_n(4)/(1.d0*Nsite/4.d0)
     
       Nav(Mpat+2)=Nav(Mpat+2)+s_1
       Nav(Mpat+3)=Nav(Mpat+3)+s_2
       Nav(Mpat+4)=Nav(Mpat+4)+s_3
       Nav(Mpat+5)=Nav(Mpat+5)+s_4
       Nav(Mpat+6)=Nav(Mpat+6)+abs(s_1-s_2)

       tempdd1=(s_1-s_2)**2+(s_1-s_3)**2+(s_1-s_4)**2
       tempdd2=(s_2-s_3)**2+(s_2-s_4)**2+(s_3-s_4)**2
       tempdd0=sqrt(tempdd1+tempdd2)
       Nav(Mpat+7)=Nav(Mpat+7)+tempdd0



       Ekin=Ekin-kinks               ! kinetic energy 
       Upot=Upot+Udiag               ! Potential energy 
       Energy=Energy+Udiag-kinks     ! Total energy           
       Energy2=Energy2+(Udiag-kinks)**2           

       DO j=1,dim;   
       if(arr(j).gt.0.0d0) then
       kw=(arr(j)+0.1d0)/Ns(j)
       else 
       kw=(arr(j)-0.1d0)/Ns(j)
       endif
       IF(ABS(kw)>Mwind) stop 'kw > Mwind; increase Mwind'    
       Wind(j,kw)=Wind(j,kw)+1; 
       END DO


       IF(mod(i1,200)==0) then

       call Equal_time_nn(3,4,1,4,2,4)
         
       do nnn=1,Nsite,4
       do mmm=1,Nsite,4  
            
       istart=(mmm-1)/4+1   
       isz=(istart-1)/(Ns(1)*Ns(2))+1
       isy=(istart-(isz-1)*Ns(1)*Ns(2)-1)/Ns(1)+1
       isx= istart-(isz-1)*Ns(1)*Ns(2)-(isy-1)*Ns(1)

       is_1=mmm
       is_2=mmm+1
       is_3=mmm+2
       is_4=mmm+3  

       ind=(nnn-1)/4+1
       iz=(ind-1)/(Ns(1)*Ns(2))+1
       iy=(ind-(iz-1)*Ns(1)*Ns(2)-1)/Ns(1)+1
       ix= ind-(iz-1)*Ns(1)*Ns(2)-(iy-1)*Ns(1)
        
       ind_1=nnn
       ind_2=nnn+1
       ind_3=nnn+2
       ind_4=nnn+4
       
       id_z=abs(isz-iz)
       if(id_z.gt.Ns(3)/2) id_z=Ns(3)-id_z
       id_y=abs(isy-iy) 
       if(id_y.gt.Ns(2)/2) id_y=Ns(2)-id_y
       id_x=abs(isx-ix)
       if(id_x.gt.Ns(1)/2) id_x=Ns(1)-id_x

       j3=id_x;j4=id_y;j5=id_z

       sum_is= Sqe2(is_1)+ Sqe2(is_2)- Sqe2(is_3)-Sqe2(is_4)
       sum_ind=Sqe2(ind_1)+Sqe2(ind_2)-Sqe2(ind_3)-Sqe2(ind_4)
       ddc1(j3,j4,j5)=ddc1(j3,j4,j5)+sum_is*sum_ind

       sum_is= Sqe2(is_1)- Sqe2(is_2) +Sqe2(is_3)-Sqe2(is_4)
       sum_ind=Sqe2(ind_1)-Sqe2(ind_2)+Sqe2(ind_3)-Sqe2(ind_4)
       ddc2(j3,j4,j5)=ddc2(j3,j4,j5)+sum_is*sum_ind

       sum_is= Sqe2(is_1)-Sqe2(is_2)- Sqe2(is_3)+Sqe2(is_4)
       sum_ind=Sqe2(ind_1)-Sqe2(ind_2)-Sqe2(ind_3)+Sqe2(ind_4)
       ddc3(j3,j4,j5)=ddc3(j3,j4,j5)+sum_is*sum_ind


       ddc_stat(j3,j4,j5)=ddc_stat(j3,j4,j5)+1.0d0

       enddo 
       enddo 
       endif

       if(mod(i1,200)==0) then
       open(199,file='ddcr_final_z.dat',status='replace')
       open(299,file='ddcr_final_y.dat',status='replace')
       open(399,file='ddcr_final_x.dat',status='replace')
       do k3=0,Ns(1)/2
       do k4=0,Ns(2)/2
       do k5=0,Ns(3)/2
       aa1=ddc1(k3,k4,k5)
       aa2=ddc2(k3,k4,k5) 
       aa3=ddc3(k3,k4,k5)
       bb=ddc_stat(k3,k4,k5)        
       write(199,*) k3,k4,k5,aa1/bb
       write(299,*) k3,k4,k5,aa2/bb
       write(399,*) k3,k4,k5,aa3/bb
       enddo 
       enddo
       enddo
       close(399)
       close(299)
       close(199)


       ENDIF      

       do i=1,Nsite            
       sub_n_et(Lat(i))=sub_n_et(Lat(i))+Sqe2(i)
       enddo

       se_1=sub_n_et(1)/(1.d0*Nsite/4.d0)
       se_2=sub_n_et(2)/(1.d0*Nsite/4.d0)
       se_3=sub_n_et(3)/(1.d0*Nsite/4.d0)
       se_4=sub_n_et(4)/(1.d0*Nsite/4.d0)

       tempe1=(se_1-se_2)**2+(se_1-se_3)**2+(se_1-se_4)**2
       tempe2=(se_2-se_3)**2+(se_2-se_4)**2+(se_3-se_4)**2
       tempe0=sqrt(tempe1+tempe2)
       
       sum_se=sum_se+tempe0
       stat_se=stat_se+1.0d0 
       endif   
       END SUBROUTINE measure
!***************************************
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
     
      subroutine Equal_time_nn(ixa,ixb,iya,iyb,iza,izb)
      INTEGER :: i,k,loop,name,n2,n11(Nsite),cou,n,nnn,lll,mmm
      DOUBLE PRECISION :: t2,b1,b2,b3,ndistr(Nsite),n0
      DOUBLE PRECISION :: t00,t0,tslice
      INTEGER :: ixa,ixb,iya,iyb,iza,izb,ivist(Nsite)
      
      ivist=0; Sqe=0.d0; Sqe2=0.d0

      cou=1
666   i=oneoutof(nmnm); name=nlist(i);
      if(Lat(site(name)).ne.1) goto 666
       
      t00=time(name); t0=0.d0; tslice=rndm()*beta;
      start=site(name); ivist(start)=ivist(start)+1
                                                                    
      e: DO; IF(t0+tau(name) .GT.  tslice ) EXIT e;
             t0=t0+tau(name); name=mr(name); enddo e 
      k=site(name); ndistr(k)=fl(name); n11(cou)=name  ! first site
      tslice=s_t(plus,t00,tslice)
      n0=ndistr(k)
      Sqe(start)=n0*n0
      Sqe2(start)=n0

!     print *,k,ndistr(k)
      
      n2=name
      do nnn=1,Ns(1)-1
      n2=ms(ixa,n2);
      n2=ms(ixb,n2); t2=-s_t(minus,tslice,time(n2))
                   ba : DO; IF(t2+tau(n2)>0) exit ba;
                   t2=t2+tau(n2); n2=mr(n2);  ENDDO ba
      k=site(n2); ndistr(k)=fl(n2); cou=cou+1; n11(cou)=n2
      ivist(site(n2))=ivist(site(n2))+1
      Sqe(site(n2))=ndistr(k)*n0
      Sqe2(site(n2))=ndistr(k)
      enddo

      do mmm=1,Ns(1) 
      n2=n11(mmm)
      do nnn=1,Ns(2)-1
      n2=ms(iya,n2);
      n2=ms(iyb,n2); t2=-s_t(minus,tslice,time(n2))
                   ca : DO; IF(t2+tau(n2)>0) exit ca;
                   t2=t2+tau(n2); n2=mr(n2);  ENDDO ca
      k=site(n2); ndistr(k)=fl(n2); cou=cou+1; n11(cou)=n2
      ivist(site(n2))=ivist(site(n2))+1
      Sqe(site(n2))=ndistr(k)*n0
      Sqe2(site(n2))=ndistr(k)
      enddo
      enddo


      do mmm=1,Ns(1)*Ns(2)
      n2=n11(mmm)
      do lll=1,Ns(3)-1
      n2=ms(iza,n2);
      n2=ms(izb,n2); t2=-s_t(minus,tslice,time(n2))
                   da : DO; IF(t2+tau(n2)>0) exit da;
                   t2=t2+tau(n2); n2=mr(n2);  ENDDO da
      k=site(n2); ndistr(k)=fl(n2); cou=cou+1; n11(cou)=n2
      ivist(site(n2))=ivist(site(n2))+1
      Sqe(site(n2))=ndistr(k)*n0
      Sqe2(site(n2))=ndistr(k)
      enddo
      enddo
      
      do mmm=1,Ns(1)*Ns(2)*Ns(3)
      n2=n11(mmm)
      n2=ms(1,n2); t2=-s_t(minus,tslice,time(n2))
                   ea : DO; IF(t2+tau(n2)>0) exit ea;
                   t2=t2+tau(n2); n2=mr(n2);  ENDDO ea
      k=site(n2); ndistr(k)=fl(n2)
      ivist(site(n2))=ivist(site(n2))+1
      Sqe(site(n2))=ndistr(k)*n0
      Sqe2(site(n2))=ndistr(k)

      n2=n11(mmm)
      n2=ms(2,n2); t2=-s_t(minus,tslice,time(n2))
                   fa : DO; IF(t2+tau(n2)>0) exit fa;
                   t2=t2+tau(n2); n2=mr(n2);  ENDDO fa
      k=site(n2); ndistr(k)=fl(n2)
      ivist(site(n2))=ivist(site(n2))+1
      Sqe(site(n2))=ndistr(k)*n0
      Sqe2(site(n2))=ndistr(k)

      n2=n11(mmm)
      n2=ms(3,n2); t2=-s_t(minus,tslice,time(n2))
                   ga : DO; IF(t2+tau(n2)>0) exit ga;
                   t2=t2+tau(n2); n2=mr(n2);  ENDDO ga
      k=site(n2); ndistr(k)=fl(n2)
      ivist(site(n2))=ivist(site(n2))+1
      Sqe(site(n2))=ndistr(k)*n0
      Sqe2(site(n2))=ndistr(k)
  
      
      enddo
      
      do i=1,Nsite
      if(ivist(i).eq.0) then 
      print *,i,'is not visited'
      stop
      else if (ivist(i).gt.1) then
      print *,i,'is visited more than one time'
      stop
      endif
      enddo
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
         DOUBLE PRECISION :: xx,yy,qq,xxn,yyn,Ecan,Egrand,temp1,temp2
         DOUBLE PRECISION ::lambda,kappa,deltaN,zu,Egrand2
!         REAL*4 :: result(11)
	 INTEGER, INTENT(IN) :: i1
	    
       print*; print* 
       print*, 'MC steps (in millions)=', i1/1000000
       print*, 'stat=', stat
       print*, '<n^2> stat       =', AINT(stat2)
       print*, 'StatSqe         =', StatSqe
       print*, '  nmnm=', nmnm, 'present=', present
       print*, 'particles now =',prtcls/beta
       
       IF(p_ther==0 . AND. stat.ne.0.d0) then
       prntout=prntout+1
       if(prntout.gt.prnt_max) stop
       result(1)=mag_field-onsite
       xxn=Nav(Mpat+1)/stat
       yyn=xxn+Nsite*fl_max*0.5d0 
       print*
       N_prn(prntout)=xxn
       call ERSTAT(N_prn,prntout,amax,tmax,amin,tmin)
       print 705, N_prn(prntout), (amax-amin)/2.
 705   format(6x,'<N>         =',g12.5,4x,'+-',g10.3)
       result(8)= N_prn(prntout); result(9)=(amax-amin)/2.
       print*
       N1_prn(prntout)=Nav(Mpat+2)/stat
       call ERSTAT(N1_prn,prntout,amax,tmax,amin,tmin)
       print 799, N1_prn(prntout), (amax-amin)/2.
 799   format(6x,'<N_1>         =',g12.5,4x,'+-',g10.3)
       result(12)= N1_prn(prntout); result(13)=(amax-amin)/2.
       print*
       N2_prn(prntout)=Nav(Mpat+3)/stat
       call ERSTAT(N2_prn,prntout,amax,tmax,amin,tmin)
       print 899, N2_prn(prntout), (amax-amin)/2.
 899   format(6x,'<N_2>         =',g12.5,4x,'+-',g10.3)
       result(14)= N2_prn(prntout); result(15)=(amax-amin)/2.
       print*
       N3_prn(prntout)=Nav(Mpat+4)/stat
       call ERSTAT(N3_prn,prntout,amax,tmax,amin,tmin)
       print 499, N3_prn(prntout), (amax-amin)/2.
499    format(6x,'<N_3>         =',g12.5,4x,'+-',g10.3)
       result(16)= N3_prn(prntout); result(17)=(amax-amin)/2.
       print*
       N4_prn(prntout)=Nav(Mpat+5)/stat
       call ERSTAT(N4_prn,prntout,amax,tmax,amin,tmin)
       print 398, N4_prn(prntout), (amax-amin)/2.
398    format(6x,'<N_4>         =',g12.5,4x,'+-',g10.3)
       result(18)= N4_prn(prntout); result(19)=(amax-amin)/2.
       print*
       N5_prn(prntout)=Nav(Mpat+6)/stat
       call ERSTAT(N5_prn,prntout,amax,tmax,amin,tmin)
       print 678, N5_prn(prntout), (amax-amin)/2.
678    format(6x,'<N_5>         =',g12.5,4x,'+-',g10.3)
       result(20)= N5_prn(prntout); result(21)=(amax-amin)/2.
       print*
       N6_prn(prntout)=Nav(Mpat+7)/stat
       call ERSTAT(N6_prn,prntout,amax,tmax,amin,tmin)
       print 638, N6_prn(prntout), (amax-amin)/2.
638    format(6x,'<N_1_2>         =',g12.5,4x,'+-',g10.3)
       result(22)= N6_prn(prntout); result(23)=(amax-amin)/2.
       print*
       N7_prn(prntout)=Nav(Mpat+8)/stat
       call ERSTAT(N7_prn,prntout,amax,tmax,amin,tmin)
       result(24)= N7_prn(prntout); result(25)=(amax-amin)/2.
       print*
       N8_prn(prntout)=sum_se/stat_se
       call ERSTAT(N8_prn,prntout,amax,tmax,amin,tmin)
       result(26)= N8_prn(prntout); result(27)=(amax-amin)/2.
       print*
       N9_prn(prntout)=Nav(Mpat+10)/stat
       call ERSTAT(N9_prn,prntout,amax,tmax,amin,tmin)
       result(28)= N9_prn(prntout); result(29)=(amax-amin)/2.
       print*


       Ecan=Energy/(stat*inv_temp)+(mag_field-onsite)*xxn
       E_prn(prntout)=Ecan/Nsite
       call ERSTAT(E_prn,prntout,amax,tmax,amin,tmin)
       result(2)= E_prn(prntout); result(3)=(amax-amin)/2.
       print*

       Egrand=Energy/stat
       Egrand_prn(prntout)=Egrand/Nsite/inv_temp
       call ERSTAT(Egrand_prn,prntout,amax,tmax,amin,tmin)
       result(30)= Egrand_prn(prntout); result(31)=(amax-amin)/2.
       print*

       Egrand2=Energy2/stat
       Egrand2_prn(prntout)=Egrand2/Nsite/Nsite/inv_temp/inv_temp
       call ERSTAT(Egrand2_prn,prntout,amax,tmax,amin,tmin)
       result(32)= Egrand2_prn(prntout); result(33)=(amax-amin)/2.
       print*

       IF(dim==1) then
       print *,'wrong dimension'
       stop
       ELSE IF(dim==2) then
       print *,'wrong dimension'
       stop
       ELSE
       
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
       result(6)= rs_prn(prntout); result(7)=(amax-amin)/2.
       ENDIF

       ENDIF

       open(1,file='result.dat')
       write(1,810) Ns(2),1/inv_temp,0.5d0*xy_coupling1,result
       close(1)       
 810   format(i4,3f10.5,16(3x,f30.20,3x,f30.20))

       open(1,file='energy.dat')
       write(1,811) Ns(2),1/inv_temp,0.5d0*xy_coupling1
       write(1, 812) result(30),result(31)
       close(1)       
 811   format(i4,3f10.5)
 812   format(f30.20, 3x, f30.20)

       temp1=result(32)-result(30)**2
       temp2=result(33)+abs(2.0d0*result(30))*result(31)
       open(1,file='result2.dat')
       write(1,815) Ns(2),1/inv_temp,0.5d0*xy_coupling1,temp1,temp2
       close(1)
 815   format(i4,2f10.5,1(3x,f30.20,3x,f30.20))

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
