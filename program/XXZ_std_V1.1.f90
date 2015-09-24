  !******************************************************************
  !This code is for XXZ model with hamiltonian:
  !H=sum (-J*S^z_i*S^z_j-Q*(S^x_i*S^x_j+S^y_i*S^y_j))
  !(1) J>0,Q>0: Ferromagnetic
  !(2) J<0,Q<0: Antiferromagnetic
  
  !*******************************************************************
  ! Defined on the square lattice with built-in error bars
  ! using the blocking technique. 
  ! Later generalized on the cubic lattice and pyroclore lattice by Yuan in 2015,Aug.

  ! Error bars are calculated using the blocking technique. 
  ! Let ' \chi^2 = \sum_{t=1}^{T} <(O_t -<O>)^2> ' be the squared chi for
  ! 'T' blocks of observable 'O'. Assuming each block of data is independent
  ! of each other, the error bar of 'O' is given by 'Dev = \sqrt{\chi^2/T(T-1)}.

  ! Reliabity of the obtained errors is monitored by t=1 correlation,
  ! for which tolerance is set by variable 'tol' (default: tol=0.15d0).

  ! To calculate the error bar of composite quantity like Binder ratios 'Q' and 
  ! specific heat 'C', an auxilliary variable is defined as derivative of
  ! 'Q' or 'C'. For instance,
  ! (1). 'Q = <O_2>/<O_1>, the auxillary variable is 
  !      AQ = (1/<O_1>^2) (<O_1> Q_2 - <Q_2> Q_1), which is again a time series
  ! (2). 'C = V(<O^2>-<O>^2),
  !      AC = V(O^2 - 2<O> O).
  ! The error bar of 'Q' ('C') is taken that of 'AQ' ('AC').

  ! Results are written into a special file 'dat.***' if the number of
  ! blocks is less than 125 or correlation is too big. Data in each 
  ! block will be also printed out in this case.

  ! Default number of extensive simulation is 'NBlck=1024'.

  ! For test purpose, for which huge amount of information will be 
  ! printed out, 'NBlck' should be set smaller but >2.

  ! Dynamical behavior is not studied.

  ! 'my_vrbls.f90', 'carlo.f90', 'monte.f90', and 'measure.f90'
  ! need to be modified for new projects.

  !  Look for 'PROJECT-DEPENDENT'.

  !  Author: Youjin Deng
  !  Date  : Jan. 25th, 2011.
  !*******************************************************************


  !============== Variable Module ====================================
  MODULE my_vrbls
    IMPLICIT NONE

    !-- common parameters and variables ------------------------------
    ! THIS IS ALMOST PROJECT-INDEPENDENT 
    double precision, parameter :: tm32   = 1.d0/(2.d0**32.d0)
    double precision, parameter :: eps    = 1.d-14            ! very small number
    double precision, parameter :: Pi     = 3.141592653       ! pi
    double precision, parameter :: tol    = 0.15d0            ! tolerance for Cor
    logical                     :: prt                        ! flag for write2file
    integer,          parameter :: Mxint  = 2147483647        ! maximum integer
    integer,          parameter :: Mnint  =-2147483647        ! minimum integer

    integer,          parameter :: MxBlck =8192              ! maximum number of blocks for statistics
    integer,          parameter :: MnBlck = 1                 ! minimum number of blocks

    integer            :: NBlck                           ! # blocks
    integer            :: Nswee                           ! # MC sweeps between measurements
    integer            :: Ntoss                           ! # MC sweeps for equilibration
    integer            :: Nsamp                           ! # samples in unit 'NBlck'
    integer            :: NmeasCorr    ! # normal measure times between correlation measurement
    integer            :: NSave
    integer            :: TotSamp                             
    integer            :: IsLoad

    integer            :: Vol                                 ! actual volumn
    double precision   :: Norm                                ! normalization factor (1/Vol)

    double precision   :: rmuv,radv                           ! auxillary variables to draw a random vertex
    double precision   :: rmue,rade                           ! auxillary variables to draw a random neighboring edge
    !-----------------------------------------------------------------

    !-- parameters and variables -------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    character(10)       :: LatticeName
    integer             :: Dim                                 ! dimensionality
    integer             :: nnb                                 ! # neighboring edges
    integer,allocatable :: Back(:,:)                           ! nnb direction revert
    integer             :: NSub                                ! # sublattices
    integer,allocatable :: Sub(:)                              ! Sublattice index for each vertex
    DOUBLE PRECISION,allocatable  :: RealVector(:, :)
    DOUBLE PRECISION,allocatable  :: LatticeVector(:, :)
    DOUBLE PRECISION,allocatable  :: SubVector(:, :)
    integer, parameter  :: MxL = 8                            ! maximum lattice size 
    integer,parameter   :: MaxKink=1024
    integer             :: MxV                                 ! maximum volumn
    integer,parameter   :: MxOmega=128                         ! maximum Matsubara frequency
    integer, allocatable:: L(:)                                ! actual linear lattice size 

    character(8 ), parameter :: ident = 'hs_sqa 0'             ! identifier
    character(12), parameter :: file1 = 'hs_sqa0.dat'          ! datafile
    character(8), parameter ::  ident2 = 'sta_corr'            ! identifier
    character(100)  :: file2        ! static correlation
    character(100)  :: file3        ! frequency correlator
    character(100)  :: file4        ! middle blocks
    character(100) :: file_config
    character(100) :: file_config2
    !-----------------------------------------------------------------

    !-- Lattice and state --------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    integer, allocatable :: Ngs(:,:)     ! neighbouring-vertex table
    integer(1), allocatable  :: dr(:,:,:)  ! auxillary variables to measure wrapping probability
    !-----------------------------------------------------------------

    !-- Observables --------------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    integer, parameter :: NPhase = 100
    integer, parameter :: NObs_b = 18                         ! #basic observables
    integer, parameter :: NObs_c = 3                          ! #composite observables
    integer, parameter :: NObs   = NObs_b+NObs_c              ! Total # observables
    !-----------------------------------------------------------------

    !-- Statistics ---------------------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    double precision, allocatable :: Quan(:)                  ! Measured quantities
    double precision, allocatable :: Obs(:,:)                 ! 1st--#quan.  2nd--#block
    double precision, allocatable :: Ave(:), Dev(:), Cor(:)   ! average, error bars, and correlation of observables
    double precision, allocatable :: StaticCorr(:)            ! static spin-spin correlations
    double precision, allocatable :: Momentum(:,:)                   ! k list
    integer :: Nk
    double precision, allocatable :: DynamicalCorr(:,:)      ! spin-spin correlations in k-omegarepresentation
    double precision, allocatable :: ReExp(:), ImExp(:)
    !-----------------------------------------------------------------

    !-- Random-number generator---------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    integer, parameter           :: mult=32781
    integer, parameter           :: mod2=2796203, mul2=125
    integer, parameter           :: len1=9689,    ifd1=471
    integer, parameter           :: len2=127,     ifd2=30
    integer, dimension(1:len1)   :: inxt1
    integer, dimension(1:len2)   :: inxt2
    integer, dimension(1:len1)   :: ir1
    integer, dimension(1:len2)   :: ir2
    integer                      :: ipnt1, ipnf1
    integer                      :: ipnt2, ipnf2
    integer, parameter           :: mxrn = 10000
    integer, dimension(1:mxrn)   :: irn(mxrn)

    integer                      :: Seed                      ! random-number seed
    integer                      :: nrannr                    ! random-number counter
    !-----------------------------------------------------------------

    !-- time-checking variables --------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    character( 8)         :: date
    character(10)         :: time
    character(5 )         :: zone
    integer, dimension(8) :: tval
    double precision      :: t_prev, t_curr, t_elap
    integer               :: h_prev, h_curr 
    double precision      :: t_simu, t_mcmc, t_meas
    !-----------------------------------------------------------------
    
    !-- data structs -------------------------------------------------

    double precision       :: Comp                             ! common variable in probability
    double precision       :: JJ                               ! common variable 
    double precision       :: Q
    double precision       :: QQ
    double precision       :: ExternalField
    double precision       :: beta                             ! common variable 
    double precision       :: Number
    double precision       :: Const
                                               
    double precision,allocatable :: KinkTime(:,:)      !Time point of every kink
    integer,allocatable          :: SegmentNum(:)             !Track the number of segments
    integer ,allocatable         :: SegmentState(:,:)  !Track the state of segments    
    integer ,allocatable         :: TailPoint(:)
    
    integer,allocatable          :: PrevSite(:,:)
    integer,allocatable          :: NextSite(:,:)
    integer,allocatable          :: NeighSite(:,:)        
    integer,allocatable          :: NeighVertexNum(:,:)
    integer,allocatable          :: NearestSite(:,:,:) 
    
    double precision      :: IraTime
    integer               :: IraVertex
    double precision      :: MashaTime
    integer               :: MashaVertex
    integer :: IraSite,MashaSite
    
    !--- auxillary data for measuring ----------------------------
    double precision     :: W0
    double precision     :: W1
    double precision     :: TotalW
    double precision     :: EnergyCheck
    integer              :: dWt, Wt, WindT
    integer    :: Direction
    integer,allocatable    :: dWR(:)
    integer,allocatable    :: WR(:)
    integer,allocatable    :: WindR(:)
    integer,allocatable    :: KinkNum(:)
    character(13)       :: RunName(13)
    double precision    :: RunNum(10)
    double precision    :: RunNum0(10)
    !double precision,parameter  :: P1=1.0/19
    !double precision,parameter :: P2=P1+5.0/19
    !double precision,parameter  :: P3=P2+5.0/19
    !double precision,parameter  :: P4=P3+1.0/19
    !double precision,parameter  :: P5=P4+1.0/19
    !double precision,parameter  :: P6=P5+5.0/19
    double precision,parameter  :: P1=1.0/4
    double precision,parameter  :: P2=P1+1.0/4
    double precision,parameter  :: P3=P2+1.0/4
  END MODULE my_vrbls

  !=====Main routine for bond percolation on square lattice ==========
  PROGRAM main
    use my_vrbls
    implicit none
    integer :: itoss,isamp,iblck,i,k,Site,Low
    character(99) :: filename

    print *, 'Dim'
    read *, Dim

    allocate(L(1:Dim))

    allocate(dWR(1:Dim))
    allocate(WR(1:Dim))
    allocate(WindR(1:Dim))
    NSub = 1
    nnb = 2*Dim

    print *, 'Lattice, L, beta, J,Q,Hz,Ntoss, Nsamp, Nswee, NSave, Seed, NBlck, IsLoad, Outputfile, file2, file3, file4'
    read  *,  LatticeName,L(:),beta,JJ,Q,ExternalField,Ntoss,Nsamp,Nswee,NSave,Seed,NBlck,IsLoad,filename,file2,file3,file4

    !!!!!!Treatments for lattices with multiple sublattices!!!!!!!!!!
    if(LatticeName=='Pyrochlore') then
	NSub = 4
    endif
    if(LatticeName=='Pyrochlore' .and. Dim/=3) then
	print *,"Pyrochlore on ", Dim, "dimension does not exist!"
	stop
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    MxV = 1 
    do i = 1, Dim
	MxV = MxV*MxL
    enddo
    MxV = MxV*NSub

    Vol=1
    do i = 1, Dim
	Vol = Vol*L(i)
    enddo
    Vol = Vol*NSub

    NmeasCorr = 50

    allocate(dr(nnb, 1:NSub, 1:Dim))
    allocate(back(nnb, 1:NSub))
    allocate(Ngs(nnb, 1:MxV))
    allocate(Sub(1:Vol))
    allocate(LatticeVector(1:Dim, 1:Dim))
    allocate(SubVector(1:NSub, 1:Dim))
    allocate(RealVector(1:Vol, 1:Dim))

    write(file_config,"(a1,a99)") "c",filename
    write(file_config2,"(a1,a99)") "w",filename

    !!!!!Jzz
    JJ=JJ/4.0

    !!!!!J_porp= 0.5*Jxx(Jyy)
    Q=Q/8.0
    QQ=abs(Q)

    !!!!!H
    ExternalField = ExternalField/2.0

    allocate(KinkNum(MxV))
    allocate(KinkTime(MxV,MaxKink))    
    allocate(SegmentNum(MxV))          
    allocate(SegmentState(MxV,MaxKink)) 
    allocate(TailPoint(MxV))
    
    allocate(PrevSite(MxV,MaxKink))
    allocate(NextSite(MxV,MaxKink))
    allocate(NeighSite(MxV,MaxKink))        
    allocate(NeighVertexNum(MxV,MaxKink))
    allocate(NearestSite(MxV,MaxKink,nnb)) 

    TotSamp = Nsamp*NBlck/1024
    RunNum(:)=0
    KinkNum(:)=0
    W0=0.0;W1=0.0
    RunName(1)="Delete worm:"
    RunName(2)="Move   worm:"
    RunName(3)="Create kink:"
    RunName(4)="Delete kink:"


    !--- Initialization ----------------------------------------------
    call t_elapse(0)         ! '0' for the 1st time.
    call carlo
    call def_lattice
    call t_elapse(1)         ! '1' for set up

    call winding_number()
    WR(:)=WindR(:)
    Wt=WindT

    !--- Equilibration -----------------------------------------------
    if(IsLoad/=1) then
        do iblck = 1, NBlck
	    do itoss = 1, Ntoss
	      call monte
	    enddo
        enddo
    endif
    call t_elapse(2)         ! '2' for equilibration
    print *, "thermalization done!"


    !--- Simulation --------------------------------------------------
    if(IsLoad==1) then
        open(18,file=file_config2,access="append")
    else
        open(18,file=file_config2)
    endif
        
    print *, "start simulation..."
    do iblck = 1, NBlck
      do isamp = 1, Nsamp
        call monte
        call measure
	if(mod(isamp,NmeasCorr)==0)  call measure_Corr
        call coll_data(iblck)
      enddo
      print *, "simulation: Block", iblck, " done!"
      call norm_Nsamp(iblck)
      call t_elapse(-1)      ! '-1' just for trace the time

      if(mod(iblck, NSave)==0) then
	  call write2file_corr(iblck)
          call midwrite2file(iblck)
          call saveconfig
	  print*, iblck,"save data and configuration"
      endif
    enddo
    close(18)

    call t_elapse(3)         ! '3' for markov-chain, including simulation and measurement

    !--- Statistics --------------------------------------------------
    call stat_alan
    call write2file
    call write2file_corr(NBlck)
    call saveconfig
    call t_elapse(4)         ! '4' total time 


  CONTAINS

  !*******************************************************************
  !            Beginning of PROJECT-DEPENDENT part 
  !*******************************************************************
  !==============Initialization ======================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE carlo
    implicit none
    integer  :: i,j, k, x, y, xp, xm, ym, yp
    integer  :: kb,temp
    double precision :: phase
    
    Comp=1.0/beta

    Number=0
    !Comp=1/(beta/<n>)^2/<n>=<n>/beta^2 
    !<n> is the average kink number per site
    ! <n>=beta*2Vol*<SxSx+SySy>/Vol~beta
    !-- define time-line list ----------------------------------------
    if(IsLoad==1) then
        open(9,file=file_config)
        60 format(i6)
        61 format(f25.17)
        do i=1,Vol
          read(9,60,advance="no") SegmentNum(i)
          read(9,60,advance="no") TailPoint(i)
          do j=1,MaxKink
            read(9,60,advance="no") SegmentState(i,j)
            read(9,60,advance="no") PrevSite(i,j)
            read(9,60,advance="no") NextSite(i,j)
            read(9,60,advance="no") NeighSite(i,j)
            read(9,60,advance="no") NeighVertexNum(i,j)
            read(9,61,advance="no") KinkTime(i,j)
            do k=1,nnb
              read(9,60,advance="no") NearestSite(i,j,k)
            enddo
          enddo
        enddo
        close(9)
    else
        do k=1,Vol
            SegmentNum(k)=1;
            SegmentState(k,1)=1
            SegmentState(k,2)=-1
            KinkTime(k,1)=0.0
            KinkTime(k,2)=beta
            NeighSite(k,1)=0
            NeighSite(k,2)=0
            PrevSite(k,1)=0
            NextSite(k,1)=2
            NearestSite(k,1,:)=1
            PrevSite(k,2)=1
            NextSite(k,2)=0
            NearestSite(k,2,:)=2
            do i=3,MaxKink-1
                NextSite(k,i)=i+1
                PrevSite(k,i)=i-1
                NearestSite(k,i,:)=0
            enddo
            PrevSite(k,MaxKink)=MaxKink-1
            NextSite(k,MaxKink)=0
            NearestSite(k,MaxKink,:)=0
            !TailPoint points to the first empty node
            TailPoint(k)=3
        enddo
    endif
    
    do i = 1, Dim
	if(L(i)<=2) then
	    write(6,*) 'L must bigger than 2'
	endif

	if(L(i)-L(i)/2*2==1) then
	    write(6,*) 'Give me even L,please!'
	endif

	if(L(i)>MxL) then
	  write(6,*) 'L < MxL?';                              stop
	endif
    enddo

    !if(Ntoss==0) then
      !write(6,*) 'to compare simulation and measurement &
      !&  time, please set Ntoss/=0';                      stop
    !endif

    if((NBlck>MxBlck).or.(NBlck<MnBlck)) then
      write(6,*) 'MnBlck <= NBlck <=MxBlck?';             stop
    endif

    !if((NBlck>200).and.(NBlck/=MxBlck)) then
      !write(6,*) '"NBlck>200" is supposed for extensive &
      !& simulation. "NBlck=MxBlk" is suggested!';         stop
    !endif

    if(Dim==2) then
	write(6,40) L(1),L(2)
	40 format(' Simulation is on L=',i6,'*',i6,2x,'square lattice.')
    else if(Dim==3) then
	write(6,46) L(1),L(2),L(3),LatticeName
	46 format(' Simulation is on L=',i6,'*',i6,'*',i6,2x,a10,1x,'lattice.')
    endif

    write(6,45) Beta
    45 format(' Simulation is on Beta=',f12.8,2x,'Time Line.')

    write(6,41) JJ, QQ*2.0, ExternalField
    41 format(' Coupling Jz is',f12.8, ', Jx is', f12.8, ', H is', f12.8)

    write(6,42) Ntoss*NBlck
    42 format(' Will throw away    ',i10,2x,'steps ')

    write(6,43) Nsamp*NBlck
    43 format(' Will simulate      ',i10,2x,'steps ')

    write(6,44) Nswee      
    44 format(' A step contains    ',i10,2x,'sweeps')


    !-- integer probability ------------------------------------------

    !-- set random-number generator-------------------------------------
    call set_RNG
    
    !-- measurement initialization -----------------------------------
    Norm = 1.d0/Vol/beta

    Nk = 2
    allocate(Momentum(1:Nk, 1:Dim))
    Momentum(1, :) = (/0.d0, 0.d0, 0.d0/)
    Momentum(2, :) = (/0.d0, 0.d0, 2.d0*Pi/)

    allocate(Quan(1:NObs_b))
    allocate(Obs(1:NObs, 1:NBlck));       Obs = 0.d0
    allocate(Ave(1:NObs), Dev(1:NObs), Cor(1:NObs))
    allocate(StaticCorr(1:Vol*NSub))
    allocate(DynamicalCorr(1:Nk, 0:MxOmega))
    StaticCorr = 0.d0
    DynamicalCorr = 0.d0

    allocate(ReExp(0:NPhase-1))
    allocate(ImExp(0:NPhase-1))

    do i = 0, NPhase-1
	phase = i/(NPhase*1.d0)*2.d0*Pi
	ReExp(i) = dcos(phase)
	ImExp(i) = dsin(phase)
    enddo
    return
  END SUBROUTINE carlo
  !===================================================================

  SUBROUTINE def_lattice
    implicit none
    integer :: x, y, z, k, xm, ym, zm, xp, yp, zp, tmp
    integer :: ix, iy, iz, ixp, iyp, izp, ixm, iym, izm, vec(3), site(4), itmp
    !!!!************DEFINITION FOR LATTICES****************************
    if(LatticeName=='Cubic') then
	Sub(:) = 1
	if(Dim==2) then
	  k = 0
	  do y = 1,L(2)
	    do x = 1,L(1)
	      k = k + 1

	      xm =  -1 ;      if(x== 1) xm = L(1)-1
	      xp =  +1 ;      if(x==L(1)) xp = 1-L(1)

	      ym =  -L(1);      if(y== 1) ym = Vol-L(1)
	      yp =  +L(1);      if(y==L(2)) yp = L(1)-Vol

	      Ngs(1,k) = k+xp;    Ngs(4,k) = k+xm
	      Ngs(2,k) = k+yp;    Ngs(3,k) = k+ym

		    !            2
		    !            |
		    !       4---   ---1
		    !            |
		    !            3
	    enddo
	  enddo 

	  Back(1, 1) = 4; Back(4, 1) = 1;
	  Back(2, 1) = 3; Back(3, 1) = 2;

	  !-- auxillary variables to measure wrapping probability-----------
	  dr(:,:,:)=0
	  dr(1,1, 1) = 1;  dr(4,1, 1) =-1; dr(2,1, 1) = 0; dr(3,1, 1) = 0
	  dr(1,1, 2) = 0;  dr(4,1, 2) = 0; dr(2,1, 2) = 1; dr(3,1, 2) =-1
	else if(Dim==3) then
	  k = 0
	  do z = 1,L(3)
	    do y = 1,L(2)
	      do x = 1,L(1)
		k = k + 1

		xm =  -1 ;      if(x== 1) xm = L(1)-1
		xp =  +1 ;      if(x==L(1)) xp = 1-L(1)

		ym =  -L(1);      if(y== 1) ym = L(1)*L(2)-L(1)
		yp =  +L(1);      if(y==L(2)) yp = L(1)-L(1)*L(2)

		zm =  -L(1)*L(2);      if(z== 1) zm = Vol-L(1)*L(2)
		zp =  +L(1)*L(2);      if(z==L(3)) zp = L(1)*L(2)-Vol

		Ngs(1,k) = k+xp;    Ngs(4,k) = k+xm
		Ngs(2,k) = k+yp;    Ngs(5,k) = k+ym
		Ngs(3,k) = k+zp;    Ngs(6,k) = k+zm
	    enddo
       	  enddo 
	enddo

        Back(1, 1) = 4; Back(4, 1) = 1;
        Back(2, 1) = 5; Back(5, 1) = 2;
        Back(3, 1) = 6; Back(6, 1) = 3;

	!-- auxillary variables to measure wrapping probability-----------
	dr(:,:,:) = 0
	dr(1,1, 1) = 1;   dr(4,1, 1) = -1
	dr(2,1, 2) = 1;   dr(5,1, 2) = -1
	dr(3,1, 3) = 1;   dr(6,1, 3) = -1
      endif
    else if(LatticeName=='Pyrochlore') then
      !RealVectors
      LatticeVector(1,:) = (/0.0, 1.0, 1.0/)
      LatticeVector(2,:) = (/1.0, 0.0, 1.0/)
      LatticeVector(3,:) = (/1.0, 1.0, 0.0/)
      SubVector(1,:) = (/0.0, 0.0, 0.0/)
      SubVector(2,:) = (/0.0, 0.5, 0.5/)
      SubVector(3,:) = (/0.5, 0.0, 0.5/)
      SubVector(4,:) = (/0.5, 0.5, 0.0/)

      do i = 1, 4
	site(i) = i
      enddo
      do iz=1,L(3)
        do iy=1,L(2)
       	  do ix=1,L(1)
	    RealVector(site(1),:)=(ix-1)*LatticeVector(1,:)+(iy-1)*LatticeVector(2,:)+(iz-1)*LatticeVector(3,:)
	    RealVector(site(2),:)=RealVector(site(1),:)+SubVector(2,:)
	    RealVector(site(3),:)=RealVector(site(1),:)+SubVector(3,:)
	    RealVector(site(4),:)=RealVector(site(1),:)+SubVector(4,:)
	    do i = 1, 4
		site(i) = site(i)+4
	    enddo
	  enddo
        enddo
      enddo

      !directions
      Back(1, 1) = 4; Back(4, 2) = 1;
      Back(2, 1) = 4; Back(4, 3) = 2;
      Back(3, 1) = 4; Back(4, 4) = 3;
      Back(4, 1) = 1; Back(1, 2) = 4;
      Back(5, 1) = 1; Back(1, 3) = 5;
      Back(6, 1) = 1; Back(1, 4) = 6;

      Back(2, 2) = 5; Back(5, 3) = 2;
      Back(3, 2) = 5; Back(5, 4) = 3;

      Back(5, 2) = 2; Back(2, 3) = 5;
      Back(6, 2) = 2; Back(2, 4) = 6;

      Back(3, 3) = 6; Back(6, 4) = 3;
      Back(6, 3) = 3; Back(3, 4) = 6;

      !lattice connections
      do i = 1, 4
        site(i) = i
	Sub(site(i)) = i
      enddo

      do iz=1,L(3)
	do iy=1,L(2)
	  do ix=1,L(1)
            ixp = ix+1
	    ixm = ix-1
	    iyp = iy+1
	    iym = iy-1
	    izp = iz+1
	    izm = iz-1

            if(ixp>L(1)) ixp = ixp-L(1)
	    if(ixm==0) ixm = L(1)
	    if(iyp>L(2)) iyp = iyp-L(2)
	    if(iym==0) iym = L(2)
	    if(izp>L(3)) izp = izp-L(3)
	    if(izm==0) izm = L(3)

	    do i = 1, 3
	      Ngs(i, site(1)) = site(i+1)
	      Ngs(back(i, 1),site(i+1)) = site(1)
	    enddo

    	    call GetSite(itmp, (/ixp, iy, iz/), 1)
	    Ngs(1, site(2))  = itmp
	    Ngs(back(1, 2), itmp) = site(2)

     	    do i = 2, 3
	      Ngs(i, site(2)) = site(i+1)
	      Ngs(back(i, 2), site(i+1)) = site(2)
	    enddo

    	    call GetSite(itmp, (/ix, iyp, iz/), 1)
	    Ngs(1, site(3)) =  itmp
	    Ngs(back(1, 3), itmp) = site(3)

    	    call GetSite(itmp, (/ixm, iyp, iz/), 2)
	    Ngs(2, site(3)) = itmp
	    Ngs(back(2, 3), itmp) = site(3)

	    Ngs(3, site(3)) = site(4)
	    Ngs(back(3, 3), site(4)) = site(3)

    	    call GetSite(itmp, (/ix, iy, izp/), 1)
	    Ngs(1, site(4)) = itmp
	    Ngs(back(1, 4), itmp) = site(4)

    	    call GetSite(itmp, (/ixm, iy, izp/), 2)
	    Ngs(2, site(4)) = itmp
	    Ngs(back(2, 4), itmp) = site(4)

	    call GetSite(itmp, (/ix, iym, izp/), 3)
	    Ngs(3, site(4)) = itmp
	    Ngs(back(3, 4), itmp) = site(4)

	    do i = 1, NSub
		site(i) = site(i)+4
		Sub(site(i)) = i
	    enddo
	  enddo
	enddo
      enddo
      
      !-- auxillary variables to measure wrapping probability-----------
      dr(:,:,:) = 0
      dr(1,1, 1) = 0;   dr(1,1, 2) =  1;   dr(1,1, 3) =  1
      dr(4,1, 1) = 0;   dr(4,1, 2) = -1;   dr(4,1, 3) = -1 
      dr(2,1, 1) = 1;   dr(2,1, 2) =  0;   dr(2,1, 3) =  1
      dr(5,1, 1) =-1;   dr(5,1, 2) =  0;   dr(5,1, 3) = -1 
      dr(3,1, 1) = 1;   dr(3,1, 2) =  1;   dr(3,1, 3) =  0    
      dr(6,1, 1) =-1;   dr(6,1, 2) = -1;   dr(6,1, 3) =  0     

      dr(1,2, 1) = 0;   dr(1,2, 2) =  1;   dr(1,2, 3) =  1
      dr(4,2, 1) = 0;   dr(4,2, 2) = -1;   dr(4,2, 3) = -1 
      dr(2,2, 1) = 1;   dr(2,2, 2) = -1;   dr(2,2, 3) =  0
      dr(5,2, 1) =-1;   dr(5,2, 2) =  1;   dr(5,2, 3) =  0 
      dr(3,2, 1) = 1;   dr(3,2, 2) =  0;   dr(3,2, 3) = -1
      dr(6,2, 1) =-1;   dr(6,2, 2) =  0;   dr(6,2, 3) =  1 

      dr(1,3, 1) = 1;   dr(1,3, 2) =  0;   dr(1,3, 3) =  1
      dr(4,3, 1) =-1;   dr(4,3, 2) =  0;   dr(4,3, 3) = -1 
      dr(2,3, 1) = 1;   dr(2,3, 2) = -1;   dr(2,3, 3) =  0
      dr(5,3, 1) =-1;   dr(5,3, 2) =  1;   dr(5,3, 3) =  0 
      dr(3,3, 1) = 0;   dr(3,3, 2) =  1;   dr(3,3, 3) = -1
      dr(6,3, 1) = 0;   dr(6,3, 2) = -1;   dr(6,3, 3) =  1 

      dr(1,4, 1) = 1;   dr(1,4, 2) =  1;   dr(1,4, 3) =  0
      dr(4,4, 1) =-1;   dr(4,4, 2) = -1;   dr(4,4, 3) =  0 
      dr(2,4, 1) = 1;   dr(2,4, 2) =  0;   dr(2,4, 3) = -1
      dr(5,4, 1) =-1;   dr(5,4, 2) =  0;   dr(5,4, 3) =  1 
      dr(3,4, 1) = 0;   dr(3,4, 2) =  1;   dr(3,4, 3) = -1
      dr(6,4, 1) = 0;   dr(6,4, 2) = -1;   dr(6,4, 3) =  1 
    endif
  END SUBROUTINE def_lattice

  SUBROUTINE GetVector(Site, Vector, Sub)
    implicit none
    integer, intent(in) :: Site
    integer, intent(out) :: Vector(1:Dim), Sub
    integer :: tmp
    if(LatticeName=='Cubic') then
	if(Dim==2) then
	    Vector(2) = (Site-1)/L(1)+1
	    Vector(1) = Site-(Vector(2)-1)*L(1)
	    Sub = 1
	else if(Dim==3) then
	    Vector(3) = (Site-1)/L(1)/L(2)+1
	    Vector(2) = (Site-1-(Vector(3)-1)*L(1)*L(2))/L(1)+1
	    Vector(1) = Site-(Vector(3)-1)*L(1)*L(2)-(Vector(2)-1)*L(1)
	    Sub = 1
	endif
    else if(LatticeName=='Pyrochlore') then
	Sub = mod(Site, NSub)
	if(Sub==0) Sub = 4
	tmp = (Site-1)/4
	Vector(3) = tmp/L(1)/L(2)+1
	Vector(2) = (tmp-(Vector(3)-1)*L(1)*L(2))/L(1)+1
	Vector(1) = tmp-(Vector(3)-1)*L(1)*L(2)-(Vector(2)-1)*L(1)+1
    endif
  END SUBROUTINE

  SUBROUTINE GetRealVector(Site, Vector)
    implicit none
    integer, intent(in) :: Site
    double precision, intent(out) :: Vector(1:Dim)
    integer :: i, sub, coord(1:Dim)
    if(LatticeName=='Pyrochlore') then
	call GetVector(Site, coord, sub)
	Vector(:) = 0.0
	do i = 1, Dim
	    Vector = Vector+coord(i)*LatticeVector(i,:)
	enddo
	Vector = Vector+ SubVector(sub,:)
    endif
    return
  END SUBROUTINE GetRealVector

  SUBROUTINE GetSite(Site, Vector, Sub)
    implicit none
    integer, intent(out) :: Site
    integer, intent(in) :: Vector(1:Dim), Sub
    if(LatticeName=='Cubic') then
	if(Dim==2) then
	    Site = Vector(1)+(Vector(2)-1)*L(1)
	else if(Dim==3) then
	    Site = (Vector(3)-1)*L(1)*L(2)+(Vector(2)-1)*L(1)+Vector(1)
	endif
    else if(LatticeName=='Pyrochlore') then
	Site = ((Vector(3)-1)*L(1)*L(2)+(Vector(2)-1)*L(1)+Vector(1)-1)*NSub+Sub
    endif
  END SUBROUTINE

  !==============Simulation ==========================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE monte
        implicit none
        integer :: i,j,k
        logical :: flag
        integer :: site, ver
        double precision :: x

        do i=1,Nswee
            dWR(:)=0
	    dWt=0
	    !print *, Number, "create_worm"
            call create_worm(flag)
            if(flag) then
                do while(flag)
                    Number=Number+1
                    call choose_Ira()
                    x=rn()
                    if(x<P1) then
			!print *, Number, "move_worm"
                        call move_worm()
                    else if(x<P2) then
			!print *, Number, "create_kink"
                        call create_kink()
                    else if(x<P3) then
			!print *, Number, "delete_kink"
                        call delete_kink()
                    else
			!print *, Number, "remove_worm"
                        call annihilate_worm(flag)
                    endif
                enddo
            endif
	    do k = 1, Dim
		WR(k)=WR(k)+(dWR(k)/L(k))
	    enddo
            Wt=Wt+dWt
        enddo
        return
    END SUBROUTINE monte

    logical FUNCTION is_in_same_segment(IraVertex,IraSite,MashaVertex,MashaSite)
        implicit none
        integer ::     IraVertex
        integer ::     MashaVertex
        integer ::  IraSite,MashaSite
        if(IraVertex/=MashaVertex) then
            is_in_same_segment=.false.
        else if(IraSite/=NextSite(MashaVertex,MashaSite).and.IraSite/=PrevSite(MashaVertex,MashaSite)) then
            is_in_same_segment=.false.
        else
            is_in_same_segment=.true.
        endif
        return
    END FUNCTION

    SUBROUTINE choose_Ira()
        implicit none
        integer :: vertex,site
        double precision :: time
        if(rn()<0.5) then
            vertex=IraVertex
            IraVertex=MashaVertex
            MashaVertex=vertex

            site=IraSite
            IraSite=MashaSite
            MashaSite=site

            time=IraTime
            IraTime=MashaTime
            MashaTime=time

            direction=-direction
        endif
        return
    END SUBROUTINE
    
    SUBROUTINE create_worm(flag)
        !Choose an Interval to create a worm with Ira and Masha
        implicit none
        integer  :: Site1,Site2,Vertex,temp,i
        double precision :: prevtime,nexttime,itime,mtime,time,Tmin,Tmax
        double precision :: Factor,total_E,TempTime
        logical :: flag
        Vertex=int(rn()*Vol)+1
        temp=int(SegmentNum(Vertex)*rn())
        Site1=1
        do i=1,temp
            Site1=NextSite(Vertex,Site1)
        enddo
        Site2=NextSite(Vertex,Site1)
        Tmin=KinkTime(Vertex,Site1)
        Tmax=KinkTime(Vertex,Site2)
        itime=rn()*(Tmax-Tmin)+Tmin
        mtime=rn()*(Tmax-Tmin)+Tmin
        if(itime>mtime) then
            TempTime=itime;itime=mtime;mtime=TempTime
        endif
        total_E=single_delta_E(vertex,Site1,itime,mtime)
        Factor=Comp*SegmentNum(vertex)*(Tmax-Tmin)**2
        if(rn()<Factor*exp(-total_E)) then

            if(SegmentState(Vertex,Site1)==-1) then
                direction=-1
            else
                direction=1
            endif

            IraVertex=Vertex
            MashaVertex=Vertex
            IraTime=itime
            MashaTime=mtime
            IraSite=insert_site(Vertex,Site1,IraTime)
            MashaSite=insert_site(Vertex,IraSite,MashaTime)

            flag=.true.
        else
            flag=.false.
        endif
        return
    END SUBROUTINE create_worm
    
    SUBROUTINE annihilate_worm(flag)
        implicit none
        double precision :: Factor,total_E,Tmin,Tmax,total_E1,factor2
        integer :: EndSite
        logical :: flag
        if(.not. is_in_same_segment(IraVertex,IraSite,MashaVertex,MashaSite)) then
            return
        endif
        if(IraTime<MashaTime) then
            EndSite=NextSite(MashaVertex,MashaSite)
            Tmin=KinkTime(IraVertex,PrevSite(IraVertex,IraSite))
            Tmax=KinkTime(IraVertex,EndSite)
            total_E=single_delta_E(IraVertex,IraSite,IraTime,MashaTime)
        else
            EndSite=NextSite(IraVertex,IraSite)
            Tmin=KinkTime(IraVertex,PrevSite(MashaVertex,MashaSite))
            Tmax=KinkTime(IraVertex,EndSite)
            total_E=single_delta_E(MashaVertex,MashaSite,MashaTime,IraTime)
        endif
        factor2=Factor
        Factor=1/Comp/(SegmentNum(IraVertex)-2)/(Tmax-Tmin)**2

        RunNum0(1)=RunNum0(1)+1
        if(rn()<Factor*exp(-total_E)) then

            RunNum(1)=RunNum(1)+1

            call remove_site(IraVertex,IraSite)
            call remove_site(MashaVertex,MashaSite)
            flag=.false.
        else
            flag=.true.
        endif
        return    
    END SUBROUTINE annihilate_worm
    
    SUBROUTINE move_worm()
        implicit none
        integer :: Site1,Site2,indicator,TempSite
        double precision :: Time1,Time2,tarTime,delta_Time,OldIraTime
        double precision :: total_E,factor
        indicator=1
        Site1=PrevSite(IraVertex,IraSite)
        Site2=NextSite(IraVertex,IraSite)
        if(Site1==1) then
            TempSite=PrevSite(IraVertex,2)
            Time1=KinkTime(IraVertex,TempSite)
            Time2=KinkTime(IraVertex,Site2)
            delta_Time=Time2+beta-Time1
            tarTime=rn()*delta_Time+Time1-beta
            if(tarTime<0) then 
                indicator=2
                tarTime=tarTime+beta
                total_E=single_delta_E(IraVertex,1,0.0d0,IraTime)+ &
                    & single_delta_E(IraVertex,TempSite,tarTime,beta)
            endif
        else if(Site2==2) then
            TempSite=NextSite(IraVertex,1)
            Time1=KinkTime(IraVertex,Site1)
            Time2=KinkTime(IraVertex,TempSite)
            delta_Time=Time2+beta-Time1
            tarTime=rn()*delta_Time+Time1
            if(tarTime>beta) then
                indicator=3
                tarTime=tarTime-beta
                total_E=single_delta_E(IraVertex,IraSite,IraTime,beta)+single_delta_E(IraVertex,1,0.d0,tarTime)
            end if
        else
            Time1=KinkTime(IraVertex,Site1)
            Time2=KinkTime(IraVertex,Site2)
            delta_Time=Time2-Time1
            tarTime=rn()*delta_Time+Time1
        endif
        if(indicator==1) then
            if(tarTime<IraTime) then
                total_E=single_delta_E(IraVertex,Site1,TarTime,IraTime)
                indicator=0
             else
                total_E=single_delta_E(IraVertex,IraSite,IraTime,TarTime)
             endif
        endif
        RunNum0(2)=RunNum0(2)+1
        if(rn()<exp(-total_E)) then

            RunNum(2)=RunNum(2)+1
            
            OldIraTime=IraTime
            IraTime=tarTime
            if(indicator==0) then
                !Ira moves in the same interval ,the simplest case
                call move_site(IraVertex,IraSite,IraTime,1)
            else if(indicator==1) then
                call move_site(IraVertex,IraSite,IraTime,0)
            else if(indicator==2)  then
                !Ira moves from the first interval to the last
                dWt=dWt-direction

                SegmentState(IraVertex,1)=SegmentState(IraVertex,IraSite)
                SegmentState(IraVertex,2)=-SegmentState(IraVertex,1)
                call remove_site(IraVertex,IraSite)
                IraSite=insert_site(IraVertex,TempSite,IraTime)
            else
                !Ira moves from the last interval to the first
                dWt=dWt+direction

                SegmentState(IraVertex,1)=-SegmentState(IraVertex,IraSite)
                SegmentState(IraVertex,2)=-SegmentState(IraVertex,1)
                call remove_site(IraVertex,IraSite)
                IraSite=insert_site(IraVertex,1,IraTime)
            endif
        endif
    END SUBROUTINE move_worm
    
    SUBROUTINE create_kink()
        implicit none
        integer :: TarVertex,OldIraVertex,OldIraSite
        integer :: Site0,Site1,Site2,newsite
        double precision    :: TargetTime,NewKinkTime
        integer :: IsPast,num,temp
        double precision    :: Factor,total_E
        double precision    :: weight

        num=int(rn()*nnb)+1
        TarVertex=Ngs(num,IraVertex)
        Site0=NearestSite(IraVertex,IraSite,num)
        if(SegmentState(IraVertex,IraSite)==SegmentState(TarVertex,Site0)) then
            IsPast=1
            Site2=Site0
            Site1=PrevSite(IraVertex,IraSite)
        else
            IsPast=0
            Site2=NextSite(TarVertex,Site0)
            Site1=NextSite(IraVertex,IraSite)
        endif
        !Site0,Site2 on TarVertex,Site1 on IraVertex

        TargetTime=find_nearest_time(IraVertex,num,IraSite,IsPast)
        !TargetTime=compare(KinkTime(IraVertex,Site1),KinkTime(TarVertex,Site2),IsPast)
        NewKinkTime=rn()*(IraTime-TargetTime)+TargetTime

        if(IsPast==1) then
            total_E=double_delta_E(IraVertex,num,Site1,Site0,NewKinkTime,IraTime)
            weight=QQ
        else
            total_E=double_delta_E(IraVertex,num,IraSite,Site0,IraTime,NewKinkTime)
            weight=QQ
        endif
        Factor=real(nnb)*abs(IraTime-TargetTime)*weight

        RunNum0(3)=RunNum0(3)+1
        if(rn()<Factor*exp(-total_E)) then

            RunNum(3)=RunNum(3)+1
            dWR(:)=dWR(:)+direction*dr(num,Sub(IraVertex),:)

            OldIraVertex=IraVertex
            OldIraSite=IraSite
            IraVertex=TarVertex
            !Must move old Ira To NewKinkTime first,then NearestSite of New Ira
            !will point to Ira

            call move_site(OldIraVertex,OldIraSite,NewKinkTime,IsPast)
            !Move IraSite to new kink

            if(IsPast==1) then
                Site2=insert_site(TarVertex,Site0,NewKinkTime)
                !Notice we use Site0=Site1 here
                IraSite=insert_site(TarVertex,Site2,IraTime)
            else
                IraSite=insert_site(TarVertex,Site0,IraTime)
                !Notice we use Site0=PrevSite(TarVertex,Site1) here
                Site2=insert_site(TarVertex,IraSite,NewKinkTime)
            endif

            NeighVertexNum(OldIraVertex,OldIraSite)=num
            NeighVertexNum(IraVertex,Site2)=Back(num, Sub(OldIraVertex))
            !Set the neighbor vertex of the new kink

            NeighSite(OldIraVertex,OldIraSite)=Site2
            NeighSite(IraVertex,Site2)=OldIraSite
            !Set the neighbor site of the new kink

            !IraTime will remain the same here.
        endif
        return
    END SUBROUTINE create_kink

    SUBROUTINE delete_kink()
        implicit none
        integer :: IsPast,i,num,OldIraVertex,OldIraSite
        integer :: TargetSite,Site1,Site2,TarVertex,NSite
        double precision    :: TargetTime,OldKinkTime
        double precision    :: Factor,total_E,weight
        IsPast=int(rn()*2)
        TargetSite=nearest_site(IraVertex,IraSite,IsPast)
        if(TargetSite==1 .or. TargetSite==2) return
        if(IraVertex==MashaVertex .and. TargetSite==MashaSite) return
        num=NeighVertexNum(IraVertex,TargetSite)
        NSite=NeighSite(IraVertex,TargetSite)
        TarVertex=Ngs(num,IraVertex)
        TargetTime=find_nearest_time(IraVertex,num,IraSite,IsPast)
        if(IsPast==1) then
            if(TargetTime>KinkTime(IraVertex,TargetSite)) return
        else
            if(TargetTime<KinkTime(IraVertex,TargetSite)) return
        endif
        TargetTime=find_nearest_time(IraVertex,num,TargetSite,IsPast)

        Site1=nearest_site(IraVertex,TargetSite,IsPast)
        Site2=nearest_site(TarVertex,NeighSite(IraVertex,TargetSite),IsPast)
        OldKinkTime=KinkTime(IraVertex,TargetSite)
        if(IsPast==1) then
            total_E=double_delta_E(IraVertex,num,TargetSite,    &
               &  NSite,OldKinkTime,IraTime)
        else
            total_E=double_delta_E(IraVertex,num,IraSite,     &
              &  PrevSite(TarVertex,NSite),IraTime,OldKinkTime)
        endif

        weight=QQ
        Factor=1.0/nnb/abs(IraTime-TargetTime)/weight

        RunNum0(4)=RunNum0(4)+1
        if(rn()<Factor*exp(-total_E)) then              

            RunNum(4)=RunNum(4)+1
            dWR(:)=dWR(:)+direction*dr(num,Sub(IraVertex),:)

            OldIraVertex=IraVertex
            OldIraSite=IraSite
            IraVertex=TarVertex
            IraSite=NeighSite(OldIraVertex,TargetSite)
            call remove_site(OldIraVertex,OldIraSite) 
            call remove_site(OldIraVertex,TargetSite)
            call move_site(IraVertex,IraSite,IraTime,1-IsPast)
        endif
        return
    END SUBROUTINE delete_kink

    integer FUNCTION insert_site(CurrentVertex,CurrentSite,NewSiteTime)
        implicit none
        integer,intent(in):: CurrentSite,CurrentVertex
        integer :: NewSite,CNSite,i,Temp,NVertex
        double precision :: NewSiteTime

        if(SegmentNum(CurrentVertex)>=MaxKink-2) then
	    print *, SegmentNum(CurrentVertex), "Kink number exceeds the upper limit!"
            stop
        endif

        NewSite=TailPoint(CurrentVertex)
        TailPoint(CurrentVertex)=NextSite(CurrentVertex,NewSite)

        PrevSite(CurrentVertex,NewSite)=CurrentSite
        CNSite=NextSite(CurrentVertex,CurrentSite)
        NextSite(CurrentVertex,NewSite)=CNSite
        KinkTime(CurrentVertex,NewSite)=NewSiteTime

        PrevSite(CurrentVertex,CNSite)=NewSite
        NextSite(CurrentVertex,CurrentSite)=NewSite

        SegmentNum(CurrentVertex)=SegmentNum(CurrentVertex)+1
        do i=1,nnb
            Temp=NearestSite(CurrentVertex,CNSite,i)
            NVertex=Ngs(i,CurrentVertex)
            do while(KinkTime(NVertex,Temp)>NewSiteTime)
                Temp=PrevSite(NVertex,Temp)
            enddo
            NearestSite(CurrentVertex,NewSite,i)=Temp
        enddo

        call change_nearest(CurrentVertex,CNSite,CurrentSite,NewSite)
        SegmentState(CurrentVertex,NewSite)=-SegmentState(CurrentVertex,CurrentSite);
        insert_site=NewSite
        return
    END FUNCTION
    
    SUBROUTINE remove_site(CurrentVertex,CurrentSite)
        implicit none
        integer :: CurrentVertex,CurrentSite,BeginSite,EndSite

        BeginSite=PrevSite(CurrentVertex,CurrentSite)
        EndSite=NextSite(CurrentVertex,CurrentSite)

        NextSite(CurrentVertex,BeginSite)=EndSite 
        PrevSite(CurrentVertex,EndSite)=BeginSite

        call change_nearest(CurrentVertex,EndSite,CurrentSite,BeginSite)

        NextSite(CurrentVertex,CurrentSite)=TailPoint(CurrentVertex)
        TailPoint(CurrentVertex)=CurrentSite
        SegmentNum(CurrentVertex)=SegmentNum(CurrentVertex)-1
        return
    END SUBROUTINE remove_site

    SUBROUTINE move_site(CurrentVertex,CurrentSite,NewSiteTime,IsPast)
        implicit none
        integer :: CurrentVertex,CurrentSite,IsPast
        double precision :: NewSiteTime
        integer :: i,Temp,NVertex,PSite
        if(IsPast==1) then
            do i=1,nnb
                NVertex=Ngs(i,CurrentVertex)
                Temp=NearestSite(CurrentVertex,CurrentSite,i)
                do while(KinkTime(NVertex,Temp)>NewSiteTime)
                    NearestSite(NVertex,Temp,Back(i, Sub(CurrentVertex)))=CurrentSite
                    Temp=PrevSite(NVertex,Temp)
                enddo
                NearestSite(CurrentVertex,CurrentSite,i)=Temp
            enddo
        else
            PSite=PrevSite(CurrentVertex,CurrentSite)
            do i=1,nnb
                NVertex=Ngs(i,CurrentVertex)
                Temp=NearestSite(CurrentVertex,CurrentSite,i)
                if(NearestSite(NVertex,Temp,Back(i, Sub(CurrentVertex)))/=CurrentSite) then
                    Temp=NextSite(NVertex,NearestSite(CurrentVertex,CurrentSite,i))
                endif
                do while(KinkTime(NVertex,Temp)<NewSiteTime)
                    NearestSite(NVertex,Temp,Back(i, Sub(CurrentVertex)))=PSite
                    Temp=NextSite(NVertex,Temp)
                enddo
                NearestSite(CurrentVertex,CurrentSite,i)=PrevSite(NVertex,Temp)
            enddo
        endif
        KinkTime(CurrentVertex,CurrentSite)=NewSiteTime
        return
    end SUBROUTINE move_site

    double precision FUNCTION compare(Time1,Time2,IsPast)
        implicit none
        double precision :: Time1,Time2
        integer :: IsPast
        if(IsPast==1) then
            if(Time1>Time2) then
                compare=Time1
            else
                compare=Time2
            endif            
        else 
            if(Time1<Time2) then
                compare=Time1
            else
                compare=Time2
            endif
        endif
        return
    END FUNCTION compare

    double precision FUNCTION find_nearest_time(CurVertex,num,CurSite,IsPast)
        implicit none
        integer :: CurVertex,num,CurSite,IsPast
        integer :: NeiVertex,NeiSite
        integer :: Vertex,Site,i,TempSite
        double precision :: Time,TempTime,TarTime
        Time=KinkTime(CurVertex,CurSite)
        NeiVertex=Ngs(num,CurVertex)
        NeiSite=NearestSite(CurVertex,CurSite,num)
        if(IsPast==1) then
            TarTime=0.0
            NeiSite=NextSite(NeiVertex,NeiSite)
            do i=1,nnb
                Vertex=Ngs(i,CurVertex)
                Site=NextSite(Vertex,NearestSite(CurVertex,CurSite,i))
                TempTime=KinkTime(Vertex,Site)
                do while(TempTime>=Time)
                    Site=PrevSite(Vertex,Site)
                    TempTime=KinkTime(Vertex,Site)
                enddo
                if(TempTime>TarTime) TarTime=TempTime

                Vertex=Ngs(i,NeiVertex)
                Site=NearestSite(NeiVertex,NeiSite,i)
                if(NeiSite/=2) then
                    Site=NextSite(Vertex,Site)
                endif
                TempTime=KinkTime(Vertex,Site)
                do while(TempTime>=Time)
                    Site=PrevSite(Vertex,Site)
                    TempTime=KinkTime(Vertex,Site)
                enddo
                if(TempTime>TarTime) TarTime=TempTime
            enddo
        else
            TarTime=beta
            do i=1,nnb
                Vertex=Ngs(i,CurVertex)
                Site=NearestSite(CurVertex,CurSite,i)
                TempTime=KinkTime(Vertex,Site)
                do while(TempTime<=Time)
                    Site=NextSite(Vertex,Site)
                    TempTime=KinkTime(Vertex,Site)
                enddo
                if(TempTime<TarTime) TarTime=TempTime

                Vertex=Ngs(i,NeiVertex)
                Site=NearestSite(NeiVertex,NeiSite,i)
                TempTime=KinkTime(Vertex,Site)
                do while(TempTime<=Time)
                    Site=NextSite(Vertex,Site)
                    TempTime=KinkTime(Vertex,Site)
                enddo
                if(TempTime<TarTime) TarTime=TempTime
            enddo
        endif
        find_nearest_time=TarTime 
        return
    end FUNCTION find_nearest_time

    SUBROUTINE change_nearest(Vertex,EndSite,OldSite,NewSite)
        implicit none
        integer ::Vertex,OldSite,NewSite,EndSite,NVertex,tarSite
        integer ::i
        double precision ::Time
        Time=KinkTime(Vertex,NewSite)
        do i=1,nnb
            NVertex=Ngs(i,Vertex)
            tarSite=NearestSite(Vertex,EndSite,i)
            if(NearestSite(NVertex,tarSite,Back(i,Sub(Vertex)))==EndSite) then
                tarSite=PrevSite(NVertex,tarSite)
            endif
            do while(KinkTime(NVertex,tarSite)>=Time .and. &
                &  NearestSite(NVertex,tarSite,Back(i,Sub(Vertex)))==OldSite .and. tarSite/=1)
                NearestSite(NVertex,tarSite,Back(i,Sub(Vertex)))=NewSite
                tarSite=PrevSite(NVertex,tarSite)
            enddo
        enddo
        return
    END SUBROUTINE

    integer FUNCTION nearest_site(CurrentVertex,CurrentSite,IsPast)
        implicit none
        integer :: CurrentVertex,CurrentSite
        integer ::IsPast
        if(IsPast/=1) then
            nearest_site=NextSite(CurrentVertex,CurrentSite)
        else
            nearest_site=PrevSite(CurrentVertex,CurrentSite)
        endif 
        return       
    END FUNCTION

    integer FUNCTION get_state(Vertex,BeginSite,Time)
        implicit none
        integer :: Vertex,BeginSite,TarSite
        double precision :: Time
        TarSite=BeginSite
        do while(KinkTime(Vertex,TarSite)<Time)
            TarSite=NextSite(Vertex,TarSite)
        enddo
        get_state=-SegmentState(Vertex,TarSite)
        return
    end FUNCTION get_state

    double precision FUNCTION single_delta_E(CurrentVertex,BeginSite,BeginTime,EndTime)
    !Measure the Energy between BeginTime and EndTime on Vertex
    !BeginSite: As we need to locate the nearest site of BeginTime,so we begin
    !the search from BeginSite(Namely,KinkTime(BeginSite)<BeginTime)
        implicit none
        integer :: CurrentVertex,BeginSite
        integer :: NSite,NVertex
        double precision :: Energy,BeginTime,EndTime
        integer :: i,Vertex(3),Site(3)
        Energy=0.0
        do i=1,nnb
            NSite=NearestSite(CurrentVertex,BeginSite,i)
            NVertex=Ngs(i,CurrentVertex)
            Energy=Energy+delta_E(NVertex,NSite,BeginTime,EndTime)
        enddo
        single_delta_E=-2.0*SegmentState(CurrentVertex,BeginSite)*Energy
	single_delta_E=single_delta_E+potential_E(CurrentVertex,BeginSite,BeginTime,EndTime)
    END FUNCTION single_delta_E

    double precision FUNCTION double_delta_E(Vertex,NeighNum,Site1,Site2,BeginTime,EndTime)
    !Measure the energy between BeginTime and EndTime on Vertex and
    !Ngs(NeighNum,Vertex)
    !Site1: the site on Vertex, where the search begin
    !Site2: the site on Ngs(NeighNum,Vertex),where the search begin
        integer :: Vertex,Vertex1,Vertex2,Site1,Site2,NSite,NVertex
        integer :: NeighNum
        double precision  :: Energy1,Energy2
        double precision  :: BeginTime,EndTime
        integer :: i
        Energy1=0.0
	Energy2=0.0
        Vertex1=Vertex
        Vertex2=Ngs(NeighNum,Vertex)
        do i=1,nnb
            NVertex=Ngs(i,Vertex1)
            if(NVertex==Vertex2) cycle
            NSite=NearestSite(Vertex1,Site1,i)
            Energy1=Energy1+delta_E(NVertex,NSite,BeginTime,EndTime)
        enddo
        Energy1=Energy1*SegmentState(Vertex1,Site1)
        do i=1,nnb
            NVertex=Ngs(i,Vertex2)
            if(NVertex==Vertex1) cycle
            NSite=NearestSite(Vertex2,Site2,i)
            Energy2=Energy2+delta_E(NVertex,NSite,BeginTime,EndTime)
        enddo
        Energy2=Energy2*SegmentState(Vertex2,Site2)
        double_delta_E=-2.0*(Energy1+Energy2)
	double_delta_E=double_delta_E+potential_E(Vertex1,Site1, BeginTime, EndTime)
	double_delta_E=double_delta_E+potential_E(Vertex2,Site2, BeginTime, EndTime)
    END FUNCTION double_delta_E

    double precision FUNCTION delta_E(NVertex,NSite,BeginTime,EndTime)
        integer :: NVertex,NSite,TarSite
        double precision :: BeginTime,EndTime,TarTime,PrevTime
        double precision :: Energy
        Energy=0
        if(NSite==2) then
            TarSite=2
        else
            TarSite=NextSite(NVertex,NSite)
            do while(KinkTime(NVertex,TarSite)<BeginTime .and. TarSite/=2)
                TarSite=NextSite(NVertex,TarSite)
            enddo
        endif
        TarTime=KinkTime(NVertex,TarSite)
        PrevTime=BeginTime
        do while(TarTime<EndTime .and. TarSite/=2)
            Energy=Energy-SegmentState(NVertex,TarSite)*(TarTime-PrevTime)
            PrevTime=TarTime
            TarSite=NextSite(NVertex,TarSite)
            TarTime=KinkTime(NVertex,TarSite)
        enddo
        Energy=Energy-SegmentState(NVertex,TarSite) &
            &  *(EndTime-PrevTime)
        delta_E=JJ*Energy
        !For ferromagnetic case, there will be a minus sign here
    END FUNCTION delta_E

    double precision FUNCTION potential_E(Vertex,Site,BeginTime,EndTime)
        integer :: Vertex,Site
        double precision :: BeginTime,EndTime
        double precision :: Energy
        Energy=-2.0*SegmentState(Vertex,Site)*(EndTime-BeginTime)
        potential_E=-ExternalField*Energy* get_Sign(Vertex)
    END FUNCTION potential_E

    double precision FUNCTION get_Sign(Vertex)
        integer :: Vertex
	integer :: i, vec(1:Dim)
	double precision :: sgn
	get_Sign = 1.0
	if(LatticeName=="Cubic") then
	    call GetVector(Vertex, vec, i)
	    sgn = 0.0
	    do i = 1, Dim
		sgn = sgn + vec(i)
	    enddo
	    get_Sign = (-1.0)**sgn
	endif
    END FUNCTION get_Sign
  !===================================================================


  !==============Measurement =========================================
  SUBROUTINE measure
    implicit none
    integer :: i, j 
    double precision :: TotalKinks
    logical :: flag

    !--Observables definition ----------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    flag = .true.
    do i = 1, Dim
	flag = flag .and. WR(i)==0
    enddo
    if(Wt==0 .and. flag) then
        W0=1;W1=0
    else
        W0=0;W1=1
    endif
    Quan( 1)= potential_energy()
    Quan( 2)= kink_number(TotalKinks)
    Quan( 3)= (Quan(1)-Quan(2))/beta/Vol
    Quan( 4)= Quan(3)**2
    Quan( 5)= KOmegaCorrelator(Momentum(1,:), 0)
    Quan( 6)= KOmegaCorrelator(Momentum(1,:), MXOmega/2)
    Quan( 7)= KOmegaCorrelator(Momentum(2,:), 0)
    Quan( 8)= KOmegaCorrelator(Momentum(2,:), MXOmega/2)
    Quan( 9)= (Wt)**2
    Quan(10)= Quan( 9)**2
    Quan(11)= Quan(10)*Vol  !Spatial Susceptibility 
    Quan(12)= SM()**2       !Staggerd Magnetization square
    Quan(13)= Quan(12)**2
    Quan(14)= Quan(13)*Vol  !Spatial Susceptibility 
    Quan(15)= Correlator(1,1)
    Quan(16)= Correlator(1,2)
    Quan(17)= Correlator(1,Vol/2)
  END SUBROUTINE measure

  SUBROUTINE measure_Corr
    implicit none
    integer :: i, j 
    do j = 1, NSub
	do i = 1, Vol
	    StaticCorr(i+(j-1)*Vol) = StaticCorr(i+(j-1)*Vol)+Correlator(j, i)
	enddo
    enddo

    do i = 0, MxOmega
	do j = 1, Nk
	    DynamicalCorr(j, i) = DynamicalCorr(j, i) + KOmegaCorrelator(Momentum(j,:), i)
	enddo
    enddo
  END SUBROUTINE measure_Corr

  double precision FUNCTION potential_energy()
    implicit none
    double precision :: BeginTime,EndTime,temp1
    integer :: i,j,BeginSite,NVertex,NSite,EndSite
    double precision :: Energy1
    Energy1=0.0
    do i=1,Vol
        BeginSite=1
        EndSite=NextSite(i,1)
        temp1=0.0
        do while(BeginSite/=2)
            BeginTime=KinkTime(i,BeginSite)
            EndTime=KinkTime(i,EndSite)
            do j=1,nnb
                NSite=NearestSite(i,BeginSite,j)
                NVertex=Ngs(j,i)
                temp1=temp1+SegmentState(i,BeginSite)*delta_E(NVertex,NSite,BeginTime,EndTime)
            enddo
            BeginSite=EndSite
            EndSite=NextSite(i,EndSite)
        enddo
        Energy1=Energy1+temp1/2.0
    enddo
    potential_energy=Energy1
  end FUNCTION
  
  double precision FUNCTION kink_number(TotalKinks)
    implicit none
    integer :: i
    double precision :: TotalKinks
    kink_number=0.0
    TotalKinks=0.0
    do i=1,Vol
        kink_number=kink_number+SegmentNum(i)
    enddo
    TotalKinks=(kink_number-Vol)/2
    kink_number=(kink_number-Vol)/2
  end FUNCTION kink_number

  SUBROUTINE winding_number()
      implicit none
      integer dist(1:Dim),dt
      integer Vertex,i,j,coord
      integer Site
      dist(:)=0
      dt=0
      do Vertex=1,Vol
          dt=dt+SegmentState(Vertex,1)
      enddo
      if(LatticeName=='Cubic') then
	  do i=1,Dim
	      do coord = 1, L(i)
		  if(i==1) then
		      Vertex = coord
		  else if(i==2) then
		      Vertex = (coord-1)*L(1)+1
		  else if(i==3) then
		      Vertex = (coord-1)*L(1)*L(2)+1
		  endif
		  Site=NextSite(Vertex,1)
		  do while(Site/=2)
		      if(dr(NeighVertexNum(Vertex,Site),Sub(Vertex),i)==-1) then
			  dist(i)=dist(i)-SegmentState(Vertex,Site)
		      endif
		      Site=NextSite(Vertex,Site)
		  enddo
	      enddo
	  enddo
      else if(LatticeName=='Pyrochlore') then
	  do i=1,Dim
	      do coord = 1, L(i)
		  if(i==1) then
		      call GetSite(Vertex, (/coord, 1, 1/), 1)
		  else if(i==2) then
		      call GetSite(Vertex, (/1, coord, 1/), 1)
		  else if(i==3) then
		      call GetSite(Vertex, (/1, 1, coord/), 1)
		  endif
		  Site=NextSite(Vertex,1)
		  do while(Site/=2)
		      do j=1, Dim
			  if(dr(NeighVertexNum(Vertex,Site),Sub(Vertex),j)==-1) then
			      dist(j)=dist(j)-SegmentState(Vertex,Site)
			  endif
		      enddo
		      Site=NextSite(Vertex,Site)
		  enddo
	      enddo
	  enddo
      endif

      WindR(:)=dist(:)
      WindT=dt/2
  end SUBROUTINE winding_number

  !!!!Staggered magnetization
  double precision FUNCTION SM()
    implicit none
    double precision :: M, Sgn
    integer :: Vertex,R(1:Dim), sub,  i
    M=0.0
    do Vertex = 1, Vol
	M=M+SegmentState(Vertex,1)*get_Sign(Vertex)
    enddo
    SM=M/Vol
  end FUNCTION SM

  double precision FUNCTION Correlator(Vertex1,Vertex2)
    !SzSz zero frequency correlator
    implicit none
    double precision :: Sz
    integer :: Vertex1,Vertex2
    integer :: EndSite,BeginSite
    double precision :: BeginTime, EndTime
    Sz=0.0
    Correlator=0.0
    BeginSite=1
    EndSite=NextSite(Vertex2,BeginSite)
    do while(BeginSite/=2)
        BeginTime=KinkTime(Vertex2,BeginSite)
        EndTime=KinkTime(Vertex2,EndSite)
        Sz=Sz+SegmentState(Vertex2,BeginSite)*(EndTime-BeginTime)
        BeginSite=EndSite
        EndSite=NextSite(Vertex2,EndSite)
    enddo
    Correlator=Sz*SegmentState(Vertex1,1)/Beta/4.0
  end FUNCTION Correlator

  SUBROUTINE FrequencyCorrelator(Vertex, omega, ReCorr, ImCorr)
    !SzSz non-zero frequency correlator
    implicit none
    integer :: Vertex,omega
    double precision :: ReCorr, ImCorr
    double precision :: ReSz, ImSz
    double precision :: domega, tmp
    integer :: iphase
    integer :: EndSite,BeginSite
    double precision :: BeginTime, EndTime
    ReSz = 0.d0
    ImSz = 0.d0
    BeginSite=1
    EndSite=NextSite(Vertex,BeginSite)
    domega = omega*2.0*Pi/Beta
    do while(BeginSite/=2)
        BeginTime=KinkTime(Vertex,BeginSite)
        EndTime=KinkTime(Vertex,EndSite)
	if(omega==0) then
	    tmp = SegmentState(Vertex, BeginSite)*(EndTime-BeginTime)
	    ReSz = ReSz + tmp
	    ImSz = ImSz + tmp
	else
	    tmp = SegmentState(Vertex, BeginSite)*(EndTime-BeginTime)/domega
	    iphase = mod(int((domega*(EndTime-BeginTime))/2.0/Pi*(NPhase*1.d0)), NPhase)
	    ReSz = ReSz + tmp*ImExp(iphase)
	    ImSz = ImSz + tmp*ReExp(iphase)
	endif
        BeginSite=EndSite
        EndSite=NextSite(Vertex,EndSite)
    enddo
    ReCorr=ReSz/Beta/2.0
    ImCorr=ImSz/Beta/2.0
  end SUBROUTINE FrequencyCorrelator

  double precision FUNCTION KOmegaCorrelator(k, omega)
    !SzSz non-zero frequency correlator in momentum space
    implicit none
    integer :: omega
    integer :: i, j, l
    double precision :: k(1:Dim), vec(1:Dim)
    double precision :: phase, ReCorr, ImCorr
    integer :: iphase
    double precision :: ReKCorr, ImKCorr
    ReKCorr = 0.d0
    ImKCorr = 0.d0
    do j = 1, Vol
	call GetRealVector(j, vec)
	phase = 0.0
	do l = 1, Dim
	    phase = phase + k(l)*vec(l)
	enddo
	call FrequencyCorrelator(j, omega, ReCorr, ImCorr)
	iphase = mod(int(phase/2.d0/Pi*(NPhase*1.d0)), NPhase)
	ReKCorr = ReKCorr +ReCorr*ReExp(iphase) +ImCorr*ImExp(iphase)
	ImKCorr = ImKCorr +ImCorr*ReExp(iphase) -ReCorr*ImExp(iphase)
    enddo
    KOmegaCorrelator = (ReKCorr**2.d0+ImKCorr**2.d0)/(Vol**2.0)
  end FUNCTION KOmegaCorrelator

  double precision FUNCTION Susceptibility_total()
    implicit none
    double precision :: M,Sz,BeginTime,EndTime
    integer :: Vertex,R(1:Dim),sub, Sgn
    integer :: BeginSite,EndSite
    M=0.0
    do Vertex = 1, Vol
	Sz=0.0
	BeginSite=1
	EndSite=NextSite(Vertex,BeginSite)
	do while(BeginSite/=2)
	    BeginTime=KinkTime(Vertex,BeginSite)
	    EndTime=KinkTime(Vertex,EndSite)
	    Sz=Sz+SegmentState(Vertex,BeginSite)*(EndTime-BeginTime)
	    BeginSite=EndSite
	    EndSite=NextSite(Vertex,EndSite)
	enddo
	call GetVector(Vertex, R, sub)
	Sgn = 1
	do i = 1, Dim
	    Sgn = Sgn + R(i)
	enddo
	M=M+Sz*(-1)**(Sgn)
    enddo
    Susceptibility_total=M**2/beta/Vol
  end FUNCTION Susceptibility_total
  
  !==============Calculate composite observables =====================
  !! THIS IS PROJECT-INDEPENDENT 
  !! call in 'stat_alan'
  SUBROUTINE cal_Obs_comp
    implicit none
    integer :: k
    double precision  :: nor

    !-- calculate the average ----------------------------------------
    Ave(NObs_b+1:NObs) = 0.d0

    !-- Q1=<C1>^2/<C1^2> ---------------------------------------------
    Ave(NObs_b+1)=1-Ave(10)/Ave( 9)**2/3
    Ave(NObs_b+2)=1-Ave(13)/Ave(12)**2/3
    Ave(NObs_b+3)=Ave(4)-Ave(3)**2

    !-- Obs(j,k) series ----------------------------------------------
    Obs(NObs_b+1,:)=1-Obs(10,:)/Obs( 9,:)**2/3
    Obs(NObs_b+2,:)=1-Obs(13,:)/Obs(12,:)**2/3
    Obs(NObs_b+3,:)=Obs(4,:)-Obs(3,:)**2

   return
  END SUBROUTINE cal_Obs_comp
  !===================================================================
  !*******************************************************************
  !            End of PROJECT-DEPENDENT part 
  !*******************************************************************


  !*******************************************************************
  !            Beginning of PROJECT-INDEPENDENT part 
  !*******************************************************************
  !============== Write to files =====================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE write2file 
    implicit none
    integer       :: i, j, k, Nwri
    double precision :: err

    open(8,file=file1, access='append') 
    write(8, *);          write(6,*)

    if(Dim==2) then
	write(8,40) ident, L(:), beta, JJ*4.0, Q*4.0, Seed, TotSamp
	40 format(a8,2i6,3f10.6,i12,i8)
    else if(Dim==3) then
	write(8,46) ident, L(:), beta, JJ*4.0, Q*4.0, Seed, TotSamp
	46 format(a8,3i6,3f10.6,i12,i8)
    endif

    do j = 1, Nobs
      write(8,41) j, Ave(j), Dev(j), Cor(j)
      41 format(i4,2f25.15,f12.5)
      write(6,41) j, Ave(j), Dev(j), Cor(j)
    enddo

    do j=1, 4
      write(8,43) RunName(j),RunNum(j)/Number, RunNum(j)/RunNum0(j)
    enddo
    43 format(a13,2f12.5)

    close(8)


    if(prt) return
    Nwri = NObs_b;                   if(Nwri>5) Nwri = 7
    write(6,*)
    do k = 1, NBlck
       write(6,42) k,(Obs(j,k),j=1,Nwri)
       42 format(i6,7f16.8) 
    end do

    return
  END SUBROUTINE write2file
  !===================================================================

  SUBROUTINE midwrite2file(iblck)
    implicit none
    integer :: iblck,i, j, start
    open (8,file=file4, access='append') 
    start=iblck-NSave+1
    do i=start,iblck
        write(8,*)
        write(8,50) ident, L(1), beta, JJ*4.0, Q*16.0, Seed, NSamp, i
        50 format(a8,i6,3f10.6,i12,i8,i8)
        do j = 1, NObs_b
          write(8,51) j,Obs(j,i)
          51 format(i4,f25.15)
        enddo
    enddo
    close(8)
  END SUBROUTINE midwrite2file
  !===================================================================

  SUBROUTINE write2file_corr(iblck)
    implicit none
    integer       :: i, j, k, Nwri, iblck
    double precision :: err,nor

    nor  = 1.d0/(iblck*1.d0)
    StaticCorr(:) = nor/(NSamp*1.d0)*(NmeasCorr*1.d0)*StaticCorr(:)
    DynamicalCorr(:,:) = nor/(NSamp*1.d0)*(NmeasCorr*1.d0)*DynamicalCorr(:,:)
    open(9,file=file2) 
    write(9, *) "{ 'Correlations': [ ["
    do j = 1, Vol
      write(9,47) StaticCorr(j), ', '
      47 format(f25.15, a2)
    enddo
    write(9, *) "], ["
    do j = Vol+1, 2*Vol
      write(9,47) StaticCorr(j), ', '
    enddo
    write(9, *) "], ["
    do j = 2*Vol+1, 3*Vol
      write(9,47) StaticCorr(j), ', '
    enddo
    write(9, *) "], ["
    do j = 3*Vol+1, 4*Vol
      write(9,47) StaticCorr(j), ', '
    enddo
    write(9, *) "]]} "
    close(9)

    open (10,file=file3) 
    write(10, *) "k=", Momentum(1,:)
    write(10,*) 0.d0, DynamicalCorr(1, 0), Dev(5)*sqrt(NmeasCorr*1.d0)
    err = Dev(6)*sqrt(NmeasCorr*1.d0)
    do j = 1, MxOmega
      write(10,*) 2.d0*Pi*(j*1.d0)/Beta, DynamicalCorr(1, j), err
    enddo
    write(10, *) 
    write(10, *) "k=", Momentum(2, :)
    write(10,*) 0.d0, DynamicalCorr(2, 0), Dev(7)*sqrt(NmeasCorr*1.d0)
    err = Dev(8)*sqrt(NmeasCorr*1.d0)
    do j = 1, MxOmega
      write(10,*) 2.d0*Pi*(j*1.d0)/Beta, DynamicalCorr(2, j), err
    enddo
    write(10, *) 
    close(10)
    return
  END SUBROUTINE write2file_corr

  !===================================================================
  SUBROUTINE saveconfig
    implicit none
    integer  :: i,j,k
    open(9,file=file_config)
    50 format(i6)
    51 format(f25.17)
    do i=1,Vol
      write(9,50,advance="no") SegmentNum(i)
      write(9,50,advance="no") TailPoint(i)
      do j=1,MaxKink
        write(9,50,advance="no") SegmentState(i,j)
        write(9,50,advance="no") PrevSite(i,j)
        write(9,50,advance="no") NextSite(i,j)
        write(9,50,advance="no") NeighSite(i,j)
        write(9,50,advance="no") NeighVertexNum(i,j)
        write(9,51,advance="no") KinkTime(i,j)
        do k=1,nnb
          write(9,50,advance="no") NearestSite(i,j,k)
        enddo
      enddo
    enddo
    close(9)
  end SUBROUTINE

  !==============Collect data ========================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE coll_data(iblck)
    implicit none
    integer, intent(in) :: iblck
    integer             :: j
    do j = 1, NObs_b
      Obs(j,iblck) = Obs(j,iblck)+ Quan(j)
    enddo
  END SUBROUTINE coll_data 
  !===================================================================

  !==============Normalize by Nsamp ==================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE norm_Nsamp(iblck)
    implicit none
    integer, intent(in) :: iblck
    integer             :: j
    double precision    :: nor
    nor = 1.d0/(Nsamp*1.d0)
    do j = 1, NObs_b
      Obs(j,iblck) = Obs(j,iblck)*nor
    enddo
  END SUBROUTINE norm_Nsamp
  !===================================================================

  !==============Statistics ==========================================
  !! THIS IS PROJECT-INDEPENDENT 
  !! In a way learned from Alan
  SUBROUTINE stat_alan 
    implicit none
    integer          :: i, j, k, k0
    double precision :: devn, devp, nor
    double precision, allocatable :: Aux(:)

    nor  = 1.d0/(NBlck*1.d0)

    ! -- calculate average -------------------------------------------
    do j = 1, NObs_b
      Ave(j) = nor*Sum(Obs(j,1:NBlck))
    enddo

    Coarsen: do
      prt = .true.
      ! -- calculate error and t=1 correlation for basics obs.--------
      DO j = 1, NObs_b
        devp = 0.d0;  Cor(j) = 0.d0;  Dev(j) = 0.d0
        do k = 1,  NBlck
          devn   = Obs(j,k)-Ave(j)
          Dev(j) = Dev(j)+devn*devn
          Cor(j) = Cor(j)+devn*devp
          devp   = devn
        enddo 
        Dev(j)   = Dev(j)*nor;        Cor(j) = Cor(j)*nor
        if(Dev(j)>eps)                Cor(j) = Cor(j)/Dev(j)
        Dev(j)   = dsqrt(Dev(j)/(NBlck-1.d0))
        if(dabs(Cor(j))>tol) prt = .false.
      ENDDO 


      IF(prt)                         EXIT Coarsen 
      IF(NBlck<64)    THEN
        prt = .false.;                EXIT Coarsen 
      ENDIF

      ! -- coarsen blocking ------------------------------------------
      NBlck = NBlck/2;    nor=nor*2.d0
      DO j = 1, NObs
        k0 = 1
        do k   = 1, NBlck
          Obs(j,k) = (Obs(j,k0)+Obs(j,k0+1))*0.5d0
          k0 = k0 +2
        enddo 
      ENDDO 
    enddo Coarsen 

    ! -- define auxillary variables and average of composite obs.-----
    call cal_Obs_comp

    ! -- calculate error and t=1 correlation for composite obs.-----
    do j = 1+NObs_b, NObs
      devp = 0.d0;  Cor(j) = 0.d0;  Dev(j) = 0.d0
      DO k = 1,  NBlck
        devn   = Obs(j,k)-Ave(j)  
        Dev(j) = Dev(j)+devn*devn
        Cor(j) = Cor(j)+devn*devp
        devp   = devn
      ENDDO
      Dev(j)   = Dev(j)*nor;        Cor(j) = Cor(j)*nor
      IF(Dev(j)>eps)                Cor(j) = Cor(j)/Dev(j)
      Dev(j)   = dsqrt(Dev(j)/(NBlck-1.d0))
    enddo
    return
  END SUBROUTINE stat_alan
  !===================================================================

  !===============Shift register random number generator =============
  !  very long period sequential version
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE set_RNG
    implicit none
    integer :: i_r,k_r,k1_r
    integer :: iseed

    nrannr = mxrn
    iseed  = iabs(Seed)+1
    k_r    = 3**18+2*iseed
    k1_r   = 1313131*iseed
    k1_r   = k1_r-(k1_r/mod2)*mod2

    do i_r = 1, len1
      k_r  = k_r *mult
      k1_r = k1_r*mul2
      k1_r = k1_r-(k1_r/mod2)*mod2
      ir1(i_r) = k_r+k1_r*8193
    enddo

    do i_r = 1, len2
      k_r  = k_r *mult
      k1_r = k1_r*mul2
      k1_r = k1_r-(k1_r/mod2)*mod2
      ir2(i_r) = k_r+k1_r*4099
    enddo

    do i_r = 1, len1
      inxt1(i_r) = i_r+1
    enddo
    inxt1(len1) = 1
    ipnt1 = 1
    ipnf1 = ifd1+1

    do i_r = 1, len2
      inxt2(i_r) = i_r+1
    enddo
    inxt2(len2) = 1
    ipnt2 = 1
    ipnf2 = ifd2 + 1
    return
  END SUBROUTINE set_RNG 
  !===================================================================

 !===============Calculate next random number =======================
  !! THIS IS ALMOST PROJECT-INDEPENDENT 
  double precision function rn()
  !integer function rn()
    implicit none
    integer   :: i_r, l_r, k_r
    nrannr = nrannr +1
    if(nrannr>=mxrn) then
      nrannr = 1
      do i_r= 1, mxrn
        l_r = ieor(ir1(ipnt1),ir1(ipnf1))
        k_r = ieor(ir2(ipnt2),ir2(ipnf2))
        irn(i_r) = ieor(k_r,l_r)
        ir1(ipnt1)=l_r
        ipnt1 = inxt1(ipnt1)
        ipnf1 = inxt1(ipnf1)
        ir2(ipnt2) = k_r
        ipnt2 = inxt2(ipnt2)
        ipnf2 = inxt2(ipnf2)
      enddo
    endif 
    !rn = irn(nrannr)
    rn = irn(nrannr)*tm32+0.5d0
  end function rn
  !===================================================================

  !==============Trace elapsed time ==================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE t_elapse(t0)
    implicit none
    integer, intent(in) :: t0   ! '0' for the 1st time
    double precision    :: dt
    
    call date_and_time(date, time, zone, tval)
    t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0  ! seconds 
    h_curr = tval(5)

    if(t0==0) then      ! first time
      t_elap = 0.d0;                          dt = 0.d0
    else
      dt = t_curr-t_prev;   if(h_curr<h_prev) dt = dt+24*3600.d0    ! across midnight
      t_elap = t_elap+dt
    endif
    t_prev = t_curr;      h_prev = h_curr 

    select case(t0)
    case(1)            ! initialization time
      write(6,40) dt
      40 format(/'        set up time:',f16.7,2x,'s')
    case(2)            ! equilibration
      write(6,41) dt
      41 format( ' equilibration time:',f16.7,2x,'s.')
      write(6,47) real(Nsamp)/Ntoss*dt/60
      47 format( ' I need more:',f16.7,2x,'minutes.')
      t_simu = dt/(Ntoss*NBlck*1.d0);
      t_mcmc = t_elap
    case(3)            ! for the whole mcmc chain
      t_mcmc = t_elap-t_mcmc
      t_mcmc = t_mcmc/(Nsamp*NBlck*1.d0)
      t_meas = t_mcmc-t_simu
      write(6,42) t_meas
      42 format( '   measurement time:',f16.7,2x,'s./step ')
      write(6,43) t_mcmc
      43 format( '  markov-chain time:',f16.7,2x,'s./step ')
      write(6,44) (t_meas/t_mcmc)
      44 format( ' rela. measure time:',f16.7)
    case(4)
      write(6,45) dt
      45 format(/' stat.-write   time:',f16.7,2x,'s.')
      write(6,46) t_elap/60.d0
      46 format( ' program ends at   :',f16.7,2x,'minutes')
    case default
      write(6,48) dt
      48 format(/' time:',f16.7,2x,'s.')
    end select
    return
  END SUBROUTINE t_elapse
  !===================================================================
  !*******************************************************************
  !            End of PROJECT-INDEPENDENT part 
  !*******************************************************************

END PROGRAM main
!sudo apt-get install xournal====================================================================
