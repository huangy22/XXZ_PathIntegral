
  !============== Variable Module ====================================
  MODULE my_vrbls
    IMPLICIT NONE

    double precision, parameter :: eps    = 1.d-14            ! very small number
    double precision, parameter :: tol    = 0.15d0            ! tolerance for Cor

    integer :: Lx
    double precision :: Beta
    integer :: NBlck=1024
    integer :: MaxBlck=1024
    integer :: MinBlck
    integer :: iblck,TotSamp,Seed
    double precision :: JJ,Q

    character(8 ), parameter :: ident = 'hs_sqa 0'            ! identifier
    character(12), parameter :: file1 = 'hs_sqa0.dat'        ! datafile
    character(100) :: file3 
    character(100) :: file4
    !-----------------------------------------------------------------

    !-- Observables --------------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    integer, parameter :: NObs_b = 23                         ! #basic observables
    integer, parameter :: NObs_c = 5                          ! #composite observables
    integer, parameter :: NObs   = NObs_b+NObs_c              ! Total # observables
    !-----------------------------------------------------------------

    !-- Statistics ---------------------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    double precision, allocatable :: Quan(:)                  ! Measured quantities
    double precision, allocatable :: Obs(:,:)                 ! 1st--#quan.  2nd--#block
    double precision, allocatable :: Ave(:), Dev(:), Cor(:)   ! average, error bars, and correlation of observables
    !-----------------------------------------------------------------
  END MODULE my_vrbls

  !=====Main routine for bond percolation on square lattice ==========
  PROGRAM main
    use my_vrbls
    implicit none
    integer :: whichone,stat,NB,ind
    character(100) :: filename
    integer :: i,j,flag

    !print *, 'Outputfile,MinBlck'
    read  *,  whichone,MinBlck

    if(whichone<10) then
        write(file3,"(a14,i1,a4)") "mid_hs_sqa0_1_",whichone,".txt"
    else if(whichone<100) then
        write(file3,"(a14,i2,a4)") "mid_hs_sqa0_1_",whichone,".txt"
    else if(whichone<1000) then
        write(file3,"(a14,i3,a4)") "mid_hs_sqa0_1_",whichone,".txt"
    endif

    allocate(Quan(1:NObs_b))
    allocate(Obs(1:NObs, 1:NBlck));       Obs = 0.d0
    allocate(Ave(1:NObs), Dev(1:NObs), Cor(1:NObs))

    open (8,file=file3) 
    ind=0
    flag=1
    do while(.true.)
        read(8, *,iostat=stat);
        if(stat/=0) exit

        ind=ind+1
        if(ind>MaxBlck) then
            write(*,*) "Too much blck!"
            exit
        else
            read(8,40) filename, Lx, beta, JJ, Q, Seed, TotSamp, NB, NBlck
            40 format(a8,i6,3f10.6,i12,i8,i8,i8)

            do j = 1, Nobs_b
	      if(j==4) then
		  read(8, *)
		  Obs(j, ind) = 0.0
	      else
		  read(8,41) i, Obs(j,ind)
		  41 format(i4,2f25.15,f12.5)
	      endif
            enddo
        endif
        if(ind>=MinBlck .and. flag==1)then
            ind=0
            flag=0
        endif
    enddo
    close(8)
    if(flag==1)stop
    NBlck=ind
    iblck=NBlck

    write(*,*) iblck

    if(ind>=1) then
        call stat_alan()
        call write2file()
    endif

  CONTAINS

  SUBROUTINE cal_Obs_comp
    implicit none
    integer :: k
    !-- calculate the average ----------------------------------------
    Ave(NObs_b+1:NObs) = 0.d0

    !-- Q1=<C1>^2/<C1^2> ---------------------------------------------
    if(abs(Ave(9))>eps) then
	Ave(NObs_b+1)=(Ave(5)+Ave(6)+Ave(7))/Ave(9)/3.0
    else
	Ave(NObs_b+1)=0
    endif
    Ave(NObs_b+2)=Ave(18)/Ave(19)
    Ave(NObs_b+3)=1.0-Ave(10)/Ave( 9)**2/3.0
    Ave(NObs_b+4)=1.0-Ave(13)/Ave(12)**2/3.0
    Ave(NObs_b+5)=(Ave(4)-Ave(3)**2)/Lx**3/4.0

    !-- Obs(j,k) series ----------------------------------------------
    do k=1,NBlck
        if(abs(Obs(9,k))>eps) then
            Obs(NObs_b+1,k)=(Obs(5,k)+Obs(6,k)+Obs(7,k))/Obs(9,k)/3.0
        else
            Obs(NObs_b+1,k)=0.0
        endif
    enddo
    Obs(NObs_b+2,:)=Obs(18,:)/Obs(19,:)
    Obs(NObs_b+3,:)=1.0-Obs(10,:)/Obs( 9,:)**2/3.0
    Obs(NObs_b+4,:)=1.0-Obs(13,:)/Obs(12,:)**2/3.0
    Obs(NObs_b+5,:)=(Obs(4,:)-Obs(3,:)**2)/Lx**3.0/4.0

   return
  END SUBROUTINE cal_Obs_comp

  SUBROUTINE stat_alan 
    implicit none
    integer          :: i, j, k, k0
    logical          :: prt
    double precision :: devn, devp, nor
    double precision, allocatable :: Aux(:)

    nor  = 1.d0/(NBlck*1.d0)
    !MinBlck=16

    ! -- calculate average -------------------------------------------
    do j = 1, NObs_b
      Ave(j) = nor*Sum(Obs(j,1:NBlck))
    enddo

    Coarsen: do
      prt = .true.
      print *, NBlck
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
      IF(NBlck<MinBlck)    THEN
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

  SUBROUTINE write2file 
    implicit none
    integer       :: i, j, k, Nwri

    open (8,file=file1, access='append') 
    write(8, *);          write(6,*)

    write(8,40) ident, Lx, beta, JJ, Q, Seed, TotSamp, NBlck, iblck
    40 format(a8,i6,3f10.6,i12,i8,i8,i8)

    do j = 1, Nobs
      write(8,41) j, Ave(j), Dev(j), Cor(j)
      41 format(i4,2f25.15,f12.5)
    enddo
    close(8)
    return
  END SUBROUTINE write2file

END PROGRAM main
!=====================================================================
