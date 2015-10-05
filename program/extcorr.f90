
  !============== Variable Module ====================================
  MODULE my_vrbls
    IMPLICIT NONE

    double precision, parameter :: eps    = 1.d-14            ! very small number
    double precision, parameter :: tol    = 0.15d0            ! tolerance for Cor
    double precision, parameter :: Pi = 3.141592653

    integer :: Lx, Dim=3, Vol, num
    double precision :: Beta
    integer :: NBlck=1024
    integer :: MaxBlck=1024
    integer :: MinBlck
    integer :: iblck,TotSamp,Seed
    double precision :: JJ,Q

    double precision,allocatable :: Momentum(:,:)
    double precision,allocatable :: StaticCorr(:,:,:), FreqCorr(:,:,:)
    double precision,allocatable :: AveStaticCorr(:,:), AveFreqCorr(:,:)
    double precision,allocatable :: DevStaticCorr(:,:), DevFreqCorr(:,:)

    character(8 ), parameter :: ident = 'hs_sqa 0'            ! identifier
    character(12), parameter :: file1 = 'hs_sqa0.dat'        ! datafile
    character(100) :: file3 
    character(100),allocatable :: corrfile(:), freqfile(:)
    !-----------------------------------------------------------------

    !-- Observables --------------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    integer, parameter :: NObs_b = 18                         ! #basic observables
    integer, parameter :: NObs_c = 3                          ! #composite observables
    integer, parameter :: NObs   = NObs_b+NObs_c              ! Total # observables
    integer, parameter :: MxOmega = 128
    integer, parameter :: Nk = 2
    integer, parameter :: NSub = 4
    !-----------------------------------------------------------------
  END MODULE my_vrbls

  !=====Main routine for bond percolation on square lattice ==========
  PROGRAM main
    use my_vrbls
    implicit none
    integer :: whichone,stat,NB,ind
    double precision :: k1, k2, k3, scorr, fcorr, err, omega
    character(100) :: filename, head
    integer :: i,j,k, flag, numstat, numfreq
    double precision :: nor
    character(8) :: string

    !print *, 'Outputfile,MinBlck'
    read  *,  num,Lx,Beta
    Vol = Lx*Lx*Lx*NSub

    allocate(corrfile(1:num))
    allocate(freqfile(1:num))

    do whichone = 1, num
      if(whichone<=10) then
          write(corrfile(whichone),"(a12,i1,a4)") "static_corr_",whichone-1,".txt"
          write(freqfile(whichone),"(a15,i1,a4)") "corr_frequency_",whichone-1,".txt"
      else
        write(corrfile(whichone),"(a12,i2,a4)")   "static_corr_",whichone-1,".txt"
        write(freqfile(whichone),"(a15,i2,a4)")   "corr_frequency_",whichone-1,".txt"
      endif
    enddo

    allocate(StaticCorr(1:NSub, 1:Vol, 1:num))
    allocate(FreqCorr(1:Nk, 0:MxOmega, 1:num))
    allocate(AveStaticCorr(1:NSub, 1:Vol))
    allocate(AveFreqCorr(1:Nk, 0:MxOmega))
    allocate(DevStaticCorr(1:NSub, 1:Vol))
    allocate(DevFreqCorr(1:Nk, 0:MxOmega))
    StaticCorr(:,:,:) = 0.d0
    FreqCorr(:,:,:) = 0.d0

    allocate(Momentum(1:Nk, 1:Dim))
    Momentum(1, :) = (/0.d0, 0.d0, 0.d0/)
    Momentum(2, :) = (/0.d0, 0.d0, 2.d0*Pi/)

    numstat = 0
    numfreq = 0
    do i = 1, num
      open (8,file=freqfile(i)) 
      ind=0
      do while(.true.)
          read(8, 40,iostat=stat) head, k1, k2, k3
          40 format(a8,3f10.6)
          if(stat/=0) exit

          ind=ind+1
          if(ind>Nk) then
              write(*,*) "Too many momentums!"
              exit
          else
              do j = 0, MxOmega
                read(8,*) omega, fcorr, err
                FreqCorr(ind, j, numstat+1) = fcorr
              enddo
              read(8, *)
          endif
      enddo
      numstat = numstat+1
      close(8)

      open(9, file=corrfile(i))
      read(9,42, iostat=stat) head
      42 format(a20)
      if(stat/=0) exit

      do k = 1, NSub
        do j = 1, Vol
          read(9,43) scorr, string
          43 format(f25.15, a3)
          StaticCorr(k, j, numfreq+1) = scorr
        enddo
        read(9, *)
      enddo
      numfreq = numfreq + 1
      close(9)
    enddo

    call stat_alan(StaticCorr(:,:,:), NSub, Vol, numstat, AveStaticCorr(:,:), DevStaticCorr(:,:))
    call stat_alan(FreqCorr(:,:,:), Nk, MxOmega+1, numfreq, AveFreqCorr(:,:), DevFreqCorr(:,:))

    call write2file()

  CONTAINS

  SUBROUTINE stat_alan(Data2d, d1, d2, num, Ave, Dev)
    implicit none
    integer          :: i, j, k, k0, d1, d2, num
    logical          :: prt
    double precision :: devn, devp, nor
    double precision :: Data2d(1:d1, 1:d2, 1:num), Ave(1:d1, 1:d2), Dev(1:d1, 1:d2)

    nor  = 1.d0/(num*1.d0)

    ! -- calculate average -------------------------------------------
    do i = 1, d1
	do j = 1, d2
	    Ave(i, j) = nor*Sum(Data2d(i, j, 1:num))
	enddo
    enddo

    ! -- calculate error and t=1 correlation for basics obs.-------- 
    DO i = 1, d1
	do j = 1, d2
	    devp = 0.d0;  Dev(i, j) = 0.d0
	    do k = 1,  num
	      devn   = Data2d(i, j, k)-Ave(i, j)
	      Dev(i, j) = Dev(i, j)+devn*devn
	    enddo 
	    Dev(i, j)   = Dev(i, j)*nor
	    Dev(i, j)   = dsqrt(Dev(i, j)/(num-1.d0))
	enddo
    ENDDO 
    return
  END SUBROUTINE stat_alan



  SUBROUTINE write2file 
    implicit none
    integer       :: i, j, k, Nwri
    open(9,file="static_corr.txt") 
    write(9, *) "{ 'Correlations': [ ["
    do j = 1, Vol
      write(9,47) AveStaticCorr(1, j), ', '
      47 format(f25.15, a2)
    enddo
    write(9, *) "], ["
    do j = 1, Vol
      write(9,47) AveStaticCorr(2, j), ', '
    enddo
    write(9, *) "], ["
    do j = 1, Vol
      write(9,47) AveStaticCorr(3, j), ', '
    enddo
    write(9, *) "], ["
    do j = 1, Vol
      write(9,47) AveStaticCorr(4, j), ', '
    enddo
    write(9, *) "]]} "
    close(9)

    open (10,file="corr_frequency.txt") 
    write(10, *) "k=", Momentum(1,:)
    do j = 0, MxOmega
      write(10,*) 2.d0*Pi*(j*1.d0)/Beta, AveFreqCorr(1, j), DevFreqCorr(1, j)
    enddo
    write(10, *) 
    write(10, *) "k=", Momentum(2, :)
    do j = 0, MxOmega
      write(10,*) 2.d0*Pi*(j*1.d0)/Beta, AveFreqCorr(2, j), DevFreqCorr(2, j)
    enddo
    write(10, *) 
    close(10)
    return
  END SUBROUTINE write2file

END PROGRAM main
!=====================================================================




