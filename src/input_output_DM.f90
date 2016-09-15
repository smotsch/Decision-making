module input_output_DM

  !-- contains :
  !--     TYPE PARAM_DM, PARAM_InitDM
  !--     Lecture, PrintInfo, TestParam
  !--     FilePrintVector, FilePrintVectorInt, BarProgress

  TYPE PARAM_DM
     Integer                            :: N
     Double Precision                   :: R_coop, R_def
     Integer                            :: choiceModel ! 1) voter model, 2) normalized, 3) threshold, 4) Best-response, 5) linear response
     Double Precision                   :: alpha
     Double Precision                   :: a00, a01, a10, a11 ! payoff matrix
     Double Precision                   :: Threshold
     Double Precision                   :: Lx, Ly
     Double Precision                   :: dt,Time
     Logical                            :: isGrid
     Logical                            :: isInitRand
     Integer                            :: nbSeed
     Character(len=80)                  :: strSeed
     Logical                            :: isTrajectorySave
     Integer                            :: jumpPrint

     Integer                            :: nCaseX, nCaseY
  end TYPE PARAM_DM

  TYPE PARAM_init
     Integer                            :: initCondX    ! square lattice, Gaussian
     Double Precision                   :: Lx,Ly
     Double Precision                   :: xMean,yMean,xVar,yVar
     Integer                            :: initCondS    ! uniform, two blocks
     Double Precision                   :: ratioCoop
  end TYPE PARAM_init


contains

  Subroutine Lecture(P,Pinit)
    !- Read the parameters of the simulation from
    !- the file PARAMETER_DM
    !
    implicit none
    TYPE(PARAM_DM), intent(out)             :: P
    TYPE(PARAM_init), intent(out)           :: Pinit
    !- temp
    character(15)                           :: temp,tpStr
    double precision                        :: R_max
    Double Precision, PARAMETER             :: PI = 3.14159265358979323846

    !----------------------------------------!
    !- 1) Digesting the parameters file...  -!
    !----------------------------------------!
    open(unit=15,file='PARAMETER_DM.txt',status='old')

    read(15,*)temp
    read(15,*)temp
    read(15,*)temp
    read(15,*)temp
    read(15,*)P%N
    read(15,*)temp
    read(15,*)P%R_coop
    read(15,*)P%R_def
    read(15,*)temp
    read(15,*)P%choiceModel
    read(15,*)temp
    read(15,*)P%alpha
    read(15,*)temp
    read(15,*)P%a00
    read(15,*)P%a01
    read(15,*)P%a10
    read(15,*)P%a11
    read(15,*)temp
    read(15,*)P%threshold
    read(15,*)temp
    read(15,*)P%Lx
    read(15,*)P%Ly
    read(15,*)temp
    read(15,*)P%dt
    read(15,*)P%Time

    read(15,*)temp
    read(15,*)temp
    read(15,*)P%isGrid
    read(15,*)temp
    read(15,*)P%isInitRand
    read(15,*)P%nbSeed
    read(15,*)temp
    read(15,*)P%isTrajectorySave
    read(15,*)temp
    read(15,*)P%jumpPrint

    close(unit=15)

    !- extra term
    R_max = max(P%R_coop,P%R_def)
    P%nCaseX  = floor( P%Lx/R_max )
    P%nCaseY  = floor( P%Ly/R_max )
    if (P%isInitRand) then
       P%strSeed = ""
    else
       write(tpStr , *) P%nbSeed
       P%strSeed = "seed"//trim(adjustl(tpStr))//"_"
    end if

    !---------------------------!
    !- 2) Digesting the IC...  -!
    !---------------------------!
    open(unit=16,file='PARAMETER_init.txt',status='old')

    read(16,*)temp
    read(16,*)temp
    read(16,*)temp
    read(16,*)temp
    read(16,*)Pinit%initCondX
    read(16,*)temp
    read(16,*)Pinit%Lx
    read(16,*)Pinit%Ly
    read(16,*)temp
    read(16,*)Pinit%xMean
    read(16,*)Pinit%yMean
    read(16,*)Pinit%xVar
    read(16,*)Pinit%yVar
    read(16,*)temp
    read(16,*)Pinit%initCondS
    read(16,*)temp
    read(16,*)Pinit%ratioCoop
    
    close(unit=16)

  end Subroutine Lecture


  Subroutine PrintInfo(P)
    !- Write the parameters in the terminal during the simulation
    !
    implicit none
    TYPE(PARAM_DM), intent(in)                      :: P
    Double Precision, PARAMETER                     :: PI = 3.14159265358979323846

    !-----  Information on the terminal  -----!
    !-----------------------------------------!
    print *,"******************************************************"
    print *,"***********  Parameters of the simulation  ***********"
    print *,"******************************************************"
    print *,"|--------  Parameters model  --------|"
    print "(A25,I8)",  " Number particles (N) = ",P%N
    print "(A25,f8.2)","   radius coop.        =  ",P%R_coop
    print "(A25,f8.2)","   radius def.         =  ",P%R_def
    select case (P%choiceModel)
    case (1)
       print *," choice model         :    Voter model"
    case (2)
       print *," choice model         :    normalized Voter model"
    case (3)
       print *," choice model         :    threshold Voter model"
    case (4)
       print *," choice model         :    best-response dynamics"
    case (5)
       print *," choice model         :    linear-response"
    end select
    print "(A25,f8.2)","  jump intensity (alpha)  =  ",P%alpha
    if (P%choiceModel==3) then
       print "(A25,f8.2)","   threshold       =  ",P%threshold
    end if
    if (P%choiceModel==4 .or. P%choiceModel==5) then
       print "(A25,f8.2)","      a00          =  ",P%a00
       print "(A25,f8.2)","      a01          =  ",P%a01
       print "(A25,f8.2)","      a10          =  ",P%a10
       print "(A25,f8.2)","      a11          =  ",P%a11
    end if
    !--  boundary shape
    print "(A25,f7.1)","         Lx            = ",P%Lx
    print "(A25,f7.1)","         Ly            = ",P%Ly
    print "(A25,f9.3)","         dt            = ",P%dt    
    print "(A25,f7.1)","     Total time        = ",P%Time

    !-- Computational stuf
    print *,"|----  Parameters computation  ------|"
    print *," Use grid             :     ",P%isGrid
    print *," Init  Random         :     ",P%isInitRand
    print *," Trajectory save      :     ",P%isTrajectorySave
    
  end subroutine PrintInfo


  Subroutine TestParameter(P,Pinit)
    !- Check if the parameters are not crazy...
    TYPE(PARAM_DM), intent(in)      :: P
    TYPE(PARAM_init), intent(in)      :: Pinit
    !- Some tests
    If (P%N>1d6) Then
       print *,"Warning : the number of particles is too large for 'FilePrint' (maximum 1 million)"
    end If
    If (Pinit%initCondX==1) then
       If (P%Lx<Pinit%Lx .or. P%Ly<Pinit%Ly) then
          print *,"Domain not large enough for the IC"
       end If
    end If
    If ( P%isGrid .and. P%nCaseX<3 ) Then
       print *,"The number of cases in x is too small (nCaseX=",P%nCaseX,")"
    end If
    If ( P%isGrid .and. P%nCaseY<3 ) Then
       print *,"The number of cases in y is too small (nCaseY=",P%nCaseY,")"
    end If
  end subroutine TestParameter


  Subroutine FilePrintVector(X,fileName,isUnFormatted,nbIteration)
    !- Write the real vector X in a file called 'fileName'.
    !- The writing can be easier formatted (i.e. human readable)
    !- or unformatted (binary format).
    implicit none
    Double Precision, Dimension(:), intent(in) :: X
    Character (len=*), intent(in)              :: fileName
    Logical, intent(in)                        :: isUnFormatted
    Integer, intent(in), optional              :: nbIteration
    Character(9)                               :: Extension
    Character (len=80)                         :: fileName_ext
    !- Init: change the fileName
    if (present(nbIteration)) Then
       write(Extension,'(I9.9)') nbIteration
       fileName_ext = trim(fileName)//trim(Extension)
    else
       fileName_ext = trim(fileName)
    endif
    if (isUnFormatted) then
       fileName_ext = trim(fileName_ext)//'.udat'
    else
       fileName_ext = trim(fileName_ext)//'.dat'
    endif
    !- Open and write
    if (isUnFormatted) then
       ! unformatted format
       Open(Unit=10,file=fileName_ext,form='UNFORMATTED')
       write(10) X
    else
       ! formatted format
30     format(100000(ES 15.7))
       Open(Unit=10,file=fileName_ext)
       write(10,30) X
    end if
    close(10)
  end Subroutine FilePrintVector

  Subroutine FilePrintVectorInt(S,fileName,isUnFormatted,nbIteration)
    !- Write the real vector X in a file called 'fileName'.
    !- The writing can be easier formatted (i.e. human readable)
    !- or unformatted (binary format).
    implicit none
    Integer, Dimension(:), intent(in)          :: S
    Character (len=*), intent(in)              :: fileName
    Logical, intent(in)                        :: isUnFormatted
    Integer, intent(in), optional              :: nbIteration
    Character(9)                               :: Extension
    Character (len=80)                         :: fileName_ext
    !- Init: change the fileName
    if (present(nbIteration)) Then
       write(Extension,'(I9.9)') nbIteration
       fileName_ext = trim(fileName)//trim(Extension)
    else
       fileName_ext = trim(fileName)
    endif
    if (isUnFormatted) then
       fileName_ext = trim(fileName_ext)//'.udat'
    else
       fileName_ext = trim(fileName_ext)//'.dat'
    endif
    !- Open and write
    if (isUnFormatted) then
       ! unformatted format
       Open(Unit=10,file=fileName_ext,form='UNFORMATTED')
       write(10) S
    else
       ! formatted format
30     format(100000(I1.1))
       Open(Unit=10,file=fileName_ext)
       write(10,30) S
    end if
    close(10)
  end Subroutine FilePrintVectorInt


  Subroutine BarProgress(i,imax)
    !- Inform the user of the progression of the simulation
    !
    implicit none
    integer          :: i,imax,k
    character(len=1) :: bar, back
    ! the subroutine
    !-init
    back = char(8)
    bar  = '='
    ! print the percentage and the bar
    write(*,'(256a1)', advance='no') (back, k =1,(30*i/imax)+9)
    write(*,'(2x,1i3,1a1,2x,1a1,256a1)', advance='no') &
         100*i/imax,'%','|', (bar, k =1,30*i/imax)
    if (i==imax) Then
       write(*,'(a)') '| done.'
    end if
  End subroutine BarProgress


end module input_output_DM
