Program DM                      ! Decision-Making

  use toolkit
  use input_output_DM
  use initial_condition_DM
  use grid_interaction
  use interaction_DM
  
  implicit none
  
  !--------------  Declaration of variables  --------------!
  !--------------------------------------------------------!
  TYPE(PARAM_DM)                                :: P
  TYPE(PARAM_init)                              :: Pinit
  Integer                                       :: N
  Double Precision, Dimension(:,:), Allocatable :: X, dX
  Integer, Dimension(:), Allocatable            :: S, coopNeigh, defectNeigh
  Double Precision, Dimension(:), Allocatable   :: rateG01, rateG10
  Integer, Dimension(:), Allocatable            :: posGrid, firstParticleGrid
  Integer, Dimension(:), Allocatable            :: verletListPrev, verletListNext
  Integer                                       :: iStep,nSteps
  Integer                                       :: i0, i0PosGrid, i0PosGrid_old
  Real                                          :: start, finish
  
  !--  Lecture of the parameters   --!
  !----------------------------------!
  Call Lecture(P,Pinit)
  N = P%N
  !- Number of steps
  nSteps = floor(P%Time/P%dt + .5)
  !- warning
  Call TestParameter(P,Pinit)
  !-- Information
  Call PrintInfo(P)
  
  !--      Initialisation       --!
  !-------------------------------!
  !- Allocate
  Allocate(X(N,2),dX(N,2))
  Allocate(S(N),coopNeigh(N),defectNeigh(N))
  Allocate(rateG01(N),rateG10(N))
  If ( P%isGrid ) Then
     Allocate(posGrid(N))
     Allocate(firstParticleGrid(P%nCaseX*P%nCaseY))
     Allocate(verletListPrev(N),verletListNext(N))
  End If
  !- Initial condition
  Call InitCond(X,S,P,Pinit)
  !- Verlet list
  if ( P%isGrid ) Then
     Call InitVerletList(X,P,posGrid,firstParticleGrid,&
          verletListPrev,verletListNext)
  endif
  !- Let's count the time it takes
  Call Cpu_time(start)
  if (P%isTrajectorySave) then
     !- Write the intial condition
     Call FilePrintVector(X(:,1),"../data/particleX_"//P%strSeed,.true.,0)
     Call FilePrintVector(X(:,2),"../data/particleY_"//P%strSeed,.true.,0)
     Call FilePrintVectorInt(S,"../data/particleS_"//P%strSeed,.true.,0)
  endif
  
     
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  !-----------------------   The big loop   ----------------------!
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  Do iStep=1,nSteps

     !- A) Switch strategy                                   -!
     !--------------------------------------------------------!
     If ( P%isGrid ) Then
        Call Neighbors_MN_Fast(X,S,P, &
             posGrid,firstParticleGrid,verletListNext, &
             coopNeigh, defectNeigh)
     else
        Call Neighbors_MN_Slow(X,S,P, coopNeigh, defectNeigh)
     End If
     Call RateSwitching(coopNeigh,defectNeigh,P, rateG01, rateG10)
     Call Switching_S(S,rateG01,rateG10,P)
     
     !- B) Compute dX for each  particle                     -!
     !--------------------------------------------------------!
     If ( P%isGrid ) Then
        Call dX_Fast(X,S,P,&
             posGrid,firstParticleGrid,verletListNext, dX)
     else
        Call dX_Slow(X,S,P, dX)
     End If
     ! update position
     X = X - P%dt*dX
     
     !- C) Update the Verlet_list -!
     !-----------------------------!
     If (P%isGrid) Then
        do i0=1,N
           !-- init
           i0PosGrid     = CellNumber(X(i0,:),P)
           i0PosGrid_old = posGrid(i0)
           !- Test if i0 has changed its position on the grid
           if (i0PosGrid /= i0PosGrid_old) then
              !- We modify the Verlet_list...
              Call ListRemove(i0,i0PosGrid_old,firstParticleGrid,&
                   verletListPrev,verletListNext)
              Call ListAdd(i0,i0PosGrid,firstParticleGrid,&
                   verletListPrev,verletListNext)
              !- ...and we update posGrid
              posGrid(i0) = i0PosGrid
           end if
        end do
     end If

     !- D) Write the data and statistical analysis -!
     !----------------------------------------------!
     if (modulo(iStep,P%jumpPrint)==0 .or. iStep==nSteps) Then
        if (P%isTrajectorySave) then
           !- print particles
           Call FilePrintVector(X(:,1),"../data/particleX_"//P%strSeed,.true.,iStep)
           Call FilePrintVector(X(:,2),"../data/particleY_"//P%strSeed,.true.,iStep)
           Call FilePrintVectorInt(S,"../data/particleS_"//P%strSeed,.true.,iStep)
        end if
     end if

     !- progress...
     Call BarProgress(iStep,nSteps)
  End do
  
  !---------------------------------------------------------------!
  !---------------------   End loop in time   --------------------!
  !---------------------------------------------------------------!
  
  !- Deallocate  
  Deallocate(X,dX,S,coopNeigh,defectNeigh)
  If ( P%isGrid ) Then
     Deallocate(posGrid,firstParticleGrid,verletListPrev,verletListNext)
  End If
  Call Cpu_time(finish)
  print *,""
  print "(A24,f11.3)"," Time to compute (s)   = ",finish-start
  print *,"******************************************************"
  print *,""

  
End Program DM
