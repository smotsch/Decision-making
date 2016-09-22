module interaction_DM

  !-- contains :
  !       phi
  !       dX_Slow, dX_Fast
  !       Neighbors_MN_Slow, Neighbors_MN_Fast
  !       RateSwitching, Switching_S

  use input_output_DM             ! for PARAM_DM
  use grid_interaction

contains

  !------------------------------------------------------!
  !-----------               phi              -----------!
  !------------------------------------------------------!

  function phi(d,s)
    ! the 'pushing' function
    implicit none
    Double precision                   :: d
    Integer                            :: s
    Double precision                   :: phi
    phi = 0
    if (s==0 .and. d<2) then
       phi = (1/d - .5)
    end if
    if (s==1 .and. d<1) then
       phi = (1/d - 1)
    end if
  end function phi
  
  !------------------------------------------------------!
  !-----------               dX               -----------!
  !------------------------------------------------------!
  
  Subroutine dX_Slow(X,S,P,dX)
    !- Compute X' with the "slow" technic of order O(N^2)
    implicit none
    Double Precision, Dimension(:,:), intent(in)   :: X
    Integer, Dimension(:), Allocatable, intent(in) :: S
    TYPE(PARAM_DM), intent(in)                     :: P
    Double Precision, Dimension(:,:), intent(out)  :: dX
    Double Precision, Dimension(2)                 :: X_i, X_j
    Double Precision                               :: d_ij
    Integer                                        :: i,j
    !- init
    dX = 0
    !----------  Look around particle i  ----------!
    Do i=1,P%N
       X_i = X(i,:)
       !- Find the neighbor of the particle i
       Do j=i+1,P%N
          !- The jth particle
          X_j  = X(j,:)
          d_ij = norm2(X_j-X_i)
          dX(i,:) = dX(i,:) + phi(d_ij,S(i))*(X_j-X_i)
          dX(j,:) = dX(j,:) + phi(d_ij,S(j))*(X_i-X_j)
       end Do
    end Do

  end subroutine dX_Slow


  Subroutine dX_Fast(X,S,P,posGrid,firstParticleGrid,verletListNext,dX)
    !- Computation X' with the "fast" method which requires the Verlet List
    !
    implicit none
    Double Precision, Dimension(:,:), intent(in)   :: X
    Integer, Dimension(:), Allocatable, intent(in) :: S
    TYPE(PARAM_DM), intent(in)                     :: P
    Integer, Dimension(:), intent(in)              :: posGrid, firstParticleGrid
    Integer, Dimension(:), intent(in)              :: verletListNext
    Double Precision, Dimension(:,:), intent(out)  :: dX
    Double Precision, Dimension(2)                 :: X_i, X_j
    Double Precision                               :: d_ij
    Integer                                        :: i,j
    Integer                                        :: k, k_i,i_grid,j_grid,i_g,j_g
    !- init
    dX = 0
    !- For each particle
    Do i=1,P%N
       X_i = X(i,:)
       k_i   = posGrid(i)
       i_grid = modulo(k_i-1,P%nCaseX) + 1
       j_grid =      (k_i-1)/P%nCaseX  + 1
       !- Compute the sum of the velocity of its neighbors
       Do i_g=i_grid-1,i_grid+1
          Do j_g=j_grid-1,j_grid+1
             !- Is the case inside the domain?
             if (i_g>=1 .and. i_g<=P%nCaseX .and. j_g>=1 .and. j_g<=P%nCaseY) Then
                !- init
                k = (j_g-1)*P%nCaseX + i_g           ! case (i_g,j_g) in the list
                j = firstParticleGrid(k)           ! first particle in the case k
                !- while there is somebody in this box...
                Do While (j/=0)
                   if (j>i) Then
                      !- The neighbor
                      X_j  = X(j,:)
                      d_ij = norm2(X_j-X_i)
                      dX(i,:) = dX(i,:) + phi(d_ij,S(i))*(X_j-X_i)
                      dX(j,:) = dX(j,:) + phi(d_ij,S(j))*(X_i-X_j)
                   end if
                   ! Next neighbour!
                   j = verletListNext(j)
                end Do
             endif

          end Do
       end Do
    end do

  end subroutine dX_Fast

  
  Subroutine Neighbors_MN_Slow(X,S,P, coopNeigh, defectNeigh)
    ! Compute the numbers of cooperating/defecting neighbors around each particle
    implicit none
    Double Precision, Dimension(:,:), intent(in)      :: X
    Integer, Dimension(:), Allocatable, intent(in)    :: S
    TYPE(PARAM_DM), intent(in)                        :: P
    Integer, Dimension(:), intent(out)                :: coopNeigh, defectNeigh
    Double Precision, Dimension(2)                    :: X_i, X_j
    Double Precision                                  :: d_ij
    Integer                                           :: i,j
    !- init
    coopNeigh = 0
    defectNeigh = 0
    !----------  Look around particle i  ----------!
    Do i=1,P%N
       X_i = X(i,:)
       !- Find the neighbor of the particle i
       Do j=1,P%N
          !- The jth particle
          X_j  = X(j,:)
          d_ij = norm2(X_j-X_i)
          if ((S(i)==0 .and. d_ij<P%R_coop) .or. (S(i)==1 .and. d_ij<P%R_def)) then
             ! j is a neighbor of i
             if (S(j)==0) then
                coopNeigh(i) = coopNeigh(i) + 1
             else            ! S(j) ==1
                defectNeigh(i) = defectNeigh(i) + 1
             end if
          end if
       end Do
    end Do
  end Subroutine Neighbors_MN_Slow
  

  Subroutine Neighbors_MN_Fast(X,S,P,posGrid,firstParticleGrid,verletListNext, coopNeigh, defectNeigh)
    ! Compute the numbers of cooperating/defecting neighbors around each particle
    ! using a 'fast' method
    implicit none
    Double Precision, Dimension(:,:), intent(in)   :: X
    Integer, Dimension(:), Allocatable, intent(in) :: S
    TYPE(PARAM_DM), intent(in)                     :: P
    Integer, Dimension(:), intent(in)              :: posGrid, firstParticleGrid
    Integer, Dimension(:), intent(in)              :: verletListNext
    Integer, Dimension(:), intent(out)             :: coopNeigh, defectNeigh
    Double Precision, Dimension(2)                 :: X_i, X_j
    Double Precision                               :: d_ij
    Integer                                        :: i,j
    Integer                                        :: k, k_i,i_grid,j_grid,i_g,j_g
    !- init
    coopNeigh = 0
    defectNeigh = 0
    !- For each particle
    Do i=1,P%N
       !-- Init
       X_i = X(i,:)
       k_i   = posGrid(i)
       i_grid = modulo(k_i-1,P%nCaseX) + 1
       j_grid =      (k_i-1)/P%nCaseX  + 1
       !- Compute the sum of the velocity of its neighbors
       Do i_g=i_grid-1,i_grid+1
          Do j_g=j_grid-1,j_grid+1
             !- Is the case inside the domain?
             if (i_g>=1 .and. i_g<=P%nCaseX .and. j_g>=1 .and. j_g<=P%nCaseY) Then
                !- init
                k = (j_g-1)*P%nCaseX + i_g           ! case (i_g,j_g) in the list
                j = firstParticleGrid(k)           ! first particle in the case k
                !- while there is somebody in this box...
                Do While (j/=0)
                   !- The neighbor
                   X_j  = X(j,:)
                   d_ij = norm2(X_j-X_i)
                   if ((S(i)==0 .and. d_ij<P%R_coop) .or. (S(i)==1 .and. d_ij<P%R_def)) then
                      ! j is a neighbor of i
                      if (S(j)==0) then
                         coopNeigh(i) = coopNeigh(i) + 1
                      else            ! S(j) ==1
                         defectNeigh(i) = defectNeigh(i) + 1
                      end if
                   end if
                   ! Next neighbour!
                   j = verletListNext(j)
                end Do
             endif
          end Do
       end Do
    end do

  end subroutine Neighbors_MN_Fast


  Subroutine RateSwitching(coopNeigh,defectNeigh,P, rateG01, rateG10)
    ! compute the rate g_01 and g_10
    implicit none
    Integer, Dimension(:), intent(in)             :: coopNeigh, defectNeigh
    TYPE(PARAM_DM), intent(in)                    :: P
    Double Precision, Dimension(:), intent(out)   :: rateG01, rateG10
    
    ! normalized voter model
    select case (P%choiceModel)
    case(1)
       ! 1) voter model
       rateG01 = defectNeigh
       rateG10 = coopNeigh
    case(2)
       ! 2) normalized voter model
       rateG01 = 1d0*defectNeigh/(coopNeigh+defectNeigh)
       rateG10 = 1d0*coopNeigh/(coopNeigh+defectNeigh)
    case(3)
       ! 3) threshold voter model
       !rateG01 = 1d0*(defectNeigh > P%Threshold)
       !rateG10 = 1d0*(coopNeigh > P%Threshold)
       rateG01 = transfer(defectNeigh > P%Threshold, 1.0)
       rateG10 = transfer(coopNeigh > P%Threshold, 1.0)
    case(4)
       ! 4) best-response: B0>B1
       !rateG01 = 1d0*((P%a10*coopNeigh + P%a11*defectNeigh) > (P%a00*coopNeigh + P%a01*defectNeigh))
       rateG01 = transfer((P%a10*coopNeigh + P%a11*defectNeigh) > (P%a00*coopNeigh + P%a01*defectNeigh), 1.0)
       rateG10 = 1-rateG01
    case(5)
       ! 5) linear response
       rateG01 = (P%a10*coopNeigh + P%a11*defectNeigh)
       rateG10 = (P%a00*coopNeigh + P%a01*defectNeigh)
    end select
  end Subroutine RateSwitching

  
  Subroutine Switching_S(S,rateG01,rateG10,P)
    ! update the strategy
    implicit none
    Integer, Dimension(:), Allocatable, intent(inout)        :: S
    Double Precision, Dimension(:), Allocatable, intent(in)  :: rateG01, rateG10
    TYPE(PARAM_DM), intent(in)                               :: P
    Double Precision, dimension(:), Allocatable              :: rand_N
    Integer, dimension(:), Allocatable                       :: shouldSwitch
    ! allocate
    Allocate(rand_N(P%N))
    Allocate(shouldSwitch(P%N))
    Call Random_number(rand_N)
    ! test if should switch: rate is alpha*g_01 or alpha*g_10 depending on S
    shouldSwitch = (rand_N < (1-exp(-P%alpha*P%dt*( rateG01*(1-S) + rateG10*S ))))
    ! update: 0 -> 1 and 1 -> 0
    S = S + (1-S)*shouldSwitch - S*shouldSwitch
    ! same as S = S + (S==0)*shouldSwtich - (S==1)*shouldSwtich 

    Deallocate(rand_N,shouldSwitch)

  end Subroutine Switching_S


  Subroutine dX_and_Neighbors_slow(X,S,P, dX,coopNeigh,defectNeigh)
    !- Compute dX and the numbers of neighbors  with the "slow" technic of order O(N^2)
    implicit none
    Double Precision, Dimension(:,:), intent(in)   :: X
    Integer, Dimension(:), Allocatable, intent(in) :: S
    TYPE(PARAM_DM), intent(in)                     :: P
    Double Precision, Dimension(:,:), intent(out)  :: dX
    Integer, Dimension(:), intent(out)             :: coopNeigh, defectNeigh
    Double Precision, Dimension(2)                 :: X_i, X_j
    Double Precision                               :: d_ij
    Integer                                        :: i,j
    !- init
    dX = 0
    coopNeigh = 0
    defectNeigh = 0
    !----------  Look around particle i  ----------!
    Do i=1,P%N
       X_i = X(i,:)
       ! we count ourself as neighbor
       if (S(i)==0) then
          coopNeigh(i) = coopNeigh(i) + 1
       else
          defectNeigh(i) = coopNeigh(i) + 1
       end if
       !- Find the neighbor of the particle i
       Do j=i+1,P%N
          !- The jth particle
          X_j  = X(j,:)
          d_ij = norm2(X_j-X_i)
          ! update dX
          dX(i,:) = dX(i,:) + phi(d_ij,S(i))*(X_j-X_i)
          dX(j,:) = dX(j,:) + phi(d_ij,S(j))*(X_i-X_j)
          ! update neighbors i
          if ((S(i)==0 .and. d_ij<P%R_coop) .or. (S(i)==1 .and. d_ij<P%R_def)) then
             ! j is a neighbor of i
             if (S(j)==0) then
                coopNeigh(i) = coopNeigh(i) + 1
             else            ! S(j) ==1
                defectNeigh(i) = defectNeigh(i) + 1
             end if
          end if
          ! update neighbors j
          if ((S(j)==0 .and. d_ij<P%R_coop) .or. (S(j)==1 .and. d_ij<P%R_def)) then
             ! i is a neighbor of j
             if (S(i)==0) then
                coopNeigh(j) = coopNeigh(j) + 1
             else            ! S(i) ==1
                defectNeigh(j) = defectNeigh(j) + 1
             end if
          end if
       end Do
    end Do
  end Subroutine dX_and_Neighbors_slow


  Subroutine dX_and_Neighbors_fast(X,S,P,posGrid,firstParticleGrid,verletListNext, dX,coopNeigh,defectNeigh)
    !- Computation dX and the numbers of neighbors with the "fast" method
    ! which requires the Verlet List
    implicit none
    Double Precision, Dimension(:,:), intent(in)   :: X
    Integer, Dimension(:), Allocatable, intent(in) :: S
    TYPE(PARAM_DM), intent(in)                     :: P
    Integer, Dimension(:), intent(in)              :: posGrid, firstParticleGrid
    Integer, Dimension(:), intent(in)              :: verletListNext
    Double Precision, Dimension(:,:), intent(out)  :: dX
    Integer, Dimension(:), intent(out)             :: coopNeigh, defectNeigh
    Double Precision, Dimension(2)                 :: X_i, X_j
    Double Precision                               :: d_ij
    Integer                                        :: i,j
    Integer                                        :: k, k_i,i_grid,j_grid,i_g,j_g
    !- init
    dX = 0
    coopNeigh = 0
    defectNeigh = 0
    !- For each particle
    Do i=1,P%N
       X_i = X(i,:)
       k_i   = posGrid(i)
       i_grid = modulo(k_i-1,P%nCaseX) + 1
       j_grid =      (k_i-1)/P%nCaseX  + 1
       ! we count ourself as neighbor
       if (S(i)==0) then
          coopNeigh(i) = coopNeigh(i) + 1
       else
          defectNeigh(i) = coopNeigh(i) + 1
       end if
       !- Compute the sum of the velocity of its neighbors
       Do i_g=i_grid-1,i_grid+1
          Do j_g=j_grid-1,j_grid+1
             !- Is the case inside the domain?
             if (i_g>=1 .and. i_g<=P%nCaseX .and. j_g>=1 .and. j_g<=P%nCaseY) Then
                !- init
                k = (j_g-1)*P%nCaseX + i_g           ! case (i_g,j_g) in the list
                j = firstParticleGrid(k)           ! first particle in the case k
                !- while there is somebody in this box...
                Do While (j/=0)
                   if (j>i) Then
                      !- The neighbor
                      X_j  = X(j,:)
                      d_ij = norm2(X_j-X_i)
                      dX(i,:) = dX(i,:) + phi(d_ij,S(i))*(X_j-X_i)
                      dX(j,:) = dX(j,:) + phi(d_ij,S(j))*(X_i-X_j)
                      ! update neighbors i
                      if ((S(i)==0 .and. d_ij<P%R_coop) .or. (S(i)==1 .and. d_ij<P%R_def)) then
                         ! j is a neighbor of i
                         if (S(j)==0) then
                            coopNeigh(i) = coopNeigh(i) + 1
                         else            ! S(j) ==1
                            defectNeigh(i) = defectNeigh(i) + 1
                         end if
                      end if
                      ! update neighbors j
                      if ((S(j)==0 .and. d_ij<P%R_coop) .or. (S(j)==1 .and. d_ij<P%R_def)) then
                         ! i is a neighbor of j
                         if (S(i)==0) then
                            coopNeigh(j) = coopNeigh(j) + 1
                         else            ! S(i) ==1
                            defectNeigh(j) = defectNeigh(j) + 1
                         end if
                      end if
                   end if
                   ! Next neighbour!
                   j = verletListNext(j)
                end Do
             endif

          end Do
       end Do
    end do

  end subroutine dX_and_Neighbors_fast

  
end module interaction_DM
