module initial_condition_DM

  !-- contains :
  !       InitCond

  use toolkit                      ! for random
  use input_output_DM              ! for PARAM_DM

contains
  
  Subroutine InitCond(X,S,P,Pinit)
    !- initial condition
    implicit none
    Double Precision, Dimension(:,:), intent(out)   :: X
    Integer, Dimension(:), intent(out)              :: S
    TYPE(PARAM_DM), intent(in)                      :: P
    TYPE(PARAM_init), intent(in)                    :: Pinit
    Integer                                         :: i0,ig,jg
    Double Precision                                :: coin
    Double Precision, Dimension(2)                  :: X_i0
    logical                                         :: isInsideDomain

    if (P%isInitRand) Then
       Call InitRandomSeed()
    else
       Call InitRandomSeed(P%nbSeed)
    endif
    
    !---------------------------------------------------!
    !-------             Position              ---------!
    !---------------------------------------------------!
    Do i0=1,P%N
       select case(Pinit%initCondX)
       case (1)
          !-  Uniform distribution
          Call Random_number(coin)
          X(i0,1) = P%Lx/2 + Pinit%Lx*(coin-.5)
          Call Random_number(coin)
          X(i0,2) = P%Ly/2 + Pinit%Ly*(coin-.5)
       case(2)
          !- Gaussian distribution
          isInsideDomain = .false.
          Do while (.not. isInsideDomain)
             X_i0(1) = Pinit%xMean + Pinit%xStd*RandNorm()
             X_i0(2) = Pinit%yMean + Pinit%yStd*RandNorm()
             if (0d0<=X_i0(1) .and. X_i0(1)<P%Lx .and. 0d0<=X_i0(2) .and. X_i0(2)<P%Ly) Then
                isInsideDomain = .true.
                X(i0,:) = X_i0
             endif
          end Do
       case (3)    
          ! square lattice
          ig = int(i0/sqrt(real(P%N)))
          jg = i0-1 - ig*int(sqrt(real(P%N)))
          X(i0,:) = (/ P%Lx/3+ig*P%R_def , P%Ly/3+jg*P%R_def /)
       end select
    end Do

    !---------------------------------------------------!
    !-------             Strategy              ---------!
    !---------------------------------------------------!
    Do i0=1,P%N
       select case(Pinit%initCondS)
       case (1)
          !-  Uniform distribution
          Call Random_number(coin)
          if (coin<Pinit%ratioCoop) then
             S(i0) = 0
          else
             S(i0) = 1
          end if
       case(2)
          !- block
          if (i0<=Pinit%ratioCoop*P%N) then
             S(i0) = 0
          else
             S(i0) = 1
          end if
       end select
    end Do


  end subroutine InitCond

end module initial_condition_DM
