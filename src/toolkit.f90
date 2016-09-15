module toolkit

  !-- contains :
  !--   . InitRandSeed, RandNorm, RandNorm_N

contains


  !--------------------------------!
  !--          Random            --!
  !--------------------------------!

  Subroutine InitRandomSeed(nbSeed)
    !- Initialise le random
    !
    implicit none
    Integer, intent(in), optional       :: nbSeed
    Integer                             :: i, n, clock
    Integer, Dimension(:), Allocatable  :: seed

    Call Random_seed(size = n)
    Allocate(seed(n))
    if (present(nbSeed)) then
       seed = nbSeed * (/ (i - 1, i = 1, n) /)
    else
       Call System_clock(COUNT=clock)
       seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    end if
    Call Random_seed(PUT = seed)
    Deallocate(seed)
  End Subroutine InitRandomSeed

  function RandNorm()
    !- Generate a Gaussian with zero mean and variance 1
    !     RandNorm : a scalar
    !
    implicit none
    Double Precision, dimension(2)     :: X_unif
    Double Precision                   :: RandNorm
    Double Precision, Parameter        :: PI = 3.14159265358979323846
    !- classic formula
    Call Random_number(X_unif)
    RandNorm = sqrt(-2*log(X_unif(1)))*cos(2*PI*X_unif(2));
  end function RandNorm

  Subroutine RandNorm_N(Noise)
    !- Generate a vector of Gaussian
    !     RandNorm_N : a N vector of Gaussian
    !
    implicit none
    Double Precision, dimension(:), intent(out)   :: Noise
    Integer                                       :: N
    Double Precision, dimension(:,:), Allocatable :: X_unif
    Double Precision, Parameter                   :: PI = 3.14159265358979323846
    !- Allocate
    N = size(Noise)
    Allocate(X_unif(N,2))
    !- the classic formula for the Gaussian
    Call Random_number(X_unif)
    Noise = sqrt(-2*log(X_unif(:,1)))*cos(2*PI*X_unif(:,2));
    !- Deallocate
    DeAllocate(X_unif)
  end Subroutine RandNorm_N



end module toolkit
