module grid_interaction

  !-- contains :
  !       CellNumber
  !       ListAdd, ListRemove
  !       InitVerletList

  use input_output_DM          ! for TYPE(PARAM_DM)
  
contains

  function CellNumber(X_i0,P)
    !- Give the position of X_i0 on the grid
    !
    implicit none
    Double Precision, Dimension(2), intent(in) :: X_i0
    TYPE(PARAM_DM), intent(in)      :: P
    Integer                                    :: i_grid, j_grid
    Double Precision                           :: dx, dy
    Integer                                    :: CellNumber

    !------!
    ! grid !
    !------!

    !  (1,nCaseY)  |  (2,nCaseY)  |     ....    |  (nCaseX-1,nCaseY)  |  (nCaseX,nCaseY)  !
    !                                                                                     !
    !       .                                                                   .         !
    !       .                                                                   .         !
    !       .                                                                   .         !
    !                                                                                     !
    !     (1,1)    |     (2,1)    |     ....    |     (nCaseX-1,1)    |    (nCaseX,1)     !
    !
    
    ! in the list   : (i_grid,j_grid)  ->  k_i0   = (j_grid-1)*nCaseX + i_grid
    !  inversely    :        k_i0      ->  i_grid = modulo(k_i0-1,nCaseX) + 1
    !                        k_i0      ->  j_grid =  floor(k_i0-1,nCaseX) + 1
    
    ! We count from the left to the right beginning from bottom.

    !- init
    dx = P%Lx/P%nCaseX
    dy = P%Ly/P%nCaseY
    
    !- Position on the grid
    i_grid = floor(X_i0(1)/dx) + 1
    j_grid = floor(X_i0(2)/dy) + 1

    !- The value
    CellNumber = (j_grid-1)*P%nCaseX + i_grid
    
  end function CellNumber


    
  Subroutine ListAdd(i0,k_i0,firstParticleGrid,verletListPrev,verletListNext)
    !- Add a new particle (i0) at the place k_i0
    !
    implicit none
    Integer, intent(in)                        :: i0, k_i0
    Integer, Dimension(:), intent(inout)       :: firstParticleGrid
    Integer, Dimension(:), intent(inout)       :: verletListPrev,verletListNext
    Integer                                    :: oldFirst
    !- Init
    oldFirst = firstParticleGrid(k_i0)
    !- Next particle after i0
    firstParticleGrid(k_i0) = i0
    verletListNext(i0)      = oldFirst
    !- Previous particle of i0 (its parent)
    verletListPrev(i0)      = 0
    if (oldFirst /= 0) Then
       verletListPrev(oldFirst) = i0
    endif
  end Subroutine ListAdd

  
  Subroutine ListRemove(i0,k_i0,firstParticleGrid,verletListPrev,verletListNext)
    !- Remove the particule i0 (at the place k_i0) from the list
    !
    implicit none
    Integer, intent(in)                        :: i0,k_i0
    Integer, Dimension(:), intent(inout)       :: firstParticleGrid
    Integer, Dimension(:), intent(inout)       :: verletListPrev,verletListNext
    Integer                                    :: i0Prev,i0Next
    !- Iinit
    i0Prev = verletListPrev(i0)
    i0Next = verletListNext(i0)
    !- Kill the links from i0
    verletListPrev(i0) = 0
    verletListNext(i0) = 0
    !----- Relink  -----!
    !-  i0Prev  -->  i0Next
    if (i0Prev /= 0) Then
       !- i0 was not the first in k_i0
       verletListNext(i0Prev)  = i0Next
    else
       !- i0 was the first in k_i0
       firstParticleGrid(k_i0) = i0Next
    end if
    !-  i0Prev  <--  i0Next
    if (i0Next /= 0) Then
       !- i0 was not the last in k_i0
       verletListPrev(i0Next)  = i0Prev
    end if
  end Subroutine ListRemove


  Subroutine InitVerletList(X,P,PosGrid,firstParticleGrid,verletListPrev,verletListNext)
    !- Initiate the verlet list(s)
    !
    implicit none
    Double Precision, Dimension(:,:), intent(in) :: X
    Type(PARAM_DM), intent(in)                   :: P
    Integer, Dimension(:), intent(out)           :: PosGrid,firstParticleGrid
    Integer, Dimension(:), intent(out)           :: verletListPrev,verletListNext
    Integer                                      :: i0, k_i0
    !- Init
    firstParticleGrid  = 0
    verletListPrev     = 0
    verletListNext     = 0
    !- For each particle
    Do i0=1,P%N
       !- Position on the grid
       k_i0 = CellNumber(X(i0,:),P)
       PosGrid(i0) = k_i0
       !- Add a new particle in the lists verletListPrev and verletListNext
       Call ListAdd(i0,k_i0,firstParticleGrid,verletListPrev,verletListNext)
    end Do
  end Subroutine InitVerletList


end module grid_interaction
