MODULE f1dim_mod !Used for communication from linmin to f1dim.
  USE nrtype
  USE fff_mod
  INTEGER(I4B) :: ncom
  REAL(SP), DIMENSION(:), POINTER :: pcom,xicom
CONTAINS
  FUNCTION f1dim(x)
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: x
    REAL(SP) :: f1dim !Used by linmin as the one-dimensional function passed to mnbrak and brent .
!!$    INTERFACE
!!$       FUNCTION func(x)
!!$         USE nrtype
!!$         REAL(SP), DIMENSION(:), INTENT(IN) :: x
!!$         REAL(SP) :: func
!!$       END FUNCTION func
!!$    END INTERFACE
    REAL(SP), DIMENSION(:), ALLOCATABLE :: xt
    allocate(xt(ncom))
    xt(:)=pcom(:)+x*xicom(:)
    f1dim=fff(xt)
    deallocate(xt)
  END FUNCTION f1dim
    

END MODULE f1dim_mod
