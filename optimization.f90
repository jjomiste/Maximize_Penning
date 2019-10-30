module optimizationmod

  use globalmod
  use fff_mod

implicit none

public

contains

!! Powell's method

SUBROUTINE linmin(p,xi,fret)
USE nrtype; USE nrutil, ONLY : assert_eq
USE nr, ONLY : mnbrak,brent
USE f1dim_mod
IMPLICIT NONE
REAL(SP), INTENT(OUT) :: fret
REAL(SP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
REAL(SP), PARAMETER :: TOL=1.0e-4_sp
!Given an N -dimensional point p and an N -dimensional direction xi , both vectors of length N , moves and resets p to where the fixed-name function func takes on a minimum along the direction xi from p , and replaces xi by the actual vector displacement that p was moved. Also returns as fret the value of func at the returned location p . This is actually all accomplished by calling the routines mnbrak and brent .Parameter: Tolerance passed to brent .
REAL(SP) :: ax,bx,fa,fb,fx,xmin,xx
ncom=assert_eq(size(p),size(xi),'linmin')
pcom=>p
!Communicate the global variables to f1dim.
xicom=>xi
ax=0.0
!Initial guess for brackets.
xx=1.0
call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
fret=brent(ax,xx,bx,f1dim,TOL,xmin)
if (xmin.gt.300) then !! in the case of not convergence in brent
   fret=333
else
xi=xmin*xi
end if
!Construct the vector results to return.
p=p+xi
END SUBROUTINE linmin

SUBROUTINE powell(p,xi,ftol,iter,fret)
USE nrtype; USE nrutil, ONLY :assert_eq, nrerror
!USE nr, ONLY :linmin
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: xi
INTEGER(I4B), INTENT(OUT) :: iter
REAL(SP), INTENT(IN) :: ftol
REAL(SP), INTENT(OUT) :: fret
!INTERFACE
!FUNCTION func(p)
!USE nrtype
!IMPLICIT NONE
!REAL(SP), DIMENSION(:), INTENT(IN) :: p
!REAL(SP) :: func
!END FUNCTION func
!END INTERFACE
INTEGER(I4B), PARAMETER :: ITMAX=200
REAL(SP), PARAMETER :: TINY=1.0e-25_sp

!Minimization of a function func of N variables. (func is not an argument, it is a fixed function name.) Input consists of an initial starting point p, a vector of length N;an initial N × N matrix xi whose columns contain the initial set of directions (usually the N unit vectors);a nd ftol, the fractional tolerance in the function value such that failure to decrease by more than this amount on one iteration signals doneness. On output, p is set to the best point found, xi is the then-current direction set, fret is the returned function value at p, and iter is the number of iterations taken. The routine linmin is used. Parameters: Maximum allowed iterations, and a small number.

INTEGER(I4B) :: i,ibig,n
REAL(SP) :: del,fp,fptt,t
REAL(SP), DIMENSION(size(p)) :: pt,ptt,xit
n=assert_eq(size(p),size(xi,1),size(xi,2),'powell')
fret=fff(p)

pt(:)=p(:) ! Save the initial point.
iter=0

do
iter=iter+1
fp=fret
ibig=0
del=0.0 ! Will be the biggest function decrease.
do i=1,n !Loop over all directions in the set.
xit(:)=xi(:,i) ! Copy the direction,
fptt=fret
call linmin(p,xit,fret) ! minimize along it
if (fret.gt.300) then !! on error in brent
   return
end if
if (fptt-fret > del) then !and record it if it is the largest decrease so
   del=fptt-fret !far
   ibig=i
end if
end do
if (2.0_sp*(fp-fret) <= ftol*(abs(fp)+abs(fret))+TINY) RETURN !Termination criterion.
if (iter == ITMAX) call nrerror('powell exceeding maximum iterations')
ptt(:)=2.0_sp*p(:)-pt(:) ! Construct the extrapolated point and the average direction moved. Save the old starting point.
xit(:)=p(:)-pt(:)
pt(:)=p(:)
fptt=fff(ptt) !Function value at extrapolated point.
if (fptt >= fp) cycle ! One reason not to use new direction.
t=2.0_sp*(fp-2.0_sp*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
if (t >= 0.0) cycle !Other reason not to use new direction.
call linmin(p,xit,fret) !Move to minimum of the new direction,
if (fret.gt.300) then !! on error in brent
   return
end if
xi(:,ibig)=xi(:,n) !and save the new direction.
xi(:,n)=xit(:)
end do !Back for another iteration.
END SUBROUTINE powell


!! Using amoeba. Not working very well

SUBROUTINE amoeba(p,y,ftol,func,iter)
IMPLICIT NONE
INTEGER(kind=4), INTENT(OUT) :: iter
REAL(kind=8), INTENT(IN) :: ftol
REAL(kind=8), DIMENSION(:), INTENT(INOUT) :: y
REAL(kind=8), DIMENSION(:,:), INTENT(INOUT) :: p
INTERFACE
   FUNCTION func(x)
   IMPLICIT NONE
   REAL(kind=8), dimension(:), INTENT(IN) :: x
   REAL(kind=8) :: func
   END FUNCTION func
END INTERFACE
INTEGER(kind=4), PARAMETER :: ITMAX=5000
REAL(kind=8), PARAMETER :: TINY=1.0D-10
!Minimization of the function func in N dimensions by the downhill simplex method of Nelder and Mead. The (N + 1) \times N matrix p is input. Its N + 1 rows are N -dimensional vectors that are the vertices of the starting simplex. Also input is the vector y of length N + 1, whose components must be preinitialized to the values of func evaluated at the N + 1 vertices (rows) of p ; and ftol the fractional convergence tolerance to be achieved in the function value (n.b.!). On output, p and y will have been reset to N + 1 new points all within ftol of a minimum function value, and iter gives the number of function evaluations taken. Parameters: The maximum allowed number of function evaluations, and a small number.
INTEGER(kind=4) :: ihi,ndim
REAL(kind=8), DIMENSION(size(p,2)) :: psum
call amoeba_private
CONTAINS

SUBROUTINE amoeba_private
IMPLICIT NONE
INTEGER(kind=4) :: i,ilo,inhi
REAL(kind=8) :: rtol,ysave,ytry,ytmp

IF((size(p,2) == size(p,1) - 1) .and. (size(p,2) == size(y) -1))THEN
    ndim = size(y) - 1
ELSE
    STOP 'ERROR: terminated in amoeba for inconsistent arr dimensions'
ENDIF
iter=0
psum(:)=sum(p(:,:),dim=1)
do
   ilo=iminloc(y(:))
   ihi=imaxloc(y(:))
   ytmp=y(ihi)
   y(ihi)=y(ilo)
   inhi=imaxloc(y(:))
   y(ihi)=ytmp
   rtol=2.0D0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
   if (rtol < ftol) then
      call swap_scalar(y(1),y(ilo))
      call swap_vector(p(1,:),p(ilo,:))
      RETURN
   end if
   if (iter >= ITMAX) STOP 'ERROR: ITMAX exceeded in amoeba'
   ytry=amotry(-1.0D0)
   iter=iter+1
   if (ytry <= y(ilo)) then
      ytry=amotry(2.0D0)
      iter=iter+1
   else if (ytry >= y(inhi)) then
      ysave=y(ihi)
      ytry=amotry(0.5D0)
      iter=iter+1
      if (ytry >= ysave) then
         p(:,:)=0.5D0*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
         do i=1,ndim+1
            if (i /= ilo) y(i)=func(p(i,:))
         end do
         iter=iter+ndim
         psum(:)=sum(p(:,:),dim=1)
      end if
   end if
end do
END SUBROUTINE amoeba_private

FUNCTION amotry(fac)
IMPLICIT NONE
REAL(kind=8), INTENT(IN) :: fac
REAL(kind=8) :: amotry
REAL(kind=8) :: fac1,fac2,ytry
REAL(kind=8), DIMENSION(size(p,2)) :: ptry
fac1=(1.0D0-fac)/ndim
fac2=fac1-fac
ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
ytry=func(ptry)
if (ytry < y(ihi)) then
   y(ihi)=ytry
   psum(:)=psum(:)-p(ihi,:)+ptry(:)
   p(ihi,:)=ptry(:)
end if
amotry=ytry
END FUNCTION amotry

FUNCTION imaxloc(arr)
REAL(kind=8), DIMENSION(:), INTENT(IN) :: arr
INTEGER(kind=4) :: imaxloc
INTEGER(kind=4), DIMENSION(1) :: imax
imax=maxloc(arr(:))
imaxloc=imax(1)
END FUNCTION imaxloc
FUNCTION iminloc(arr)
REAL(kind=8), DIMENSION(:), INTENT(IN) :: arr
INTEGER(kind=4) :: iminloc
INTEGER(kind=4), DIMENSION(1) :: imax
imax=minloc(arr(:))
iminloc=imax(1)
END FUNCTION iminloc
SUBROUTINE swap_scalar(a,b)
REAL(kind=8), INTENT(INOUT) :: a,b
REAL(kind=8) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_scalar
SUBROUTINE swap_vector(a,b)
REAL(kind=8), DIMENSION(:), INTENT(INOUT) :: a,b
REAL(kind=8), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_vector
END SUBROUTINE amoeba

end module optimizationmod
