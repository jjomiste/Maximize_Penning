program maximize_penning

use globalmod
use optimizationmod

implicit none

integer :: nn, ij,jk, ik
integer :: iloop1,iloop2, iloop3, iloop4, iloop5, iloop6, iloop7, iloop8
integer :: ngrid
integer (kind=4) :: iteracion
real(SP) :: froot, tolerancia,sigma,pop0,pop1
real(SP), allocatable :: xx(:,:), yy(:)
real(SP), allocatable :: pp(:)
complex(sp) :: as(1:3), bs(1:3)
integer :: method

!! if method is 1, then go to the powell method
!! if method is 2, then go to the self consistent method

method = 2

ngrid= 3 !! ngrid must be larger than 1

nn=8 !! number of variables of the function to minimize
tolerancia=1e-8

if (ngrid.le.1) then
   write(*,*) 'ngrid must be larger than 1'
   stop
end if

!! physical parameters

aij0=1.0d0 !! AI for J=0
aij1=2.0d0 !! AI for J=1
pij0=3.d0  !! PI for J=0
pij1=5.0d0 !! PI for J=1

allocate(xx(1:nn,1:nn))
allocate(yy(1:nn+1))
allocate(pp(1:nn))

!! write headings

do ij=7,10

   select case (method)
   case(1)
      write(ij,*) '#Powell Method'
   case(2)
      write(ij,*) '#Self-consistent Method'
   case default
      write(*,*)
      write(*,*) 'Select a method to obtain the optimal points'
      write(*,*) 
end select

write(ij,*) '#Columns:'
write(ij,*) '#1: \rho'
write(ij,*) '#2: \eta'
write(ij,*) '#3: \alpha_0'
write(ij,*) '#4: \alpha_1'
write(ij,*) '#5: \psi'
write(ij,*) '#6: \chi'
write(ij,*) '#7: \beta_0'
write(ij,*) '#8: \beta_1'
write(ij,*) '#9: Value of the functions'
write(ij,*) '#10: Population J=0'
write(ij,*) '#11: Population J=1'
write(ij,*) '#12: Population J=2'

end do

write(19,*) 'Coefficients'
write(19,*)

do iloop1=0, ngrid-1
do iloop2=0, ngrid-1
do iloop3=0, ngrid-1
do iloop4=0, ngrid-1
do iloop5=0, ngrid-1
do iloop6=0, ngrid-1
do iloop7=0, ngrid-1
do iloop8=0, ngrid-1

!! initialize the directional states

xx=zero

   do ij=1,nn
      xx(ij,ij)=one/10.d0
   end do
   
   !! set the initial point
   pp=one

   pp(1)=dble(iloop1)*pi*half/dble(ngrid-1) !rho
   pp(2)=dble(iloop2)*pi*half/dble(ngrid-1) !eta
   pp(3)=dble(iloop3)*two*pi/dble(ngrid-1)  !alpha0
   pp(4)=dble(iloop4)*two*pi/dble(ngrid-1)  !alpha1

   pp(5)=dble(iloop5)*pi*half/dble(ngrid-1) !psi
   pp(6)=dble(iloop6)*pi*half/dble(ngrid-1) !chi
   pp(7)=dble(iloop7)*two*pi/dble(ngrid-1) !beta0
   pp(8)=dble(iloop8)*two*pi/dble(ngrid-1) !beta1

   select case(method)

      case(1) !! using a numerical minimization method
         
         call powell(pp,xx,tolerancia,iteracion,froot)
         
         if (froot.gt.300) then
            continue
         else
            write(7,'(12E21.9)') pp(:), fff(pp), population0(pp), population1(pp), one-population0(pp)-population1(pp)
      if (population0(pp).lt.1.d-5) then
         write(8,'(12E21.9)') pp(:), fff(pp), population0(pp), population1(pp), one-population0(pp)-population1(pp)
      else if (population1(pp).lt.1d-5) then
         write(9,'(12E21.9)') pp(:), fff(pp), population0(pp), population1(pp), one-population0(pp)-population1(pp)
         call a_bs(pp,as,bs)
         write(19,*) as,bs,fff(pp)
      else 
         write(10,'(12E21.9)') pp(:), fff(pp), population0(pp), population1(pp), one-population0(pp)-population1(pp)
      end if
   end if
   
case(2) !! using self consistent method


call a_bs_selfconsistent(pp,tolerancia,aij0,aij1,as,bs,sigma,pop0,pop1)

if (sigma.gt.0.0d0) then
   
   write(7,'(12E21.9)') pp(:), sigma, pop0, pop1, one-pop0-pop1
   if (pop0.lt.1.d-5) then
      write(8,'(12E21.9)') pp(:), sigma, pop0, pop1, one-pop0-pop1
      write(80,*) as,bs,sigma
   else if (pop1.lt.1d-5) then
      write(9,'(12E21.9)') pp(:), sigma, pop0, pop1, one-pop0-pop1
      write(90,*) as,bs,sigma
   else 
      write(10,'(12E21.9)') pp(:), sigma, pop0, pop1, one-pop0-pop1
      write(100,*) as,bs,sigma
   end if
   
end if

case default
   write(*,*)
   write(*,*) 'Select a method to obtain the optimal points'
   write(*,*) 
   stop
   
end select
   
end do
end do
end do
end do
end do
end do
end do
end do

end program maximize_penning
