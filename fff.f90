MODULE fff_mod !! used for the function to minimize

use globalmod
USE nrtype

public

CONTAINS

!! the function fff(x) gives the ratio sigma_AI/sigma_PI for HeHe

  function fff(x)
    real(sp), dimension(:), intent(in) :: x
    real(sp) :: fff
    complex(sp) :: am1,a0,a1,bm1,b0,b1, aux
    real(sp) :: pop0, pop1 !! pop0 and pop1 are the populations of the molecular
                           !! channels with J=0 and 1, respectively
    !! x(1)=rho
    !! x(2)=eta
    !! x(3)=alpha0
    !! x(4)=alpha1
    !! x(5)=psi
    !! x(6)=chi
    !! x(7)=beta0
    !! x(8)=beta1

    a0=sin(x(1))*sin(x(2))*exp(ci*x(3))
    a1=sin(x(1))*sqrt(one-sin(x(2))**2.)*exp(ci*x(4))
    am1=a0*conjg(a0)+a1*conjg(a1)
    am1=sqrt(one-am1)


    b0=sin(x(5))*sin(x(6))*exp(ci*x(7))
    b1=sin(x(5))*sqrt(one-sin(x(6))**2.)*exp(ci*x(8))
    bm1=b0*conjg(b0)+b1*conjg(b1)
    bm1=sqrt(one-bm1)

    !population for J=0

    aux=a1*bm1-a0*b0+am1*b1
    pop0=third*aux*conjg(aux)

    !population for J=1

    aux=a1*b0-a0*b1
    pop1=aux*conjg(aux)

    aux=a1*bm1-am1*b1
    pop1=pop1+aux*conjg(aux)
    
    aux=a0*bm1-am1*b0
    pop1=pop1+aux*conjg(aux)

    pop1=pop1*half

    fff=-(aij0*pop0+aij1*pop1)!/(pij0*pop0+pij1*pop1)

 end function fff

subroutine a_bs(x,aa,bb) !! this subroutine gives the value of the coefficients a's and b's
   real(sp), dimension(:), intent(in) :: x
   complex(sp), intent(inout) :: aa(1:3), bb(1:3)
   complex(sp):: aux

   aa=(0.0d0,0.0d0)
   bb=(0.0d0,0.0d0)

   !! aa(1)-> a_-1
   !! aa(2)-> a_0
   !! aa(3)-> a_1
   !! bb(1)-> b_-1
   !! bb(2)-> b_0
   !! bb(3)-> b_1

    aa(2)=sin(x(1))*sin(x(2))*exp(ci*x(3))
    aa(3)=sin(x(1))*sqrt(one-sin(x(2))**2.)*exp(ci*x(4))
    aa(1)=aa(2)*conjg(aa(2))+aa(3)*conjg(aa(3))
    aa(1)=sqrt(one-aa(1))

    bb(2)=sin(x(5))*sin(x(6))*exp(ci*x(7))
    bb(3)=sin(x(5))*sqrt(one-sin(x(6))**2.)*exp(ci*x(8))
    bb(1)=bb(2)*conjg(bb(2))+bb(3)*conjg(bb(3))
    bb(1)=sqrt(one-bb(1))
   
 end subroutine a_bs

 function population0(x)
   real(sp), dimension(:), intent(in) :: x
   complex(sp) :: am1,a0,a1,bm1,b0,b1, aux
   real(sp) :: population0 !! population0 is the populations of the molecular
                          !! channels with J=0
   
    a0=sin(x(1))*sin(x(2))*exp(ci*x(3))
    a1=sin(x(1))*sqrt(one-sin(x(2))**2.)*exp(ci*x(4))
    am1=a0*conjg(a0)+a1*conjg(a1)
    am1=sqrt(one-am1)

    b0=sin(x(5))*sin(x(6))*exp(ci*x(7))
    b1=sin(x(5))*sqrt(one-sin(x(6))**2.)*exp(ci*x(8))
    bm1=b0*conjg(b0)+b1*conjg(b1)
    bm1=sqrt(one-bm1)

    !population for J=0

    aux=a1*bm1-a0*b0+am1*b1
    population0=third*aux*conjg(aux)
      
 end function population0

 function population1(x)
   real(sp), dimension(:), intent(in) :: x
   complex(sp) :: am1,a0,a1,bm1,b0,b1, aux
   real(sp) :: population0 !! population0 is the populations of the molecular
                           !! channels with J=0

   
    a0=sin(x(1))*sin(x(2))*exp(ci*x(3))
    a1=sin(x(1))*sqrt(one-sin(x(2))**2.)*exp(ci*x(4))
    am1=a0*conjg(a0)+a1*conjg(a1)
    am1=sqrt(one-am1)

    b0=sin(x(5))*sin(x(6))*exp(ci*x(7))
    b1=sin(x(5))*sqrt(one-sin(x(6))**2.)*exp(ci*x(8))
    bm1=b0*conjg(b0)+b1*conjg(b1)
    bm1=sqrt(one-bm1)

    !population for J=1

    aux=a1*b0-a0*b1
    population1=aux*conjg(aux)

    aux=a1*bm1-am1*b1
    population1=population1+aux*conjg(aux)
    
    aux=a0*bm1-am1*b0
    population1=population1+aux*conjg(aux)

    population1=population1*half

  end function population1

!! function to calculate the solutions by solving the self-consistent equations

  subroutine a_bs_selfconsistent(x,tolerance,ai0,ai1,aa,bb,sigma_total,pop0,pop1)
    real(sp), dimension(:), intent(in) :: x
    real(sp), intent(in) :: ai0, ai1 !! AI for channels J=0 and 1
    complex(sp), intent(inout) :: aa(1:3), bb(1:3)
    real(sp), intent(inout) :: sigma_total,pop0,pop1
    real(sp) :: tolerance
    complex(sp) :: aanew(1:3), bbnew(1:3)
    complex(sp) :: aux
    complex(sp) :: c11, c10, c1m1, c00 !! coupling coefficients
    complex(sp) :: dc11a1, dc10a1, dc00a1 !! derivative of c11 with respect to a1 and so on
    complex(sp) :: dc11b1, dc10b1, dc00b1 !! derivative of c11 with respect to b1 and so on
    complex(sp) :: dc11a0, dc1m1a0, dc00a0 
    complex(sp) :: dc11b0, dc1m1b0, dc00b0 
    complex(sp) :: dc10am1, dc1m1am1, dc00am1
    complex(sp) :: dc10bm1, dc1m1bm1, dc00bm1
    complex(sp) :: suma(1:3)
    complex(sp) :: sumb(1:3)
    real(sp) :: errora, errorb

    integer :: ij, jk, ik

    aanew=(0.0d0, 0.0d0)
    bbnew=(0.0d0, 0.0d0)
    sigma_total=-300.0d0 !! If sigma_total is smaller than 0 at the end of the subroutine, then it is not converged

    !! Set the initial coefficients

    call a_bs(x,aa,bb)
    
   !! aa(1)-> a_-1
   !! aa(2)-> a_0
   !! aa(3)-> a_1
   !! bb(1)-> b_-1
   !! bb(2)-> b_0
   !! bb(3)-> b_1

    do jk=1,10000

       !! initialize

       aanew=czero
       bbnew=czero

       !! calculation of numerator
       !! calculation of the coupling coefficients
       
       c11=aa(3)*bb(2)-aa(2)*bb(3)
       c11=c11*sqrt(half)
       
       c10=aa(3)*bb(1)-aa(1)*bb(3)
       c10=c10*sqrt(half)

       c1m1=aa(2)*bb(1)-aa(1)*bb(2)
       c1m1=c1m1*sqrt(half)
       
       c00=aa(3)*bb(1)-aa(2)*bb(2)+aa(1)*bb(3)
       c00=c00*sqrt(third)
      
       !! calculation of the derivatives of the coefficients
       
       dc11a1=bb(2)*sqrt(half)
       dc10a1=bb(1)*sqrt(half)
       dc00a1=bb(1)*sqrt(third)
       
       suma(3)=(conjg(c11)*dc11a1+conjg(c10)*dc10a1)*ai1
       suma(3)=suma(3)+conjg(c00)*dc00a1*ai0 !! denominator for a1
       
       dc11b1=-aa(2)*sqrt(half)
       dc10b1=-aa(1)*sqrt(half)
       dc00b1=aa(1)*sqrt(third)
       
       sumb(3)=(conjg(c11)*dc11b1+conjg(c10)*dc10b1)*ai1
       sumb(3)=sumb(3)+conjg(c00)*dc00b1*ai0 !! denominator for b1
       
       dc11a0=-bb(3)*sqrt(half)
       dc1m1a0=bb(1)*sqrt(half)
       dc00a0=-bb(2)*sqrt(third)
       
       suma(2)=(conjg(c11)*dc11a0+conjg(c1m1)*dc1m1a0)*ai1
       suma(2)=suma(2)+conjg(c00)*dc00a0*ai0
       
       dc11b0=aa(3)*sqrt(half)
       dc1m1b0=-aa(1)*sqrt(half)
       dc00b0=-aa(2)*sqrt(third)
       
       sumb(2)=(conjg(c11)*dc11b0+conjg(c1m1)*dc1m1b0)*ai1
       sumb(2)=sumb(2)+conjg(c00)*dc00b0*ai0
       
       dc10am1=-bb(3)*sqrt(half)
       dc1m1am1=-bb(2)*sqrt(half)
       dc00am1=bb(3)*sqrt(third)
       
       suma(1)=(conjg(c10)*dc10am1+conjg(c1m1)*dc1m1am1)*ai1
       suma(1)=suma(1)+conjg(c00)*dc00am1*ai0
       
       dc10bm1=aa(3)*sqrt(half)
       dc1m1bm1=aa(2)*sqrt(half)
       dc00bm1=aa(3)*sqrt(third)
       
       sumb(1)=(conjg(c10)*dc10bm1+conjg(c1m1)*dc1m1bm1)*ai1
       sumb(1)=sumb(1)+conjg(c00)*dc00bm1*ai0
       
       !! Calculation of denominator
       !! Population for J=0
       
       aux=aa(3)*bb(1)-aa(2)*bb(2)+aa(1)*bb(3)
       pop0=third*aux*conjg(aux)
       
       !! Population for J=1
       
       aux=aa(3)*bb(2)-aa(2)*bb(3)
       pop1=aux*conjg(aux)

       aux=aa(3)*bb(1)-aa(1)*bb(3)
       pop1=pop1+aux*conjg(aux)
       
       aux=aa(2)*bb(1)-aa(1)*bb(2)
       pop1=pop1+aux*conjg(aux)
       
       pop1=pop1*half
       
       !! Compute the total cross section

       sigma_tot=pop0*ai0+pop1*ai1

       !! Obtain new values
       do ij=1,3
          aanew(ij)=suma(ij)/sigma_tot          
          bbnew(ij)=sumb(ij)/sigma_tot
       end do
       
       !!normalize

       aanew(:)=conjg(aanew(:)/sqrt(dot_product(aanew,aanew)))
       bbnew(:)=conjg(bbnew(:)/sqrt(dot_product(bbnew,bbnew)))

       !! compute error
       
       errora=abs(dot_product(aanew,aa))-one
       errorb=abs(dot_product(bbnew,bb))-one      

       if (abs(errora).lt.tolerance.and.abs(errorb).lt.tolerance) then

          !! Calculate cross sections (TODO: this should be in a subroutine)

       !! Population for J=0
       
       aux=aa(3)*bb(1)-aa(2)*bb(2)+aa(1)*bb(3)
       pop0=third*aux*conjg(aux)
       
       !! Population for J=1
       
       aux=aa(3)*bb(2)-aa(2)*bb(3)
       pop1=aux*conjg(aux)

       aux=aa(3)*bb(1)-aa(1)*bb(3)
       pop1=pop1+aux*conjg(aux)
       
       aux=aa(2)*bb(1)-aa(1)*bb(2)
       pop1=pop1+aux*conjg(aux)
       
       pop1=pop1*half
          
       print*,
          print*, 'SUCCESS',jk
          sigma_total=pop1*ai1+pop0*ai0
          print*, 'Cross section', sigma_total
          print*, pop0, pop1, one-pop0-pop1
          exit
       end if
       
       !! set new values of the coefficients
       
       aa=aanew
       bb=bbnew
              
    end do

    if (errora.lt.tolerance.and.errorb.gt.tolerance) then 
       write(*,*) 'Not converged in a_bs_selfconsistent'
    end if

  end subroutine a_bs_selfconsistent

end MODULE fff_mod
