! global definition of parameters
! and general purpose routines
module globalmod

   implicit none

! global constants
  integer, parameter :: idp=SELECTED_real_KIND(12,200) ! idp=8	
  integer, parameter :: isp=4 ! single precision
  integer, parameter :: iint = SELECTED_INT_KIND(8) !standard integer
  integer, parameter :: ilng = SELECTED_INT_KIND(18) !largest possible integer
  integer(ilng), parameter :: imega = 1024*1024, igiga = imega*1024, itera = 1024*igiga

  integer, parameter :: nwin = 3 ! number of window operators 
  real(idp), parameter :: num_zero = 1d-15 !! numerical zero
  real(idp), parameter :: zero = 0.0_idp, one = 1.0_idp, two = 2.0_idp
  real(idp), parameter :: half = 0.5_idp, quart = 0.25_idp, three = 3.0_idp, four = 4.0_idp
  real(idp), parameter :: third = one/three
  real(idp), parameter :: Pi = 3.14159265358979323846_idp
  real(idp), parameter :: Pihalf = 1.5707963267948965600_idp
  real(idp), parameter :: FourPi = 4.0_idp * Pi

  real(idp), parameter :: kilo = 1.0e3_idp, mega = 1.0e6_idp, giga = 1.0e9_idp
  real(idp), parameter :: tera = 1.0e12_idp, peta = 1.0e15_idp 
 
  real(idp), parameter :: au2fs= 2.4189211d-2, au2eV=27.2113834_idp
  real(idp), parameter :: au2cm=2.194746313710d5, au2K=3.1577465d5
  real(idp), parameter :: au2Hz = 6.5796839207e15_idp, nm2cm = 1.0e7_idp
  real(idp), parameter :: au2MHz = 6.5796839207e9_idp !MHz in energy
  real(idp), parameter :: au2Vpcm=5.14225e+9_idp  ! atomic unit of electric field
  real(idp), parameter :: au2kVpcm=5.14225e+6_idp  ! atomic unit of electric field
  real(idp), parameter :: wpcm22au=7.123633018179838e-18_idp  ! atomic units for the whole factor multiplying the coupling with the laser
  real(idp), parameter :: eps_0=8.854187817e-12_idp !permittivity of vacuum
  real(idp), parameter :: mu_0=1.2566370614e-06_idp !permeability of vacuum
  real(idp), parameter :: c_light=2.99792458e+08_idp !speed of light
  real(idp), parameter :: debye2au=0.39343032_idp ! debye to atomic units
  real(idp), parameter :: a0=0.5291772108_idp !bohr radius 1/a0^3 polar2au
  real(idp), parameter :: Rb5s5p=12737.3499_idp  ! in cm^1
  real(idp), parameter :: RbD1=12578.9510_idp ! in cm^1
  real(idp), parameter :: RbD2=12816.5494_idp ! in cm^1
  real(idp), parameter :: Rbion=33690.94_idp ! in cm^1

  real(idp), parameter :: fac_fwhm = 4.0_idp

  complex(idp), parameter :: czero = (0.0_idp, 0.0_idp), &
       cone = (1.0_idp,0.0_idp),  ci= (0.0_idp, 1.0_idp), &
       cid = (1.0_idp, 1.0_idp)

! global output parameter
  character (LEN=4) :: runid 

! global string for error messages
  character(len=200) :: string

! global integers to keep track of required memory
  integer(ilng) :: imem, imemtemp

!! global variables for the associative and penning ionization

  real(kind=4) :: aij0, aij1, pij0, pij1

  public
  
contains

  subroutine allocerror (i,nn)
    integer :: i
    integer(ilng):: nn
    if (i > 0) then
       if (nn > igiga) then
          write (*,'(A,F8.2,A)') 'could not allocate', real(nn)/real(igiga), &
               ' GB of memory'
       else if (nn > imega) then
          write (*,'(A,F8.2,A)') 'could not allocate', real(nn)/real(imega), &
               ' MB of memory'
       else
          write (*,*) 'could not allocate memory', nn, ' B'
       end if
       if (imem > igiga) then
          write(*,'(A,F8.2,A)') 'we have already allocated', &
               real(imem)/real(igiga), ' GB of memory'
       else if (imem > imega) then
          write(*,'(A,F8.2,A)') 'we have already allocated', &
               real(imem)/real(imega), ' MB of memory'
       else
          write(*,*) 'we have already allocated', imem, ' B of memory'
       end if
       stop
    else
       imem = imem + nn
    end if
  end subroutine allocerror

  subroutine writemem (nn)
    integer(ilng):: nn
    if (nn > igiga) then
       write (*,'(F8.2,A)') real(nn)/real(igiga), ' GB of memory'
    else if (nn > imega) then
       write (*,'(F8.2,A)') real(nn)/real(imega), ' MB of memory'
    else 
       write (*,*) nn, ' B of memory'
    end if
  end subroutine writemem

  subroutine writemsg(string)
    character(len=*), intent(in) :: string
    write (6,*)  string
  end subroutine writemsg

  subroutine writeerror(string)
    character(len=*), intent(in) :: string
    write (6,*) 'error: ', string
    STOP !'program terminated by writeerror'
  end subroutine writeerror

  subroutine writetime(t)
    real :: t
    real :: sec
    integer:: days,hours,min
    call DHMS(t,days,hours,min,sec)
    write (*,'(A,I4,A,I2,A,I2,A,F5.2,A)') &
         ' Total Time  (CPU): ', days,' D ', hours,' H ', min,' M ', sec,' S'
    write(*,*)
  end subroutine writetime

  subroutine dhms(ss,d,h,m,s)
    integer :: d, h, m
    real :: ss, s
    d = int( ss / 86400.0 )
    ss = mod( ss, 86400.0 )
    h = int( ss / 3600.0 )
    ss = mod( ss, 3600.0 )
    m = int( ss / 60.0 )
    s = mod( ss, 60.0 )
  end subroutine dhms

  subroutine getfilenameending(cend,itime,istate)
    character(len=13) :: cend
    integer , optional :: itime, istate
    character (len=3) :: cstate
    character (len=2) :: ctim1
    character (len=3) :: ctim2
    character (len=4) :: ctim3
    character (len=5) :: ctim4
    character (len=6) :: ctim5
    character (len=7) :: ctim6
    character (len=8) :: ctim7
    character (len=9) :: ctim8
    character (len=10) :: ctim9
    if (present(istate)) then
       if (istate < 10) then;  write (cstate,'(A1,I1.1)') '_', istate
       else ;  write (cstate,'(A1,I2.2)') '_', istate
       end if
    end if
    if (present(itime)) then   
       if (itime < 10) then          
          write (ctim1,'(A1,BNI1.1)') '_',itime
          if (present(istate)) then;        cend = ctim1//cstate
          else ; cend = ctim1
          end if
       elseif (itime < 100) then
          write (ctim2,'(A1,I2.2)') '_',itime
          if (present(istate)) then;        cend = ctim2//cstate
          else ; cend = ctim2
          end if
       elseif (itime < 1000) then
          write (ctim3,'(A1,I3.3)') '_',itime
          if (present(istate)) then;        cend = ctim3//cstate
          else ; cend = ctim3
          end if
       elseif (itime < 10000) then
          write (ctim4,'(A1,I4.4)') '_',itime
          if (present(istate)) then;        cend = ctim4//cstate
          else ; cend = ctim4
          end if
       elseif (itime < 100000) then
          write (ctim5,'(A1,I5.5)') '_',itime
          if (present(istate)) then;        cend = ctim5//cstate
          else ; cend = ctim5
          end if
       else if (itime < 1000000) then
          write (ctim6,'(A1,I6.6)') '_',itime
          if (present(istate)) then;        cend = ctim6//cstate
          else ; cend = ctim6
          end if
       else if (itime < 10000000) then
          write (ctim7,'(A1,I7.7)') '_',itime
          if (present(istate)) then;        cend = ctim7//cstate
          else ; cend = ctim7
          end if
       else if (itime < 100000000) then
          write (ctim8,'(A1,I8.8)') '_',itime
          if (present(istate)) then;        cend = ctim8//cstate
          else ; cend = ctim8
          end if
       else if (itime < 1000000000) then
          write (ctim9,'(A1,I9.9)') '_',itime
          if (present(istate)) then;        cend = ctim9//cstate
          else ; cend = ctim9
          end if
       else
          write(*,*) 'this should not happen in getfilenamending!'
       end if
    else
       if (present(istate)) then;        cend = cstate
       else ; cend = ""
       end if
    end if
  end subroutine getfilenameending

  function check4zero (xx) result (yy)
    real(idp) :: xx, yy 
    real(idp), parameter :: eps = 0.1e-99_idp
    if (abs(xx) < eps) then
       yy = zero
    else
       yy = xx
    end if
  end function check4zero

  subroutine writegnufiles(nsurf,nstart,nend,nout,nend1,nout1,wig)
    ! write files to automize plotting with gnuplot
    ! in gnuplot type e.g.:  load 'r001/wig.gnu'
    implicit none
    integer :: nsurf, nstart, nend, nout, nend1, nout1
    logical, optional :: wig
    integer :: itime, isurf, iout
    character(13) :: cend
    character(1) :: csurf
    do isurf = 1 , nsurf
       write (csurf,'(I1.1)') isurf       
       open(10,file=runid//'/phi'//csurf//'r.gnu')
       open(11,file=runid//'/phi'//csurf//'k.gnu')
       if (present(wig).and.wig) open(12,file=runid//'/wig'//csurf//'.gnu')
       iout = nout1
       do itime = nstart, nend
          if (mod(itime,iout) == 0) then
             call getfilenameending(cend,itime) 
             write(10,'(A,A,A,A,A,A)') 'plot "', &
                  runid//'/phi'//csurf//'r'//trim(cend), '" w l'
             write(10,'(A)') 'pause -1'
             write(11,'(A,A,A,A,A,A)') 'plot "', &
                  runid//'/phi'//csurf//'k'//trim(cend), '" w l'
             write(11,'(A)') 'pause -1'
             if (present(wig).and.wig) then
                write(12,'(A,A,A,A,A,A,A)') 'splot "', &
                     runid//'/wig'//csurf//trim(cend), '" w l'
                write(12,'(A)') 'pause -1'             
             end if
          end if
          if (itime == nend1) then
             iout = nout
          end if
       end do
       close(10); close(11)
       if (present(wig).and.wig) close(12)
    end do
  end subroutine writegnufiles

  function powerof2(nn) result(res)
    implicit none
    integer :: nn, res
    res = 2**ceiling(log(real(nn,idp))/log(two))
  end function powerof2

  function phase(cc) result (xx)
    implicit none
    complex(idp) :: cc
    real(idp) :: xx
    if (real(cc,idp) > 1.0e-50_idp) then
       xx = atan(aimag(cc)/real(cc,idp))
    else
       xx = zero
    end if
  end function phase

  subroutine append_real(array,val) !! need to be checked
!! this subroutine appends a double precion value (val) to an array of 
!! double precision quantities
    implicit none
!!INPUT 
    real(idp), allocatable, intent(inout) :: array(:) !! also output
    real(idp), intent(in):: val
!!AUXILIAR
    real(idp), allocatable :: array_auxiliar(:)
    integer :: size_array
    
    if (.not.allocated(array)) then
       
       allocate(array(1))   
       array(1)=val
       
    else
       
       !!define the auxiliar array
       size_array=size(array)
       allocate(array_auxiliar(1:size_array))       
       array_auxiliar=array

       !!resize the output
       deallocate(array)
       allocate(array(1:size_array+1))
       
       !!restore the array
       array(1:size_array)=array_auxiliar
       
       !!add the new value 
       array(size_array+1)=val       

       !!deallocate the auxiliars
       deallocate(array_auxiliar)
    end if
       
  end subroutine append_real


  subroutine append_integer(array,val) !! need to be checked
!! this subroutine appends an integer value (val) to an array of 
!! double precision quantities
    implicit none
!!INPUT 
    integer, allocatable, intent(inout) :: array(:) !! also output
    integer, intent(in):: val
!!AUXILIAR
    integer, allocatable :: array_auxiliar(:)
    integer :: size_array
    
    if (.not.allocated(array)) then
       
       allocate(array(1))   
       array(1)=val
       
    else
       
       !!define the auxiliar array
       size_array=size(array)
       allocate(array_auxiliar(1:size_array))       
       array_auxiliar=array

       !!resize the output
       deallocate(array)
       allocate(array(1:size_array+1))
       
       !!restore the array
       array(1:size_array)=array_auxiliar
       
       !!add the new value 
       array(size_array+1)=val       

       !!deallocate the auxiliars
       deallocate(array_auxiliar)
    end if
       
  end subroutine append_integer

  !! resize an array to a given size (size_new)

  subroutine resize_real_array(array,size_new)
    !!INPUT
    real (idp),allocatable, intent(inout) :: array(:) !! also OUTPUT
    integer, intent(in) :: size_new
    !!OUTPUT
    !! the array
    !!AUXILIAR
    real (idp), allocatable :: array_auxiliar(:)
    integer :: ij,ik

    if (size(array).eq.size_new) return !! not need to resize
    
    allocate(array_auxiliar(size_new))

    !! copy the datas
    do ij=1,min(size_new,size(array))
       array_auxiliar(ij)=array(ij)
    end do

    deallocate(array)
    allocate(array(size_new))

    array=array_auxiliar !! copy the auxiliar array to the new one

    deallocate(array_auxiliar)
    
  end subroutine resize_real_array


  !! resize an array to a given size (size_new)

  subroutine resize_integer_array(array,size_new)
    !!INPUT
    integer,allocatable, intent(inout) :: array(:) !! also OUTPUT
    integer, intent(in) :: size_new
    !!OUTPUT
    !! the array
    !!AUXILIAR
    integer, allocatable :: array_auxiliar(:)
    integer :: ij,ik

    if (size(array).eq.size_new) return !! not need to resize
    
    allocate(array_auxiliar(size_new))

    !! copy the datas
    do ij=1,min(size_new,size(array))
       array_auxiliar(ij)=array(ij)
    end do

    deallocate(array)
    allocate(array(size_new))

    array=array_auxiliar  !! copy the auxiliar array to the new one

    deallocate(array_auxiliar)
    
  end subroutine resize_integer_array


end module globalmod
