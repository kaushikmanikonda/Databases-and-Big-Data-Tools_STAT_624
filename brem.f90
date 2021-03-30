program brem
!
!     calculates bremsstrahlung straggling for electrons thru
!     rather thick (0.01 radiation lengths and greater) targets
!     or slabs.  Nothing fancy.  Uses formalism from Tsai,
!     Rev. Mod. Phys. 46, 815.  Particularly:
!     uses formulas 3.67 L_rad and 3.68 for L_rad'
!     uses formula 3.84 for the bremsstrahlung spectrum
!        NOTE this formula is not exact; gives about 2.5% error
!        in distribution for small energy losses
!     uses formula 4.3 to calculate the "b" parameter
!     uses formula 4.7 to actually calculate the straggling.
!
!     This formula is said to be accurate to within 1% (possibly
!     not accounting for the 2.5% above) for electron energies
!     greater than 0.2 times the incident energy and for absorbers
!     thinner than 0.05 radiation lengths.  For targets thicker
!     than 0.05 r.l., the straggling is iterated in steps of
!     0.05 r.l.
!
!     There are several approximations made which are only good
!     for Z >= 5, so don't expect good results for anything
!     lighter than carbon.
!
!     input:
!     card 1: the Z of the material and the A
!     card 2: the thickness of the material in cm
!     card 3: the density of the material in g/cm^3
!     card 4 ... the incident electron energy spectrum like:
!         0.05 0.001
!         0.15 0.02
!         0.25 0.5
!         ...
!         -10.0  -10.      (signaling a stop)
!
!     The bins are at the moment fixed with a 0.1 MeV step size
!     and centered as shown above.
!
!     output: file fort.4 receives some summary information
!     std. output receives a file which is in the format of
!     an energy spectrum suitable for direct inclusion into a
!     following instance of the program (without, of course,
!     the absorber data)
!
      implicit none
!
!     parameters for program brem
!
      integer EMAX              ! maximum incident energy
      integer BINMEV            ! bins per MeV
      parameter (EMAX=2002,BINMEV=10)
!
      integer i, j              ! loop variables
      real A                    ! nucleon number of the absorber
      real e_in, e_out          ! energies passed to distribution functions
      real array1(EMAX*BINMEV)  ! array of edist
      real array2(EMAX*BINMEV)  ! array to hold result of convolution
      real ebins(EMAX*BINMEV)   ! array holding bin centers (MeV)
      real debins(EMAX*BINMEV)  ! array holding bin widths  (MeV)
      real lrad, lradprime      ! radiation logarithms
      real b                    ! parameter needed for dist func.
      real x0                   ! radiation length
      real z                    ! atomic number of absorber
      real dens                 ! density of material in g/cm3
      real thick                ! thickness of absorber in cm
      real ept, distpt          ! energy and probab of input spectrum point
      real fcoul                ! coulomb correction in X0
      real littlez              ! used in calculation of f
      real rth                  ! thickness in units of X0
      logical looping           ! tells if we need to iterate spectrum
      integer nloops            ! number of times we need to iterate
      integer loopi             ! loop index for above
      real rem_thick            ! the remainder thickness after looping
                                ! so, the last iteration
      integer maxbin            ! max value necessary for loop index

      read (5,*) Z, A
      read (5,*) thick
      read (5,*) dens
!

      array1 = 0.0              ! initialize everything to zero for compilers
                                ! that don't do this.

 10   read (5,*) ept, distpt
      if (ept .ge. 0.0) then
         i = nint((ept + 0.5/BINMEV)*BINMEV)
         array1(i) = distpt
         go to 10
      endif
      maxbin = i
!
!     this is for the current setup: fixed bins of BINMEV size
!
      debins(:maxbin) = 1.0/BINMEV
      ebins( :maxbin) = (/ (I-0.5,I=1,maxbin) /)/BINMEV

!     compute Lrad and Lrad'

      lrad      = log(184.15/(Z**(1./3.)))
      lradprime = log(1194.0/(Z**(2./3.)))
      littlez   = (Z/137.0)**2
      fcoul     = 1.202*littlez - 1.0369*(littlez**2) + &
          1.008*(littlez**3)/(1.0 + littlez)

!     compute X0

      X0 = 716.405*A/( (Z**2)*(lrad - fcoul) + Z*lradprime )

      open (unit=4,file='brem.log')
      write(4,*) 'Results for Radiation Lengths'
      write(4,*) 'Lrad = ', lrad, ', Lrad(prime) = ', lradprime
      write(4,*) 'Coulomb correction f = ', fcoul
      write(4,*) 'Radiation Length X0 = ', X0
!
!     calculate b parameter
!
      b = 4*(1.0 + (1.0/12.0)*(Z**2 + Z)/( (Z**2)*lrad + Z*lradprime))/3.
!
      write(4,*) 'Radiation parameter b = ', b
!
!     calculate thickness of absorber in rad lengths
!
      rth = thick*dens/X0
      write (4,*) 'Thickness = ', rth, ' radiation lengths'
!
!     make decision about looping or not
!
      if (rth .gt. 0.05) then
         looping = .true.
         nloops    = nint(rth/0.05)
         rem_thick = rth - nloops*0.05
      else
         looping = .false.
         rem_thick = rth
      endif
!
!     DEBUG
!
      write (4,*) 'looping = ', looping, ', nloops = ', nloops
      write (4,*) 'remaining thickness = ', rem_thick
!
!     do the work
!
      if (looping) then
         do loopi = 1, nloops
            write (4,*) 'Starting Loop Number ', loopi
            call conv_spec(array1(:maxbin),ebins(:maxbin), &
                 debins(:maxbin),b,0.05)
         end do
      endif
!
!     last (or perhaps only) step
!
      write (4,*) 'Starting Last Iteration'
      call conv_spec(array1,ebins,debins,b,rem_thick)
!
      do i = 1, maxbin
         write (*,'(5x, F8.3, 3x, ES12.6)') ebins(i), array1(i)
      end do
!
      stop
      CONTAINS
!
      subroutine conv_spec(edist_in,ebins,debins,b,rth)
!
      implicit none
!
      real, dimension(:), intent(INOUT) :: edist_in  ! incident energy distribution
      real, dimension(:), intent(IN)    :: ebins     ! bin centers in E for edists
      real, dimension(:), intent(IN)    :: debins    ! bin widths for edists
      real,    intent(IN) :: b           ! radiation parameter
      real,    intent(IN) :: rth         ! thickness in radiation-length units


      real, dimension(SIZE(edist_in)) :: edist_out ! work array for energy dist to be returned
      real, dimension(SIZE(edist_in)) :: delta  ! change in probability for each step
      real prob_res             ! remaining probability in bin (amount
!                                 that did not radiate
      integer i, j              ! loop indices

!     clear out the temporary array

      edist_out = 0.0

!     loop over incident energies, starting from the largest
!     stop at index 2 since we "can't" radiate to zero energy

      do i = maxbin, 2, -1
         prob_res = edist_in(i) ! to keep track of how much did not radiate
         if (edist_in(i) .gt. 1.0E-10) then
            delta(:i-1) = edist_in(i) * &
                 intens(Ebins(i),Ebins(:i-1),b,rth) * dEbins(:i-1)
            delta(i:) = 0.0
            edist_out(:i-1) = edist_out(:i-1) + delta(:i-1)
            prob_res = prob_res - sum(delta)
         end if

!     add to output energy dist the amount of flux that did not radiate
!     in the input energy dist.  Note that this gets added outside
!     the if; thus if there is no radiation calculated due to low
!     probability content, then we just leave it there.

         edist_out(i) = edist_out(i) + prob_res
      end do
!
      edist_in = edist_out
!
      return
      end subroutine conv_spec
!
      function intens(ein,eout,b,rth)
!
      implicit none
!
      real, intent(IN) :: ein                  ! incident energy
      real, dimension(:), intent(IN) :: eout   ! exit energy
      real, intent(IN) :: b                    ! radiation parameter
      real, intent(IN) :: rth                  ! thickness of absorber
      real, dimension(size(eout))    :: intens ! array-valued function result

      real, dimension(size(eout)) :: rhok, &   ! variable which need
                                     y         ! dimension of eout
      real gammln
      integer i, j

!xx      if (eout .ge. ein) then
!xx         write(*,*) 'INTENS -- Error - Cannot radiate to equal'
!xx         write(*,*) '       -- or higher E!!'
!xx         write(*,*) '       -- Exiting'
!xx         stop
!xx      endif

      y = (ein - eout)/ein
      rhok = (1./(ein-eout))*(4/3. - 4*y/3 + y**2)
      intens = (( (ein-eout)/ein )**(b*rth)) * rhok*rth / &
           exp(gammln(1.0+b*rth))
!
      return
      end function intens

      end program brem

