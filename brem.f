C
      program brem
C
C     calculates bremsstrahlung straggling for electrons thru
C     rather thick (0.01 radiation lengths and greater) targets
C     or slabs.  Nothing fancy.  Uses formalism from Tsai,
C     Rev. Mod. Phys. 46, 815.  Particularly:
C     uses formulas 3.67 L_rad and 3.68 for L_rad'
C     uses formula 3.84 for the bremsstrahlung spectrum
C        NOTE this formula is not exact; gives about 2.5% error
C        in distribution for small energy losses
C     uses formula 4.3 to calculate the "b" parameter
C     uses formula 4.7 to actually calculate the straggling.
C
C     This formula is said to be accurate to within 1% (possibly
C     not accounting for the 2.5% above) for electron energies
C     greater than 0.2 times the incident energy and for absorbers
C     thinner than 0.05 radiation lengths.  For targets thicker
C     than 0.05 r.l., the straggling is iterated in steps of
C     0.05 r.l.
C
C     There are several approximations made which are only good
C     for Z >= 5, so don't expect good results for anything
C     lighter than carbon.
C
C     input:
C     card 1: the Z of the material and the A
C     card 2: the thickness of the material in cm
C     card 3: the density of the material in g/cm^3
C     card 4 ... the incident electron energy spectrum like:
C         0.05 0.001
C         0.15 0.02
C         0.25 0.5
C         ...
C         -10.0  -10.      (signaling a stop)
C
C     The bins are at the moment fixed with a 0.1 MeV step size
C     and centered as shown above.
C
C     output: file fort.4 receives some summary information
C     std. output receives a file which is in the format of
C     an energy spectrum suitable for direct inclusion into a
C     following instance of the program (without, of course,
C     the absorber data)
C
      implicit none
C
C     parameters for program brem
C
      integer EMAX              ! maximum incident energy
      integer BINMEV            ! bins per MeV
      parameter (EMAX=2002,BINMEV=10)
C
      integer i, j              ! loop variables
      real A                    ! nucleon number of the absorber
      real array1(EMAX*BINMEV)  ! array of edist
      real array2(EMAX*BINMEV)  ! array to hold result of convolution
      real ebins(EMAX*BINMEV)   ! array holding bin centers (MeV)
      real debins(EMAX*BINMEV)  ! array holding bin widths  (MeV)
      real lrad, lradprime      ! radiation logarithms
      real b                    ! parameter needed for dist func.
      real X0                   ! radiation length
      real Z                    ! atomic number of absorber
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

C to produce correct results with HPUX Fortran compiler reset array1
C (side benefit: also stifles complaint about unused variable j)

      do  j = 1, EMAX*BINMEV
         array1(j) = 0.0
      end do

 10   read (5,*) ept, distpt
      if (ept .ge. 0.0) then
         i = nint((ept + 0.5/BINMEV)*BINMEV)
         array1(i) = distpt
         go to 10
      endif
      maxbin = i
C
C     this is for the current setup: fixed bins of BINMEV size
C
      do i = 1, maxbin
         debins(i) = 1.0/BINMEV
         ebins(i)  = (i - 0.5)/BINMEV
      end do
C
C     compute Lrad and Lrad'
C
      lrad      = log(184.15/(Z**(1./3.)))
      lradprime = log(1194.0/(Z**(2./3.)))
      littlez   = (Z/137.0)**2
      fcoul     = 1.202*littlez - 1.0369*(littlez**2) +
     $     1.008*(littlez**3)/(1.0 + littlez)
C
C     compute X0
C
      X0 = 716.405*A/( (Z**2)*(lrad - fcoul) + Z*lradprime )
C
      open (unit=4,file='brem.log')
      write(4,*) 'Results for Radiation Lengths'
      write(4,*) 'Lrad = ', lrad, ', Lrad(prime) = ', lradprime
      write(4,*) 'Coulomb correction f = ', fcoul
      write(4,*) 'Radiation Length X0 = ', X0
C
C     calculate b parameter
C
      b = 4*(1.0 + (1.0/12.0)*(Z**2 + Z)/( (Z**2)*lrad + Z*lradprime))/3.
C
      write(4,*) 'Radiation parameter b = ', b
C
C     calculate thickness of absorber in rad lengths
C
      rth = thick*dens/X0
      write (4,*) 'Thickness = ', rth, ' radiation lengths'
C
C     make decision about looping or not
C
      if (rth .gt. 0.05) then
         looping = .true.
         nloops    = nint(rth/0.05)
         rem_thick = rth - nloops*0.05
      else
         looping = .false.
         rem_thick = rth
      endif
C
C     DEBUG
C
      write (4,*) 'looping = ', looping, ', nloops = ', nloops
      write (4,*) 'remaining thickness = ', rem_thick
C
C     do the work
C
      if (looping) then
         do loopi = 1, nloops
            write (4,*) 'Starting Loop Number ', loopi
            call conv_spec(array1,array2,ebins,debins,maxbin,b,0.05)
         end do
      endif
C
C     last (or perhaps only) step
C
      write (4,*) 'Starting Last Iteration'
      call conv_spec(array1,array2,ebins,debins,maxbin,b,rem_thick)
C
      do i = 1, maxbin
         write (*,'(5x, F8.3, 3x, E12.6)') ebins(i), array1(i)
      end do
C
      stop
      end
C
      subroutine conv_spec(edist_in,edist_out,ebins,debins,
     $     maxbin,b,rth)
C
      implicit none
C
      real edist_in(*)          ! incident energy distribution
      real edist_out(*)         ! energy dist to be returned
      real ebins(*)             ! bin centers in E for edists
      real debins(*)            ! bin widths for edists
      integer maxbin               ! maximum bin number that we need
      real delta                ! change in probability for each step
      real prob_res             ! remaining probability in bin (amount
C                                 that did not radiate
      real b                    ! radiation parameter
      real intens               ! function to return rad. prob. dist.
      real rth                  ! thickness in radiation-length units
      integer i, j              ! loop indices
C
C     clear out the temporary array
C
      do i = 1, maxbin
         edist_out(i) = 0.0
      end do
C
C     loop over incident energies, starting from the largest
C     stop at index 2 since we "can't" radiate to zero energy
C
      do i = maxbin, 2, -1
         prob_res = edist_in(i) ! to keep track of how much did not radiate
         if (edist_in(i) .gt. 1.0E-10) then
            do j = i - 1, 1, -1
               delta = edist_in(i)*intens(ebins(i),ebins(j),b,rth)*
     $              debins(j)
               edist_out(j) = edist_out(j) + delta
               prob_res = prob_res - delta ! this flux gets removed from i
            end do
         end if
C
C     add to output energy dist the amount of flux that did not radiate
C     in the input energy dist.  Note that this gets added outside
C     the if; thus if there is no radiation calculated due to low
C     probability content, then we just leave it there.
C
         edist_out(i) = edist_out(i) + prob_res
      end do
C
      do i = 1, maxbin
         edist_in(i) = edist_out(i)
      end do
C
      return
      end
C
      real function intens(ein,eout,b,rth)
C
      implicit none
C
      real ein                  ! incident energy
      real eout                 ! exit energy 
      real b                    ! radiation parameter
      real rth                  ! thickness of absorber
      real rhok
      real gammln
C
      real y
C
C
      if (eout .ge. ein) then
         write(*,*) 'INTENS -- Error - Cannot radiate to equal'
         write(*,*) '       -- or higher E!!'
         write(*,*) '       -- Exiting'
         stop
      endif
      y = (ein - eout)/ein
Cxx      write (16,*) 'ein = ', ein, ', eout = ', eout, ', y = ', y
      rhok = (1./(ein-eout))*(4/3. - 4*y/3 + y**2)
Cxx      write (16,*) 'rhok = ', rhok
      intens = (( (ein-eout)/ein )**(b*rth)) * rhok*rth /
     $     exp(gammln(1.0+b*rth))
C
      return
      end
