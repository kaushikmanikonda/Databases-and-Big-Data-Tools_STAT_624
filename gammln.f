      real function gammln(xx)

C     returns the value ln(Gamma(xx)) for xx > 0.  Full accuracy is
C     obtained for XX > 1.  For 0 < XX < 1, the reflection formula
C     can be used first.

      implicit none
      real xx
      integer j

      real*8 cof(6), stp, half, one, fpf, x, tmp, ser

C     Internal arithmetic will be done in double precision, a nicety
C     that you can omit if five-figure accuracy is good enough.

      data cof,stp /76.18009173D0,-86.50532033D0,24.01409822D0,
     $     -1.231739516D0,0.120858003D-2,-0.536382D-5,2.50662827465D0/
      data half,one,fpf /0.5D0, 1.0D0, 5.5D0/

      x = xx - one
      tmp = x + fpf
      tmp = (x + half)*log(tmp) - tmp
      ser = one
      do 11 j = 1, 6
         x = x + one
         ser = ser + cof(j)/x
 11   continue

      gammln = tmp + log(stp*ser)

      return
      end
