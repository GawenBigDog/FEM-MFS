module GAUSSM
    IMPLICIT NONE
    PUBLIC gauleg
!
    INTEGER,PARAMETER :: DBP = SELECTED_REAL_KIND(15,307)
    REAL(DBP) :: NEWV
    REAL(DBP) :: EPS,M_PI
    PARAMETER (EPS=3.0D-15)         ! EPS IS THE RELATIVE PRECISION
    PARAMETER (M_PI=3.141592654D0)  ! PI VALUE


    CONTAINS
!
    SUBROUTINE  GAULEG(ngp, xabsc, weig)
!{{{
!********************************************************************************
!* Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!********************************************************************************
      implicit none
      INTEGER  i, j, m
      REAL(dbp)  p1, p2, p3, pp, z, z1
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp), INTENT(OUT) :: xabsc(ngp), weig(ngp)
!
      m = (ngp + 1) / 2
!* Roots are symmetric in the interval - so only need to find half of them  */

      do i = 1, m       ! Loop over the desired roots */

         z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
!*   Starting with the above approximation to the ith root
!*          we enter the main loop of refinement by NEWTON'S method   */
100      p1 = 1.0d0
         p2 = 0.0d0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

         do j = 1, ngp
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
         enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
         pp = ngp*(z*p1-p2)/(z*z-1.0d0)
         z1 = z
         z = z1 - p1/pp             ! Newton's Method  */

         if (dabs(z-z1) .gt. EPS) GOTO  100

         xabsc(i) =  - z                        ! Roots will be bewteen -1.0 & 1.0 */
         xabsc(ngp+1-i) =  + z                  ! and symmetric about the origin  */
         weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
         weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

      end do     ! i loop
!
   RETURN
!}}}
   End subroutine GAULEG
!********************************************************************************
!*     Returns the SINGLE integral of the function (of ONE VARIABLE) "func"
!* between x1 and x2 by N-point Gauss-Legendre integration. The function
!* is evaluated exactly N times at interior points in the range of
!* integration.       */
!********************************************************************************
   recursive function qgauss(func, x1, x2, ngp) RESULT(intgrl)
!{{{
      implicit none
      REAL(dbp)  intgrl, x1, x2, func
      REAL(dbp)  xm, xl
      INTEGER j
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp) :: xabsc(ngp), weig(ngp)
!
      call gauleg(ngp, xabsc, weig)
!
      intgrl = 0.0d0
      xm = 0.5 * (x2 + x1)
      xl = 0.5 * (x2 - x1)
      do j = 1, ngp
        intgrl = intgrl + weig(j) * func( xm + xl*xabsc(j) )
      END do

      intgrl = intgrl * xl;    !Scale the answer to the range of integration  */
!}}}
   END function qgauss
!
end module GAUSSM
