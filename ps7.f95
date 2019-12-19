MODULE PS7PARA
    IMPLICIT NONE
    REAL(KIND=8), PARAMETER:: BETA=0.96, SIGMA=1, ALPHA=0.33, DELTA=0.06
    INTEGER, PARAMETER:: NK=40, NKBAR=10
    INTEGER:: I
    REAL(KIND=8), PARAMETER:: K_MIN=0.01, K_MAX=10., STEP=(10.-0.01)/40.
    REAL(KIND=8), DIMENSION(NK), PARAMETER:: K_GRID=(/(I*STEP, I=1, NK)/) +K_MIN - STEP
    REAL(KIND=8), PARAMETER:: KBAR_MIN=0.01, KBAR_MAX=10., BAR_STEP=4./9.
    REAL(KIND=8), DIMENSION(NK), PARAMETER:: KBAR_GRID=(/(I*BAR_STEP, I=1, NKBAR)/) +KBAR_MIN - BAR_STEP
    REAL(KIND=8):: A0=0.05, A1=0.9
END MODULE

MODULE PS7RES
    PS7PARA
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(NKBAR):: INTEREST, WAGE
    REAL(KIND=8), DIMENSION(NK):: K2BAR
    REAL(KIND=8), DIMENSION(NK):: UTIL=EXP(K_GRID)
    REAL(KIND=8), DIMENSION(NK):: B,C,D
    CONTAINS
        SUBROUTINE INIT_RW()
            DO I=1, NKBAR
                INTEREST(I)= ALPHA*(KBAR_GRID(I)**(ALPHA-1))-DELTA
                WAGE(I)= (1-ALPHA)*(KBAR_GRID(I)**ALPHA)
                K2BAR(I) = EXP(A0+A1*LOG(KBAR_GRID(I)))
            ENDDO
        END SUBROUTINE
END MODULE

PROGRAM PS7
    USE PS7PARA
    USE PS7RES
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(NKBAR):: INTEREST, WAGE
    INTEGER:: KBIDX, KIDX

    CALL SPLINE(K_GRID, K2BAR)
    DO KBIDX=1, NKBAR
        DO KIDX=1, NK
            CALL ISPLINE(KBIDX, KIDX)
        ENDDO
    ENDDO
        
    ENDDO

END PROGRAM PS7
! COPYING SOMEONE'S CODE
subroutine spline (X,Y)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i = 1, 2, ...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2+d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n >= 2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
    USE PS7PARA
    USE PS7RES
    implicit none
    double precision x(n), y(n)
    integer i, j, gap
    double precision h

    gap = n-1
    ! check input
    if ( n < 2 ) return
    if ( n < 3 ) then
      b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
    end if
    !
    ! step 1: preparation
    !
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
      d(i) = x(i+1) - x(i)
      b(i) = 2.0*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
    end do
    !
    ! step 2: end conditions
    !
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0.0
    c(n) = 0.0
    if(n /= 3) then
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if
    !
    ! step 3: forward elimination
    !
    do i = 2, n
      h = d(i-1)/b(i-1)
      b(i) = b(i) - h*d(i-1)
      c(i) = c(i) - h*c(i-1)
    end do
    !
    ! step 4: back substitution
    !
    c(n) = c(n)/b(n)
    do j = 1, gap
      i = n-j
      c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do
    !
    ! step 5: compute spline coefficients
    !
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
    do i = 1, gap
      b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
      d(i) = (c(i+1) - c(i))/d(i)
      c(i) = 3.*c(i)
    end do
    c(n) = 3.0*c(n)
    d(n) = d(n-1)
end subroutine spline

 function ispline(u)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points ! kbar and k in this context
! b, c, d = arrays of spline coefficients computed by spline ! in modules
! n       = the number of data points !
! output:
! ispline = interpolated value at point u
! Same Author as above
!=======================================================================
implicit none
double precision ispline
integer n
double precision  u
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
        j = k
    else
        i = k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u-x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline


SUBROUTINE MYPLINE(KBIDX, KIDX)
    USE PS7PARA
    USE PS7RES
    IMPLICIT NONE
    INTEGER, INTENT(IN):: KBIDX, KIDX
    REAL(KIND=8):: CONSUM

    CONSUM = (1+INTEREST(KBIDX))*K_GRID(KIDX) + WAGE(KIDX) -

END SUBROUTINE
