MODULE PS7PARA
    IMPLICIT NONE
    REAL(KIND=8), PARAMETER:: BETA=0.99, SIGMA=1, ALPHA=0.36
    INTEGER, PARAMETER:: NK=40, NAK=40
    INTEGER:: I
    REAL(KIND=8), PARAMETER:: K_MIN=0.01, K_MAX=0.5, STEP=(0.5-0.01)/FLOAT(NK-1)
    REAL(KIND=8), DIMENSION(NK), PARAMETER:: K_GRID=(/(I*STEP, I=1, NK)/) +K_MIN - STEP
    REAL(KIND=8), PARAMETER:: KBAR_MIN=0.15, KBAR_MAX=0.25, BAR_STEP=0.1/FLOAT(NAK-1)
    REAL(KIND=8), DIMENSION(NAK), PARAMETER:: AK_GRID=(/(I*BAR_STEP, I=1, NAK)/) +KBAR_MIN - BAR_STEP
    REAL(KIND=8):: A0=0.01, A1=0.99
    REAL(KIND=8), PARAMETER:: CRIT=1E-3
END MODULE

MODULE PS7RES
    USE PS7PARA
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(NAK):: INTEREST, WAGE
    INTEGER, DIMENSION(NK, NAK):: PFUNC
    REAL(KIND=8), DIMENSION(NK, NAK):: VFUNC,  VFUNC_NEW
    REAL(KIND=8), DIMENSION(NK):: KPR_EST
    REAL(KIND=8), DIMENSION(NK, NAK, NK):: UTIL, CONSUM
    CONTAINS
        SUBROUTINE INIT_RW()
            DO I=1, NAK
                INTEREST(I)= ALPHA*(AK_GRID(I)**(ALPHA-1.))
                WAGE(I)= (1.-ALPHA)*(AK_GRID(I)**ALPHA)
            ENDDO
        END SUBROUTINE
END MODULE

PROGRAM PS7
    USE PS7PARA
    USE PS7RES
    USE OMP_LIB
    IMPLICIT NONE
    INTEGER:: KIDX, AKIDX, KPIDX
    REAL(KIND=8):: ERROR_VFI

        CALL INIT_RW()

        ! FOR PARALLEL CALCULATION, REDUCE PRIVATE VARIABLES IN MAIN LOOP
        ! PRE CALCULATE THE KPR, CONSUM, UTILTIY
        DO AKIDX=1,NAK
            KPR_EST(AKIDX) = AK_GRID(AKIDX) ! NO MOTION
            PRINT*, KPR_EST(AKIDX)
        ENDDO

        UTIL = -1e12
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(AKIDX, KIDX, KPIDX)
        DO AKIDX=1,NAK
        DO KIDX=1,NK
        DO KPIDX=1,NK ! NEXT PERIOD K CHOICE
            CONSUM(KIDX, AKIDX, KPIDX) = (INTEREST(AKIDX))*K_GRID(KIDX)&
            +WAGE(AKIDX)-K_GRID(KPIDX)
            IF (CONSUM(KIDX, AKIDX, KPIDX)>0.) THEN
                UTIL(KIDX, AKIDX, KPIDX)= LOG(CONSUM(KIDX, AKIDX, KPIDX))
            ENDIF
        ENDDO
        ENDDO
        ENDDO
        !$OMP END PARALLEL DO

    ERROR_VFI = 100.
    DO AKIDX=1,NAK ! INIT VALUE
    DO KIDX=1,NK
        VFUNC(KIDX, AKIDX)= MAXVAL(UTIL(KIDX, AKIDX, :),1)
    ENDDO
    ENDDO
    DO WHILE (ERROR_VFI > CRIT)
        CALL VFI()
        ERROR_VFI = MAXVAL(ABS(VFUNC-VFUNC_NEW))
        VFUNC = VFUNC_NEW
    ENDDO

    CALL WRITE_ALL()
END PROGRAM PS7

SUBROUTINE VFI()
    USE PS7PARA
    USE PS7RES
    USE OMP_LIB
    IMPLICIT NONE
    INTEGER:: AKIDX, KIDX, KPIDX
    REAL(KIND=8), DIMENSION(1):: NEXTU, KPR, SMALL_KPR
    REAL(KIND=8), DIMENSION(NK, NAK, NK):: VFUNC_TMP


    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(AKIDX, KIDX, KPIDX)
    DO AKIDX=1,NAK
        KPR(1) = KPR_EST(AKIDX)
    DO KIDX=1,NK
    DO KPIDX=1,NK ! NEXT PERIOD K CHOICE
        SMALL_KPR(1) = K_GRID(KPIDX)
        CALL PWL_INTERP_2D(NK, NAK, K_GRID, AK_GRID, VFUNC, 1, SMALL_KPR, KPR, NEXTU)
        IF (NEXTU(1)>1E10 .OR. -1E10>NEXTU(1)) THEN
            PRINT*, SMALL_KPR(1),KPR(1), AK_GRID(NAK), K_GRID(NK),AK_GRID(1), K_GRID(1)
        ENDIF
        VFUNC_TMP(KIDX, AKIDX, KPIDX) = UTIL(KIDX, AKIDX, KPIDX)+BETA*NEXTU(1)
    ENDDO
        VFUNC_NEW(KIDX, AKIDX) = MAXVAL(VFUNC_TMP(KIDX, AKIDX, :),1)
        PFUNC(KIDX,AKIDX) = MAXLOC(VFUNC_TMP(KIDX, AKIDX, :),1)
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO
END SUBROUTINE

subroutine pwl_interp_2d ( nxd, nyd, xd, yd, zd, ni, xi, yi, zi )

!*****************************************************************************80
!
!! PWL_INTERP_2D: piecewise linear interpolant to data defined on a 2D grid.
!
!  Discussion:
!
!    Thanks to Adam Hirst for pointing out an error in the formula that
!    chooses the interpolation triangle, 04 February 2018.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NXD, NYD, the number of X and Y data values.
!
!    Input, real ( kind = 8 ) XD(NXD), YD(NYD), the sorted X and Y data.
!
!    Input, real ( kind = 8 ) ZD(NXD,NYD), the Z data.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), YI(NI), the coordinates of the
!    interpolation points.
!
!    Output, real ( kind = 8 ) ZI(NI), the value of the interpolant.
!
  implicit none

  integer ( kind = 4 ):: ni
  integer ( kind = 4 ):: nxd
  integer ( kind = 4 ):: nyd

  real ( kind = 8 ):: alpha1
  real ( kind = 8 ):: beta1
  real ( kind = 8 ):: det
  real ( kind = 8 ):: dxa
  real ( kind = 8 ):: dxb
  real ( kind = 8 ):: dxi
  real ( kind = 8 ):: dya
  real ( kind = 8 ):: dyb
  real ( kind = 8 ):: dyi
  real ( kind = 8 ):: gamma
  integer ( kind = 4 ):: i1
  integer ( kind = 4 ):: j
  integer ( kind = 4 ):: k
  real ( kind = 8 ):: r8_huge= 1.79769313486231571D+308
  integer ( kind = 4 ):: r8vec_bracket5
  real ( kind = 8 ):: xd(nxd)
  real ( kind = 8 ):: xi(ni)
  real ( kind = 8 ):: yd(nyd)
  real ( kind = 8 ):: yi(ni)
  real ( kind = 8 ):: zd(nxd,nyd)
  real ( kind = 8 ):: zi(ni)

  do k = 1, ni
!
!  For interpolation point (xi(k),yi(k)), find data intervals I and J so that:
!
!    xd(i) <= xi(k) <= xd(i+1),
!    yd(j) <= yi(k) <= yd(j+1).
!
!  But if the interpolation point is not within a data interval,
!  assign the dummy interpolant value zi(k) = infinity.
!
    i1 = r8vec_bracket5 ( nxd, xd, xi(k) )
    if ( i1 == -1 ) then
      zi(k) = r8_huge
      cycle
    end if

    j = r8vec_bracket5 ( nyd, yd, yi(k) )
    if ( j == -1 ) then
      zi(k) = r8_huge
      cycle
    end if
!
!  The rectangular cell is arbitrarily split into two triangles.
!  The linear interpolation formula depends on which triangle
!  contains the data point.
!
!    (I,J+1)--(I+1,J+1)
!      |\       |
!      | \      |
!      |  \     |
!      |   \    |
!      |    \   |
!      |     \  |
!    (I,J)---(I+1,J)
!
    if ( yi(k) < yd(j+1) &
      + ( yd(j) - yd(j+1) ) * ( xi(k) - xd(i1) ) / ( xd(i1+1) - xd(i1) ) ) then

      dxa = xd(i1+1) - xd(i1)
      dya = yd(j)   - yd(j)

      dxb = xd(i1)   - xd(i1)
      dyb = yd(j+1) - yd(j)

      dxi = xi(k)   - xd(i1)
      dyi = yi(k)   - yd(j)

      det = dxa * dyb - dya * dxb

      alpha1 = ( dxi * dyb - dyi * dxb ) / det
      beta1 =  ( dxa * dyi - dya * dxi ) / det
      gamma = 1.0D+00 - alpha1 - beta1

      zi(k) = alpha1 * zd(i1+1,j) + beta1 * zd(i1,j+1) + gamma * zd(i1,j)

    else

      dxa = xd(i1)   - xd(i1+1)
      dya = yd(j+1) - yd(j+1)

      dxb = xd(i1+1) - xd(i1+1)
      dyb = yd(j)   - yd(j+1)

      dxi = xi(k)   - xd(i1+1)
      dyi = yi(k)   - yd(j+1)

      det = dxa * dyb - dya * dxb

      alpha1 = ( dxi * dyb - dyi * dxb ) / det
      beta1 =  ( dxa * dyi - dya * dxi ) / det
      gamma = 1.0D+00 - alpha1 - beta1

      zi(k) = alpha1 * zd(i1,j+1) + beta1 * zd(i1+1,j) + gamma * zd(i1+1,j+1)

    end if
  end do
end subroutine

function r8vec_bracket5 ( nd, xd, xi )

!*****************************************************************************80
!
!! R8VEC_BRACKET5 brackets data between successive entries of a sorted R8VEC.
!
!  Discussion:
!
!    We assume XD is sorted.
!
!    If XI is contained in the interval [XD(1),XD(N)], then the returned
!    value B indicates that XI is contained in [ XD(B), XD(B+1) ].
!
!    If XI is not contained in the interval [XD(1),XD(N)], then B = -1.
!
!    This code implements a version of binary search which is perhaps more
!    understandable than the usual ones.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data values.
!
!    Input, real ( kind = 8 ) XD(N), the sorted data.
!
!    Input, real ( kind = 8 ) XD, the query value.
!
!    Output, integer ( kind = 4 ) R8VEC_BRACKET5, the bracket information.
!
  implicit none

  integer ( kind = 4 ) nd

  integer ( kind = 4 ) b
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  integer ( kind = 4 ) r8vec_bracket5
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi

  if ( xi < xd(1) .or. xd(nd) < xi ) then

    b = -1

  else

    l = 1
    r = nd

    do while ( l + 1 < r )
      m = ( l + r ) / 2
      if ( xi < xd(m) ) then
        r = m
      else
        l = m
      end if
    end do

    b = l

  end if

  r8vec_bracket5 = b

  return
end

SUBROUTINE WRITE_ALL()
    USE PS7RES
    USE PS7PARA
    IMPLICIT NONE
    INTEGER:: SROWIDX
    CHARACTER(LEN=130):: PATH="/Users/chek_choi/Downloads/fortran/"
    CHARACTER(LEN=150):: FILE_NAME
    REAL(KIND=8), DIMENSION(1):: NEXTU, SMALL_KPR

        FILE_NAME = TRIM(PATH)//"VFUNC"
        OPEN(UNIT=1, FILE=FILE_NAME, STATUS='REPLACE') ! START WITH THE TWO VALUE FUNCTIONS
        DO SROWIDX=1, NAK
            WRITE(UNIT=1,FMT=*) VFUNC(16,SROWIDX)
        ENDDO
        CLOSE(UNIT=1)

        FILE_NAME = TRIM(PATH)//"VFUNC_SS"
        SMALL_KPR(1) = 0.1995
        OPEN(UNIT=3, FILE=FILE_NAME, STATUS='REPLACE') ! START WITH THE TWO VALUE FUNCTIONS
        DO SROWIDX=1, NAK
            CALL PWL_INTERP_2D(NK, NAK, K_GRID, AK_GRID, VFUNC, 1, SMALL_KPR, AK_GRID(SROWIDX), NEXTU)
                IF (NEXTU(1)> 1E12) THEN
                    PRINT*, AK_GRID(SROWIDX), K_GRID(NK), AK_GRID(NAK)
                ENDIF
            WRITE(UNIT=3,FMT=*) NEXTU(1)
        ENDDO
        CLOSE(UNIT=3)

        FILE_NAME = TRIM(PATH)//"PFUNC"
        OPEN(UNIT=2, FILE=FILE_NAME, STATUS='REPLACE') ! ALSO SAVE POLICY FUNCTIONS
        DO SROWIDX=1, NK
            WRITE(UNIT=2,FMT=*) K_GRID(PFUNC(SROWIDX,NAK/2))
        ENDDO
        CLOSE(UNIT=2)

        FILE_NAME = TRIM(PATH)//"KGRID"
        OPEN(UNIT=4, FILE=FILE_NAME, STATUS='REPLACE') ! ALSO SAVE POLICY FUNCTIONS
        DO SROWIDX=1, NAK
            WRITE(UNIT=4,FMT=*)   K_GRID(SROWIDX)
        ENDDO
        CLOSE(UNIT=4)


END SUBROUTINE

