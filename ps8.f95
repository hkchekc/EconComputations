MODULE PS8PARA
    IMPLICIT NONE
    INTEGER:: I ! FOR ITER
    REAL(KIND=8), PARAMETER:: BETA=0.99, ALPHA=0.36, DELTA=0.025, EFF=0.3271
    REAL(KIND=8), DIMENSION(2), PARAMETER:: PROD_SHOCK=(/1.01, 0.99/), EMP_STATE=(/1.,0./)
    ! DROP: THE NUMBER OF PERIODS TO DROP INIT DEPENDENCE
    INTEGER, PARAMETER:: N=50, T=200, DROP=50

    REAL(KIND=8), PARAMETER:: K_MIN=0.01, K_MAX=15.0 ! DEFINE SMALL K GRIDS
    REAL(KIND=8), PARAMETER:: AGG_K_MIN=0.01, AGG_K_MAX=15.0
    INTEGER, PARAMETER:: NZ=SIZE(PROD_SHOCK), NY=SIZE(EMP_STATE), NK=100, NAK=10
    REAL(KIND=8), PARAMETER:: KSTEP=(K_MAX-K_MIN)/FLOAT(NK), AGG_KSTEP=(AGG_K_MAX-AGG_K_MIN)/FLOAT(NAK)
    REAL(KIND=8), DIMENSION(NK),PARAMETER:: K_GRID =(/(I*KSTEP, I=1,NK)/) + K_MIN - KSTEP
    REAL(KIND=8), DIMENSION(NAK), PARAMETER:: AK_GRID =(/(I*AGG_KSTEP, I=1,NAK)/) + AGG_K_MIN - AGG_KSTEP
    INTEGER, DIMENSION(T):: AGG_SHOCK_VEC
    INTEGER, DIMENSION(N, T):: IND_SHOCK_VEC
    ! ASD
    REAL(KIND=8),PARAMETER:: U_G=0.04, U_B=0.1
    REAL(KIND=8), DIMENSION(2):: LAB_SUP=(/1.-U_G, 1.-U_B/)*EFF
    REAL(KIND=8), DIMENSION(4,4):: TRANS_MAT=&
    RESHAPE((/0.8507,0.1159,0.0243,0.0091,0.1229,0.8361,0.0021,&
    0.0389,0.5833,0.0313,0.2917,0.0938,0.0938,0.3500,0.0313,0.5250/), (/4,4/))

    REAL(KIND=8), DIMENSION(NAK, NZ):: INTEREST, WAGE
    CONTAINS
        SUBROUTINE CAL_PRICE()
            INTEGER:: ZIDX, AKIDX
            REAL(KIND=8):: NORM_K
            ! UNRELATED TO INDIVIDUAL PROBLEMS
            DO AKIDX=1, NAK
                DO ZIDX=1, NZ
                    NORM_K = AK_GRID(AKIDX)/LAB_SUP(ZIDX)
                    INTEREST(AKIDX, ZIDX) = ALPHA*PROD_SHOCK(ZIDX)*NORM_K**(ALPHA-1)
                    WAGE(AKIDX, ZIDX) = (1-ALPHA)*PROD_SHOCK(ZIDX)*NORM_K**ALPHA
                ENDDO
            ENDDO
        END SUBROUTINE
    END MODULE

MODULE PS8RES
    USE PS8PARA
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(NK, NY, NAK, NZ):: VFUNC, VFUNC_NEW
    INTEGER, DIMENSION(NK, NY, NAK, NZ):: PFUNC, PFUNC_NEW
    REAL(KIND=8), DIMENSION(2):: INTERCEPT=(/0.095, 0.085/), SLOPE=(/0.99, 0.99/)
    REAL(KIND=8):: R_SQ=0., ERROR=100.
    INTEGER, DIMENSION(T-DROP):: AGG_SHOCK_VEC_TRIM, K_VEC_TRIM
    REAL(KIND=8), DIMENSION(NAK, NZ):: KPR_EST
    REAL(KIND=8), DIMENSION(NK, NY, NAK, NZ, NK):: CONSUM, UTIL
    REAL(KIND=8), PARAMETER:: CRIT=1e-3
END MODULE

PROGRAM PS8
    USE PS8PARA
    USE PS8RES
    IMPLICIT NONE
    REAL(KIND=8):: ERROR_VFI

    CALL INIT_SHOCKS()
    DO WHILE (ERROR>CRIT)
        CALL VFI()
        CALL PSEUDO_PANEL()
        CALL CAL_ERRORS()
        PRINT*, ERROR
    ENDDO
END PROGRAM

SUBROUTINE VFI()
        USE PS8PARA
        USE PS8RES
        USE OMP_LIB
        IMPLICIT NONE
        INTEGER:: KIDX, AKIDX, ZIDX, YIDX, KPIDX
        REAL(KIND=8):: ERROR_VFI

        ! FOR PARALLEL CALCULATION, REDUCE PRIVATE VARIABLES IN MAIN LOOP
        ! PRE CALCULATE THE KPR, CONSUM, UTILTIY
        DO AKIDX=1,NAK !TODO: MOVVE THIS OUT OF VFI LOOP/ STATIC!
        DO ZIDX=1, NZ
            KPR_EST(AKIDX, ZIDX) = EXP(INTERCEPT(ZIDX)+SLOPE(ZIDX)*LOG(K_GRID(AKIDX)))
        ENDDO
        ENDDO

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(AKIDX, ZIDX, KIDX, YIDX)
        UTIL = -1e12
        DO AKIDX=1,NAK
        DO ZIDX=1, NZ
            DO KIDX=1,NK
            DO YIDX=1,NY
            DO KPIDX=1,NK ! NEXT PERIOD K CHOICE
                CONSUM(KIDX, YIDX, AKIDX, ZIDX, KPIDX) = (1+INTEREST(AKIDX, ZIDX)-DELTA)*K_GRID(KIDX)&
                +WAGE(AKIDX, ZIDX)*EMP_STATE(YIDX)-K_GRID(KPIDX)
                IF (CONSUM(KIDX, YIDX, AKIDX, ZIDX, KPIDX)>0.) THEN
                    UTIL(KIDX, YIDX, AKIDX, ZIDX, KPIDX)= LOG(CONSUM(KIDX, YIDX, AKIDX, ZIDX, KPIDX))
                ENDIF
            ENDDO
            ENDDO
            ENDDO
        ENDDO
        ENDDO
        !$OMP END PARALLEL DO

        ERROR_VFI=100.
        VFUNC_NEW = -1e12
        VFUNC_NEW = -1e12
        DO WHILE (ERROR_VFI> CRIT) ! VFI
            CALL BELLMAN()
            ERROR_VFI = MAXVAL(ABS(VFUNC_NEW-VFUNC))
            VFUNC = VFUNC_NEW
        ENDDO
END SUBROUTINE

SUBROUTINE BELLMAN()
    USE PS8PARA
    USE PS8RES
    USE OMP_LIB
    IMPLICIT NONE
    INTEGER:: AKIDX, ZIDX, KIDX,YIDX, ZZI,YYI, ROWIDX, NEXTROW
    REAL(KIND=8), DIMENSION(1):: KPR, SMALL_KPR !, NEXTUGE, NEXTUGN, NEXTUBE, NEXTUBN
    REAL(KIND=8), DIMENSION(1):: NEXTU
    REAL(KIND=8):: EXPU

    ! DEFINE SOME INTERPOLATIONS
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(AKIDX, ZIDX, KIDX, YIDX, KPR, SMALL_KPR)
    DO AKIDX=1,NAK
    DO ZIDX=1, NZ
        DO KIDX=1,NK
        DO YIDX=1,NY
            KPR(1) = KPR_EST(AKIDX, ZIDX)
            SMALL_KPR(1) = K_GRID(KIDX)
            EXPU = 0.
            ROWIDX = ZIDX + 2*(YIDX-1)
            DO ZZI=1, NZ
            DO YYI=1,NY
                NEXTROW = ZZI + 2*(YYI-1)
                CALL PWL_INTERP_2D(NK, NK, K_GRID, K_GRID, VFUNC(:,YYI, :, ZZI), 1, KPR, SMALL_KPR, NEXTU)
                EXPU = EXPU+TRANS_MAT(ROWIDX,NEXTROW)*NEXTU(1)
            ENDDO
            ENDDO
            VFUNC_NEW(KIDX, YIDX, AKIDX, ZIDX) = MAXVAL(UTIL(KIDX, YIDX, AKIDX, ZIDX, :))+ BETA*EXPU
            ! EXPU NOT BASED ON CHOICE, BECAUSE IT IS AUTOMATICALLY MAXIMIZED
            PFUNC(KIDX, YIDX, AKIDX, ZIDX) = MAXLOC(UTIL(KIDX, YIDX, AKIDX, ZIDX, :),1)
            ! EXPU IS OUT OF LOOP, INDEPENDENT OF CHOICE
        ENDDO
        ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE INIT_SHOCKS()
    USE PS8PARA
    USE PS8RES
    IMPLICIT NONE
    INTEGER:: TIDX, NIDX
    REAL(KIND=8):: IND_PREV, Z_PREV, RAND_TMP
    REAL(KIND=8), DIMENSION(NZ*NY):: TRANS_PROB
    INTEGER, PARAMETER:: RSEED=4200
    INTEGER:: NEXT_STATE, OTHER_STATE, CUR_SHOCK

    CALL SRAND(RSEED)

    ! FIRST ONE, DETERMINISTIC EVERY RUN (WILL DROP ANYWAY, IS OK)
    AGG_SHOCK_VEC(1) = 1 ! STORE INDEX INSTEAD OF VALUE OF SHOCK HERE
    ! AGGREGATE SHOCKS
    DO TIDX=2,T
        Z_PREV = AGG_SHOCK_VEC(TIDX-1)
        CALL RANDOM_NUMBER(RAND_TMP)
        NEXT_STATE = Z_PREV
        IF (RAND_TMP<1./8.) THEN ! RANDOM NEXT STATE (SYMMETRIC)
            OTHER_STATE = FINDLOC((/1,2/), (/1,2/)/=Z_PREV)
            NEXT_STATE = PROD_SHOCK(OTHER_STATE)
        ENDIF

    ENDDO


    IND_SHOCK_VEC(:,1) = 2 ! EVERYONE UNEMPLOYED STATE
    ! INDIVIDUAL SHOCKS
    DO TIDX=2,T
        CUR_SHOCK = AGG_SHOCK_VEC(TIDX)
        DO NIDX=1,N
            CALL RANDOM_NUMBER(RAND_TMP)
            TRANS_PROB = RAND_TMP
            IND_SHOCK_VEC(NIDX, TIDX) = COUNT(TRANS_MAT(IND_SHOCK_VEC(NIDX, TIDX-1),:)<TRANS_PROB)+1.
        ENDDO
    ENDDO
END SUBROUTINE

SUBROUTINE PSEUDO_PANEL()
    USE PS8PARA
    USE PS8RES
    IMPLICIT NONE
    REAL(KIND=8):: INIT_K=5.7163
    REAL(KIND=8), DIMENSION(T):: K_VEC ! DON'T OUTPUT
    INTEGER:: TIDX, ZIDX, NIDX

    ! NEED TRIMMED TO DROP INITIAL DEPENDENCE
    AGG_SHOCK_VEC_TRIM = AGG_SHOCK_VEC((DROP+1):T)
    K_VEC_TRIM = K_VEC(DROP+1:T)
END SUBROUTINE


SUBROUTINE CAL_ERRORS()
    USE PS8PARA
    USE PS8RES
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION((NZ*2)):: ERROR_PARA
    REAL(KIND=8), DIMENSION(2)::INTERCEPT_NEW, SLOPE_NEW
    INTEGER:: EIDX

    ! CALCUALTE NEW INTERCEPT AND SLOPE FOR FUNCTIONAL FORMS


    ! CALCULATE ERROR
    R_SQ =
    DO EIDX=1, SIZE(ERROR_PARA)
        ERROR_PARA(IDX) =
    ENDDO

    ERROR = MAXVAL(ERROR_PARA)
    ERROR=MAXVAL((/ERROR, R_SQ/))
    INTERCEPT=INTERCEPT_NEW
    SLOPE=SLOPE_NEW
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

  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nxd
  integer ( kind = 4 ) nyd

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) det
  real ( kind = 8 ) dxa
  real ( kind = 8 ) dxb
  real ( kind = 8 ) dxi
  real ( kind = 8 ) dya
  real ( kind = 8 ) dyb
  real ( kind = 8 ) dyi
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_huge
  integer ( kind = 4 ) r8vec_bracket5
  real ( kind = 8 ) xd(nxd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nyd)
  real ( kind = 8 ) yi(ni)
  real ( kind = 8 ) zd(nxd,nyd)
  real ( kind = 8 ) zi(ni)

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
    i = r8vec_bracket5 ( nxd, xd, xi(k) )
    if ( i == -1 ) then
      zi(k) = r8_huge ( )
      cycle
    end if

    j = r8vec_bracket5 ( nyd, yd, yi(k) )
    if ( j == -1 ) then
      zi(k) = r8_huge ( )
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
      + ( yd(j) - yd(j+1) ) * ( xi(k) - xd(i) ) / ( xd(i+1) - xd(i) ) ) then

      dxa = xd(i+1) - xd(i)
      dya = yd(j)   - yd(j)

      dxb = xd(i)   - xd(i)
      dyb = yd(j+1) - yd(j)

      dxi = xi(k)   - xd(i)
      dyi = yi(k)   - yd(j)

      det = dxa * dyb - dya * dxb

      alpha = ( dxi * dyb - dyi * dxb ) / det
      beta =  ( dxa * dyi - dya * dxi ) / det
      gamma = 1.0D+00 - alpha - beta

      zi(k) = alpha * zd(i+1,j) + beta * zd(i,j+1) + gamma * zd(i,j)

    else

      dxa = xd(i)   - xd(i+1)
      dya = yd(j+1) - yd(j+1)

      dxb = xd(i+1) - xd(i+1)
      dyb = yd(j)   - yd(j+1)

      dxi = xi(k)   - xd(i+1)
      dyi = yi(k)   - yd(j+1)

      det = dxa * dyb - dya * dxb

      alpha = ( dxi * dyb - dyi * dxb ) / det
      beta =  ( dxa * dyi - dya * dxi ) / det
      gamma = 1.0D+00 - alpha - beta

      zi(k) = alpha * zd(i,j+1) + beta * zd(i+1,j) + gamma * zd(i+1,j+1)

    end if

  end do

  return
end subroutine
