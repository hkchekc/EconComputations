MODULE PS8PARA
    IMPLICIT NONE
    INTEGER:: I ! FOR ITER
    REAL(KIND=8):: BETA=0.99, ALPHA=0.36, DELTA=0.025, EFF=0.3271
    REAL(KIND=8), DIMENSION(2), PARAMETER:: PROD_SHOCK=(/1.01, 0.99/), EMP_STATE=(/1.,0./)

    REAL(KIND=8), PARAMETER:: K_MIN=0.01, K_MAX=15.0 ! DEFINE SMALL K GRIDS
    REAL(KIND=8), PARAMETER:: AGG_K_MIN=0.01, AGG_K_MAX=15.0

    INTEGER, PARAMETER:: NZ=SIZE(PROD_SHOCK), NY=SIZE(EMP_STATE), NK=100, N_AGGK=10
    REAL(KIND=8), PARAMETER:: KSTEP=(K_MAX-K_MIN)/FLOAT(NK), AGG_KSTEP=(AGG_K_MAX-AGG_K_MIN)/FLOAT(N_AGGK)
    REAL(KIND=8), DIMENSION(NK),PARAMETER:: K =(/(I*KSTEP, I=1,NK)/) + K_MIN - KSTEP
    REAL(KIND=8), DIMENSION(N_AGGK), PARAMETER:: AGG_K =(/(I*AGG_KSTEP, I=1,N_AGGK)/) + AGG_K_MIN - AGG_KSTEP
    REAL(KIND=8), DIMENSION(T):: AGG_SHOCK_VEC
    REAL(KIND=8), DIMENSION(T, N):: IND_SHOCK_VEC
    ! ASD
    REAL(KIND=8),PARAMETER:: U_G=0.04, U_B=0.1
    REAL(KIND=8), DIMENSION(2):: LAB_SUP=(/1.-U_G, 1.-U_B/)*EFF
    REAL(KIND=8), DIMENSION(2,2):: TRANS_MAT=( &
    (/0.850394, 0.115885, 0.024306, 0.009115/), &
    (/0.122917, 0.836111, 0.002083, 0.038889/), &
    (/0.583333, 0.031250, 0.291667, 0.093750/), &
    (/0.093750, 0.350000, 0.031250, 0.525000/))

    ! DROP: THE NUMBER OF PERIODS TO DROP INIT DEPENDENCE
    INTEGER, PARAMETER:: N=50, T=200, DROP=50

    REAL(KIND=8), DIMENSION(N_AGGK, NZ) = INTEREST, WAGE
    CONTAINS
        SUBROUTINE CAL_PRICE()
            INTEGER:: ZIDX, AKIDX
            REAL(KIND=8):: NORM_K
            ! UNRELATED TO INDIVIDUAL PROBLEMS
            DO AKIDX=1, N_AGGK
                DO ZIDX=1, NZ
                    NORM_K = AGG_K(AKIDX)/LAB_SUP(ZIDX)
                    INTEREST(AKIDX, ZIDX) = ALPHA*PROD_SHOCK(ZIDX)*NORM_K**(ALPHA-1)
                    WAGE(AKIDX, ZIDX) = (1-ALPHA)*PROD_SHOCK(ZIDX)*NORM_K**ALPHA
                ENDDO
            ENDDO
        END SUBROUTINE
    END MODULE

MODULE PS8RES
    USE PS8PARA
    USE PS8RES
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(NK, NY, N_AGGK, NZ):: VFUNC, VFUNC_NEW
    INTEGER, DIMENSION(NK, NY, N_AGGK, NZ):: PFUNC, PFUNC_NEW
    REAL(KIND=8), DIMENSION(2):: INTERCEPT=(/0.095, 0.085/), SLOPE=(/0.99, 0.99/),
    REAL(KIND):: R_SQ=0., ERROR=100.
    REAL(KIND=8), DIMENSION(T-DROP):: AGG_SHOCK_VEC_TRIM, K_VEC_TRIM
END MODULE

PROGRAM PS8
    USE PS8PARA
    USE PS8RES
    IMPLICIT NONE
    REAL(KIND=8):: ERROR=100., ERROR_VFI
    REAL(KIND=8), PARAMETER:: CRIT=1e-3, CRIT_VFI=1e-3

    CALL INIT_SHOCKS()
    DO WHILE (ERROR>CRIT)
        ERROR_VFI=100.
        DO WHILE (ERROR_VFI> CRIT_VFI) ! VFI
            CALL BELLMAN(ERROR_VFI)
        ENDDO
        CALL PSEUDO_PANEL()
        CALL CAL_ERRORS(ERROR)
        PRINT*, ERROR
    ENDDO
END PROGRAM

SUBROUTINE BELLMAN(ERROR_VFI)
    USE PS8PARA
    USE PS8RES
    IMPLICIT NONE
    REAL(KIND=8), INTENT(INOUT):: ERROR_VFI
    INTEGER:: KIDX, AKIDX, ZIDX, YIDX
    REAL(KIND=8):: KPR_EST

    ! DEFINE SOME INTERPOLATIONS

    DO AKIDX=1,N_AGGK
        DO ZIDX=1, NZ
            KPR_EST = EXP(INTERCEPT[ZIDX]+SLOPE[ZIDX]*LOG(K[AKIDX]))
            CALL FMINSEARCH()
            DO KIDX=1,NK
                DO YIDX=1,NY

                ENDDO
            ENDDO
        ENDDO
    ENDDO

    ERROR_VFI = MAXVAL(ABS(VFUNC_NEW-VFUNC))
    VFUNC = VFUNC_NEW
END SUBROUTINE

SUBROUTINE INIT_SHOCKS()
    USE PS8PARA
    USE PS8RES
    IMPLICIT NONE
    INTEGER:: TIDX, NIDX
    INTEGER:: SEED_Z, SEED_AGG
    REAL(KIND=8):: IND_PREV, Z_PREV
    REAL(KIND=8), DIMENSION(NZ):: TRANS_PROB

    AGG_SHOCK_VEC(1) = RAND ! STORE INDEX INSTEAD OF VALUE OF SHOCK HERE
    ! AGGREGATE SHOCKS
    DO TIDX=2,T
        Z_PREV = AGG_SHOCK_VEC(TIDX-1)
        TRANS_PROB = PI*RAND()
        !
        ! TODO: NORMALIZE TO 1,2
    ENDDO

    IND_SHOCK_VEC(:,1) = RAND
    ! INDIVIDUAL SHOCKS
    DO TIDX=2,T
        DO NIDX=1,N

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
    SHOCK_VEC_TRIM = SHOCK_VEC(DROP:T)
    K_VEC_TRIM = K_VEC(DROP:T)
END SUBROUTINE


SUBROUTINE CAL_ERRORS(ERROR)
    USE PS8PARA
    USE PS8RES
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT):: ERROR
    REAL(KIND=8):: R_SQ, ERROR_PARA(NZ*2)
    REAL(KIND=8), DIMENSION(2)::INTERCEPT_NEW, SLOPE_NEW
    INTEGER:: EIDX

    ! CALCUALTE NEW INTERCEPT AND SLOPE FOR FUNCTIONAL FORMS


    ! CALCULATE ERROR
    R_SQ =
    DO EIDX=1, SIZE(ERROR_PARA)
        ERROR_PARA(IDX) =
    ENDDO

    ERROR=MAXVAL(R_SQ, ERROR_PARA)
    INTERCEPT=INTERCEPT_NEW
    SLOPE=SLOPE_NEW
END SUBROUTINE

FUNCTION INTERPOLATION()
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN):: KPR
    REAL(KIND=8), INTENT(OUT):: INTE_VAL
END FUNCTION

FUNCTION CAL_KPR(KPR) RESULT(MINIMUM)
    IMPLICIT NONE
    REAL(KIND=8):: INTE_VAL
    REAL(KIND=8), INTENT(IN):: KPR
    REAL(KIND=8), INTENT(OUT):: MINIMUM

    CALL INTERPOLATION(KPR, INTE_VAL)
    MINIMUM = ABS(INTE_VAL-KPR_EST)
END FUNCTION

SUBROUTINE FMINSEARCH() ! WRAPPER
    USE TOOLBOX
    IMPLICIT NONE

END SUBROUTINE

