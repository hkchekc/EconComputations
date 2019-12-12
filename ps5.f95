MODULE PS5PARA
    IMPLICIT NONE
    INTEGER, PARAMETER:: N=66, RETIRE_AGE=46 ! LIFE EXPECTANCY AND RETIREMENT
    REAL(KIND=8), PARAMETER:: A_MAX=2.0, A_MIN=-0.525, STEP=0.005
    INTEGER, PARAMETER:: NZ=2, NA=(A_MAX-A_MIN)/STEP
    REAL(KIND=8):: ALPHA=0.5, BETA=0.8
    REAL(KIND=8):: THETA=0.11, GAMMA=0.42, SIGMA=2.0, DELTA=0.06, POP_GROWTH= 0.011
    REAL(KIND=8), DIMENSION(2):: PROD_STATES = (/3.0, 0.5/)
    REAL(KIND=8), DIMENSION(2,2):: PROD_MARKOV = TRANSPOSE(RESHAPE((/0.9261,(1-0.9261),(1-0.9811),0.9811/),(/NZ,NZ/)))
END MODULE

PROGRAM PS5
    USE PS5PARA
    IMPLICIT NONE
    REAL(KIND=8):: WAGE, AGG_LABOR=1., AGG_CAP=1., BENEFIT, INTEREST, OPT_LABOR, NOM
    REAL(KIND=8), DIMENSION(NA):: A_GRID
    REAL(KIND=8), DIMENSION(N):: COHORT_POP
    REAL(KIND=8), DIMENSION(N,NZ):: AGE_EFF
    REAL(KIND=8), DIMENSION(NA,NZ,N):: VFUNC, LFUNC, STAT_DIST
    INTEGER, DIMENSION(NA,NZ,N):: PFUNC=0
    REAL(KIND=8):: ERROR=100. , CRIT_P=1e-3, LAST_CHOICE, CUR_CHOICE
    INTEGER:: SROWIDX, SCOLIDX, AGEIDX, I
    ! INIT ARRAYS
    A_GRID = (/(I*STEP, I=1,NA)/) + A_MIN - STEP
    COHORT_POP = (/((1.0/((1.0+POP_GROWTH)**I)), I=1,N)/)
    COHORT_POP(:) = COHORT_POP(:)/SUM(COHORT_POP)
    PRINT*, "cohort population",COHORT_POP
    AGE_EFF(:,:)=0.
    OPEN(UNIT=13, FILE='ef.txt', STATUS='OLD')
    DO SROWIDX=1,RETIRE_AGE-1
        READ(UNIT=13, FMT=*) AGE_EFF(SROWIDX,1)
    ENDDO
    CLOSE(UNIT=13)
    AGE_EFF(:,1) = AGE_EFF(:,1)*0.5
    AGE_EFF(:,2) = AGE_EFF(:,1)*3.0

    DO WHILE(ERROR>CRIT_P)
        ! UPDATE VARS
        WAGE = (1-ALPHA)*(AGG_CAP**ALPHA)*(AGG_LABOR**(-ALPHA))! FOC WRT LABOR
        INTEREST = ALPHA*(AGG_CAP**(ALPHA-1))*(AGG_LABOR**(1-ALPHA))-DELTA! FOC WRT CAPITAL
        BENEFIT = THETA*WAGE*AGG_LABOR/SUM(COHORT_POP(RETIRE_AGE:N)) ! TAX REV/RETIRED_POP
        PRINT*, "=================================================="
        PRINT*, "R, B,W:", INTEREST, BENEFIT, WAGE

        ! FIND OPTIMAL LABOR AND POPULATE LFUNC
        DO AGEIDX=1, RETIRE_AGE
            DO SROWIDX=1,NA
                DO SCOLIDX=1,NZ
                    CUR_CHOICE = A_GRID(PFUNC(SROWIDX,SCOLIDX,AGEIDX))
                    IF (AGEIDX/=1) THEN
                        LAST_CHOICE = (1+INTEREST)*A_GRID(PFUNC(SROWIDX,SCOLIDX,AGEIDX-1)) ! DISCOUNTED
                    ELSE
                        LAST_CHOICE = 0.
                    ENDIF
                    NOM = (GAMMA*(1-THETA)*PROD_STATES(SCOLIDX)*AGE_EFF(AGEIDX,SCOLIDX)*WAGE-(1-GAMMA)*(LAST_CHOICE-CUR_CHOICE))
                    OPT_LABOR = NOM/((1-THETA)*WAGE*PROD_STATES(SCOLIDX)*AGE_EFF(AGEIDX,SCOLIDX))
                    IF (OPT_LABOR > 1.0) THEN
                        OPT_LABOR = 1.0
                    ELSE IF (OPT_LABOR < 0.) THEN
                        OPT_LABOR = 0.0
                    ENDIF
                LFUNC(SROWIDX, SCOLIDX, AGEIDX) = OPT_LABOR
                ENDDO
            ENDDO
        ENDDO
        ! BACKWARD INDUCTION USING NEW VARS (BELLMAN)
        CALL BACKWARD_INDUCTION(INTEREST, BENEFIT, A_GRID, VFUNC, LFUNC, PFUNC)
        ! UPDATE STAT_DIST
        CALL FIND_STAT_DIST(A_GRID,COHORT_POP, PFUNC, STAT_DIST)
        ! OBSERVE ERROR AND UPDATE PRICE
        CALL UPDATE_PRICE(A_GRID, STAT_DIST, LFUNC, ERROR, AGG_CAP, AGG_LABOR)
    ENDDO

    WAGE = (1-ALPHA)*(AGG_CAP**ALPHA)*(AGG_LABOR**(-ALPHA))! FOC WRT LABOR
    INTEREST = ALPHA*(AGG_CAP**(ALPHA-1))*(AGG_LABOR**(1-ALPHA))-DELTA! FOC WRT CAPITAL
    BENEFIT = THETA*WAGE*AGG_LABOR/SUM(COHORT_POP(RETIRE_AGE:N)) ! TAX REV/RETIRED_POP

END PROGRAM PS5

! SUBROUTINES
SUBROUTINE BACKWARD_INDUCTION(INTEREST, BENEFIT, A_GRID, VFUNC, LFUNC, PFUNC)
    USE PS5PARA
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN):: INTEREST, BENEFIT
    REAL(KIND=8), DIMENSION(NA), INTENT(IN):: A_GRID
    REAL(KIND=8), DIMENSION(NA,NZ,N), INTENT(OUT):: VFUNC
    REAL(KIND=8), DIMENSION(NA,NZ,N), INTENT(INOUT):: LFUNC
    INTEGER, DIMENSION(NA,NZ,N), INTENT(INOUT):: PFUNC
    INTEGER:: AGEIDX, AGE, SROWIDX
    REAL(KIND=8):: CONSUM, UTIL, ERROR_VFI=100.
    REAL(KIND=8), DIMENSION(NA,NZ,N):: VFUNC_NEW

    ! FOR LAST AGE, ! SHOULD BE DETERMINISTIC, DO ONCE IS OK (NO CHOICE)
    DO SROWIDX=1,NA ! THIS IS THE CHOICE FROM LAST PERIOD/ CURRENT ASSET
            CONSUM = (1+INTEREST)*A_GRID(SROWIDX) + BENEFIT ! SPEND ALL IN LAST PERIOD
            PFUNC(:,:,N) = 1 ! ACTUALLY DOES NOT MATTER
            UTIL = CONSUM**((1-SIGMA)*GAMMA)/(1-SIGMA) ! NO NEXT PERIOD VALUE
            VFUNC_NEW(SROWIDX,:,N) = UTIL
    ENDDO

    DO AGEIDX = 2, N ! FOR AGE 20-65
        AGE = N- AGEIDX+1
        CALL BELLMAN(BENEFIT,INTEREST,AGE,A_GRID,VFUNC,LFUNC,VFUNC_NEW,PFUNC)
    ENDDO
    ERROR_VFI = MAXVAL(ABS(VFUNC_NEW-VFUNC))
    VFUNC = VFUNC_NEW
END SUBROUTINE

SUBROUTINE BELLMAN(BENEFIT,INTEREST,AGE,A_GRID,VFUNC,LFUNC,VFUNC_NEW,PFUNC)
    USE PS5PARA
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN):: BENEFIT, INTEREST
    INTEGER, INTENT(IN):: AGE
    REAL(KIND=8), DIMENSION(NA), INTENT(IN):: A_GRID
    REAL(KIND=8), DIMENSION(NA,NZ,N), INTENT(IN):: VFUNC, LFUNC
    REAL(KIND=8), DIMENSION(NA,NZ,N), INTENT(INOUT):: VFUNC_NEW
    INTEGER, DIMENSION(NA,NZ,N), INTENT(INOUT):: PFUNC
    REAL(KIND=8):: COND_MAX_UTIL, CONSUM, UTIL, WORKING, PROD, NEXTU
    INTEGER:: SROWIDX, SCOLIDX, CHOICEIDX

    IF (AGE< RETIRE_AGE) THEN ! WOKRING BELLMAN
        DO SROWIDX=1, NA
            DO SCOLIDX=1, NZ
                COND_MAX_UTIL = -1e12
                WORKING = LFUNC(SROWIDX, SCOLIDX, AGE)
                PROD = PROD_STATES(SCOLIDX)
                DO CHOICEIDX=1,NA ! LOOP OVER CHOICE OF ASSET PRIME
                    CONSUM = A_GRID(SROWIDX) + PROD*WORKING - INTEREST* A_GRID(CHOICEIDX)
                    IF (CONSUM > 0.) THEN
                        NEXTU = BETA*SUM(PROD_MARKOV(SCOLIDX,:)*VFUNC(CHOICEIDX,:,AGE+1))
                        UTIL = (CONSUM**GAMMA*(1-WORKING)**(1-GAMMA))**(1-SIGMA)/(1-SIGMA)+NEXTU
                        IF (UTIL>COND_MAX_UTIL) THEN
                            PFUNC(SROWIDX, SCOLIDX, AGE) = CHOICEIDX
                            COND_MAX_UTIL = UTIL
                        ENDIF
                    ENDIF
                ENDDO ! END LOOP CHOICE SPACE FOR ONE STATE
            ENDDO
        ENDDO
        VFUNC_NEW(SROWIDX, SCOLIDX, AGE) = COND_MAX_UTIL
    ELSE ! BELLMAN FOR RETIRED PEOPLE
        DO SROWIDX=1, NA
            COND_MAX_UTIL = -1e12
            DO CHOICEIDX=1, NA
                CONSUM = A_GRID(SROWIDX) + BENEFIT
                IF (CONSUM>0.) THEN
                    UTIL = CONSUM**((1-SIGMA)*GAMMA)/(1-SIGMA)+ BETA*(VFUNC(CHOICEIDX,1,AGE+1))
                    IF (UTIL>COND_MAX_UTIL) THEN
                        PFUNC(SROWIDX,:, AGE) = CHOICEIDX
                        COND_MAX_UTIL = UTIL
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        VFUNC_NEW(SROWIDX, :, AGE) = COND_MAX_UTIL ! PROD_STATE DO NOT AFFECT RETIRED PEOPLE
    ENDIF
END SUBROUTINE

SUBROUTINE FIND_STAT_DIST(A_GRID, COHORT_POP, PFUNC, STAT_DIST)
    USE PS5PARA
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(NA), INTENT(IN):: A_GRID
    REAL(KIND=8), DIMENSION(N), INTENT(IN):: COHORT_POP
    INTEGER, DIMENSION(NA,NZ,N), INTENT(IN):: PFUNC
    REAL(KIND=8), DIMENSION(NA,NZ,N), INTENT(OUT):: STAT_DIST
    INTEGER:: SROWIDX, SCOLIDX, AGEIDX, LAST_AGE, LAST_CHOICE
    REAL(KIND=8), DIMENSION(NA,NZ,N):: STAT_DIST_NEW
    REAL(KIND=8):: ERROR_STAT=100., CRIT_STAT=1e-2

    STAT_DIST(:,:,:) = 1.0/SIZE(STAT_DIST) ! START WITH UNIFORM
    ! EVERYONE HAVE ZERO ASSET AT STARTING AGE
    DO WHILE (ERROR_STAT>CRIT_STAT)
        STAT_DIST_NEW(:,:,:) = 0.
        STAT_DIST_NEW(:,1,1) = 0.2037* COHORT_POP(1)/SIZE(A_GRID) ! ALWAYS NO ASSET AT FIRST AGE, BUT 2 PROD STATES
        STAT_DIST_NEW(:,2,1) = 0.7963*COHORT_POP(1)/SIZE(A_GRID)
        DO AGEIDX=2, N
            DO SROWIDX=1, NA ! LAST PERIOD ASSET
                DO SCOLIDX=1, NZ !LAST PERIOD PRODUCTIVITY
                    LAST_CHOICE = PFUNC(SROWIDX, SCOLIDX, AGEIDX-1) ! TWO POSSIBLE
                    LAST_AGE = AGEIDX-1
                    STAT_DIST_NEW(LAST_CHOICE, 1, AGEIDX) = STAT_DIST_NEW(SROWIDX,1,LAST_AGE)*PROD_MARKOV(1,SCOLIDX)
                    STAT_DIST_NEW(LAST_CHOICE, 2, AGEIDX) =  STAT_DIST_NEW(SROWIDX,2,LAST_AGE)*PROD_MARKOV(2,SCOLIDX)
               ENDDO
            ENDDO
            STAT_DIST_NEW(:,:,AGEIDX) = STAT_DIST_NEW(:,:,AGEIDX)*COHORT_POP(AGEIDX)
        ENDDO
        ERROR_STAT = ABS(MAXVAL(STAT_DIST_NEW-STAT_DIST))
        STAT_DIST = STAT_DIST_NEW
    ENDDO
END SUBROUTINE

SUBROUTINE UPDATE_PRICE(A_GRID, STAT_DIST, LFUNC, ERROR, AGG_CAP, AGG_LABOR)
    ! COMPUTE ERROR AND UPDATE PRICE
    USE PS5PARA
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(NA), INTENT(IN):: A_GRID
    REAL(KIND=8), DIMENSION(NA,NZ,N), INTENT(IN):: STAT_DIST, LFUNC
    REAL(KIND=8), INTENT(OUT):: ERROR
    REAL(KIND=8), INTENT(INOUT):: AGG_CAP, AGG_LABOR
    INTEGER:: SROWIDX, SCOLIDX, AGEIDX
    REAL(KIND=8):: AGG_LABOR_NEW=0, AGG_CAP_NEW=0, PORTION, ERROR_LAB, ERROR_CAP

    DO AGEIDX=1, N ! CAL ERRORS
        DO SROWIDX=1,NA
            DO SCOLIDX=1,NZ
                PORTION = STAT_DIST(SROWIDX, SCOLIDX, AGEIDX)
                AGG_CAP_NEW = AGG_CAP_NEW+ A_GRID(SROWIDX)*PORTION
                IF (AGEIDX<RETIRE_AGE) THEN
                    AGG_LABOR_NEW = AGG_LABOR_NEW + LFUNC(SROWIDX, SCOLIDX, AGEIDX)*PORTION
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    ERROR_LAB = ABS(AGG_LABOR_NEW-AGG_LABOR)
    ERROR_CAP = ABS(AGG_CAP_NEW-AGG_CAP)
    PRINT*, "UPDATE_PRICE ERROR", ERROR_LAB, ERROR_CAP
    ERROR = MAX(ERROR_CAP, ERROR_LAB)
    AGG_LABOR = 0.05*AGG_LABOR_NEW + 0.95*AGG_LABOR
    AGG_CAP = 0.05*AGG_CAP_NEW + 0.95*AGG_CAP
    ! UPDATE INTEREST, WAGE, BENEFIT AT THE START OF NEXT PERIOD (IT IS THE SAME)
END SUBROUTINE


