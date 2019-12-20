MODULE PS4PARA
    IMPLICIT NONE
    INTEGER:: I ! FOR ITER
    REAL(KIND=8), PARAMETER:: ALPHA=1.5, BETA=0.9932
    REAL(KIND=8), DIMENSION(2), PARAMETER:: STATES=(/1.0, 0.5/)
    INTEGER, PARAMETER:: NZ=SIZE(STATES), NA=700
    REAL(KIND=8), DIMENSION(NZ, NZ), PARAMETER:: PI=TRANSPOSE(RESHAPE((/0.97, 0.03,0.5,0.5/),(/2,2/)))
    REAL(KIND=8), PARAMETER:: A_MIN=-2., A_MAX=5., STEPS=(A_MAX-A_MIN)/(FLOAT(NA)-1.)
    REAL(KIND=8), DIMENSION(NA), PARAMETER:: A_GRID=(/(I*STEPS, I=1, NA)/) + A_MIN - STEPS
    LOGICAL:: COMPLETE=.FALSE. ! COMPLETE MARKET IF 1 OR INCOMPLETE MARKET IF 0
END MODULE

MODULE PS4RES
    USE PS4PARA
    IMPLICIT NONE
    REAL(KIND=8):: HIGH_Q=1., LOW_Q=0.9932, Q
    REAL(KIND=8), DIMENSION(NA,NZ):: VFUNC=0., VFUNC_NEW=0., CONSUM_ARR, LAMBDA
    REAL(KIND=8), DIMENSION(NA*NZ):: STAT_DIST
    REAL(KIND=8), DIMENSION(NA*NZ,NA*NZ):: A_CHANGE_MAT
    REAL(KIND=8), DIMENSION(700):: LORENZ
    INTEGER, DIMENSION(NA,NZ):: PFUNC=0
END MODULE

PROGRAM PS4
    USE PS4PARA
    USE PS4RES
    IMPLICIT NONE
    INTEGER:: SROWIDX
    REAL(KIND=8):: ERROR_Q=100., ERROR_VFI, CRIT=1e-3
    REAL(KIND=8), DIMENSION(NA,NZ):: VFUNC_INCOMPLETE=0.

    DO WHILE (ERROR_Q> CRIT)
        PRINT*, "Q", Q
        Q = (HIGH_Q+LOW_Q)/2.
        ERROR_VFI=100.
        IF (COMPLETE) THEN
                ! MUST RUN THE INCOMPLETE VERSION PRIOR TO
                OPEN(UNIT=13, FILE='VFUNC', STATUS='OLD')
                DO SROWIDX=1,NA
                    READ(UNIT=13, FMT=*) VFUNC_INCOMPLETE(SROWIDX,:)
                ENDDO
                CLOSE(UNIT=13)
            DO WHILE (ERROR_VFI> CRIT)
                CALL BELLMAN_COMPLETE(ERROR_VFI, VFUNC_INCOMPLETE)
            ENDDO
        ELSE
            DO WHILE (ERROR_VFI> CRIT)
                CALL BELLMAN(ERROR_VFI)
            ENDDO
        ENDIF
        CALL POP_A_CHANGE_MAT()
        CALL FIND_STAT_DIST()
        CALL COMPUTE_ERROR(ERROR_Q)
    ENDDO
    PRINT*, "Q", Q
    CALL FIND_LORENZ()
    IF (COMPLETE) THEN
        CALL WELFARE(VFUNC_INCOMPLETE)
    ENDIF
    CALL WRITE_ALL()
END PROGRAM PS4

SUBROUTINE BELLMAN(ERROR_VFI)
    USE PS4PARA
    USE PS4RES
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT):: ERROR_VFI
    REAL(KIND=8):: CONSUM, UTIL, COND_UTIL, CU, NU
    REAL(KIND=8), DIMENSION(NA,NZ):: ABS_DIFF
    INTEGER:: RIDX, CIDX, CHOICE

    DO RIDX=1, NA
        DO CIDX=1, NZ
            COND_UTIL=-1e10
            DO CHOICE=1, NA
                CONSUM = STATES(CIDX) + A_GRID(RIDX) - Q*A_GRID(CHOICE)
                IF (CONSUM>0.) THEN
                    CU = CONSUM**(1-ALPHA)/(1-ALPHA) ! CURRENT UTIL
                    NU = SUM(PI(CIDX,:)*VFUNC(CHOICE,:))! NEXT EXPECTED UTIL
                    UTIL = CU + BETA*NU
                    IF (UTIL>COND_UTIL) THEN
                        VFUNC_NEW(RIDX, CIDX)=UTIL
                        PFUNC(RIDX, CIDX)=CHOICE
                        COND_UTIL=UTIL
                        CONSUM_ARR(RIDX, CIDX) = CONSUM
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    ABS_DIFF = ABS(VFUNC_NEW-VFUNC)
    ERROR_VFI = MAXVAL(ABS_DIFF)
    VFUNC = VFUNC_NEW
    
END SUBROUTINE

SUBROUTINE BELLMAN_COMPLETE(ERROR_VFI, VFUNC_INCOMPLETE)
    USE PS4PARA
    USE PS4RES
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT):: ERROR_VFI
    REAL(KIND=8):: CONSUM, UTIL, COND_UTIL, CU, NOM, DENOM
    REAL(KIND=8), DIMENSION(NA,NZ):: ABS_DIFF
    REAL(KIND=8), DIMENSION(NA,NZ), INTENT(IN):: VFUNC_INCOMPLETE
    INTEGER:: RIDX, CIDX, CHOICE

    DO RIDX=1, NA
        DO CIDX=1, NZ
            COND_UTIL=-1e10
            DO CHOICE=1, NA
                CONSUM = STATES(CIDX) + A_GRID(RIDX) - Q*A_GRID(CHOICE)
                IF (CONSUM>0.) THEN
                    CU = CONSUM**(1-ALPHA)/(1-ALPHA) ! CURRENT UTIL
                    UTIL = CU ! COMPLETE MARKET NEXT UTIL= CURRENT UTIL
                    IF (UTIL>COND_UTIL) THEN
                        VFUNC_NEW(RIDX, CIDX)=UTIL
                        PFUNC(RIDX, CIDX)=CHOICE
                        COND_UTIL=UTIL
                        CONSUM_ARR(RIDX, CIDX) = CONSUM
                    ENDIF
                ENDIF
            ENDDO
        NOM = VFUNC_NEW(RIDX, CIDX) + 1/((1-ALPHA)*(1-BETA))
        DENOM = VFUNC_INCOMPLETE(RIDX, CIDX) + 1/((1-ALPHA)*(1-BETA))
        LAMBDA(RIDX, CIDX) = (NOM/DENOM)**(1/(1-ALPHA))-1
        ENDDO
    ENDDO
    ABS_DIFF = ABS(VFUNC_NEW-VFUNC)
    ERROR_VFI = MAXVAL(ABS_DIFF)
    VFUNC = VFUNC_NEW

END SUBROUTINE

SUBROUTINE POP_A_CHANGE_MAT()
    USE PS4PARA
    USE PS4RES
    IMPLICIT NONE
    INTEGER:: RIDX, CIDX, CUR_IDX, NEXT_IDX, CHOICE, NEXT_Z

    A_CHANGE_MAT=0.
    DO RIDX=1, NA
        DO CIDX=1, NZ
            CUR_IDX = 2*RIDX -2 + CIDX
            CHOICE=PFUNC(RIDX, CIDX)
            DO NEXT_Z=1, NZ
                NEXT_IDX= 2*CHOICE -2 + NEXT_Z
                A_CHANGE_MAT(CUR_IDX, NEXT_IDX) = A_CHANGE_MAT(CUR_IDX, NEXT_IDX)+ PI(CIDX, NEXT_Z)
            ENDDO
        ENDDO
    ENDDO
END SUBROUTINE

SUBROUTINE FIND_STAT_DIST()
    USE PS4PARA
    USE PS4RES
    IMPLICIT NONE
    REAL(KIND=8):: ERROR_STAT, CRIT_STAT=1e-3
    REAL(KIND=8), DIMENSION(NA*NZ) :: STAT_DIST_NEW, ABS_DIFF

    STAT_DIST= 1/(FLOAT(NA)*2.)
    ERROR_STAT=100.
    DO WHILE (ERROR_STAT>CRIT_STAT)
        STAT_DIST_NEW = MATMUL(STAT_DIST, A_CHANGE_MAT)
        ABS_DIFF = ABS(STAT_DIST_NEW-STAT_DIST)
        ERROR_STAT = MAXVAL(ABS_DIFF)
        STAT_DIST = STAT_DIST_NEW
    ENDDO
END SUBROUTINE

SUBROUTINE COMPUTE_ERROR(ERROR_Q)
    USE PS4PARA
    USE PS4RES
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT):: ERROR_Q
    INTEGER:: RIDX, CIDX
    REAL(KIND=8):: NET_ASSET! MARKET CLEARING

    NET_ASSET = 0.
    DO RIDX=1, NA
        DO CIDX=1, NZ
            NET_ASSET = NET_ASSET + A_GRID(PFUNC(RIDX, CIDX))*STAT_DIST(2*RIDX -2 + CIDX)
        ENDDO
    ENDDO

    IF (NET_ASSET>0.) THEN
        HIGH_Q=Q
    ELSE
        LOW_Q = Q
    ENDIF
    ERROR_Q = ABS(HIGH_Q-LOW_Q)
END SUBROUTINE

SUBROUTINE FIND_LORENZ()
    USE PS4PARA
    USE PS4RES
    IMPLICIT NONE
    INTEGER:: RIDX
    REAL(KIND=8):: GINI
    REAL(KIND=8), DIMENSION(NA):: DEGREE

    DEGREE=(/(I*(1./FLOAT(NA)), I=1,NA)/)

    GINI= 0.
    LORENZ(1) = STAT_DIST(1)
      DO RIDX = 2, 700
          LORENZ(RIDX)= LORENZ(RIDX-1) + STAT_DIST(RIDX*2-2+1) + STAT_DIST(RIDX*2-2+2)
          GINI = GINI + ABS(DEGREE(RIDX)-LORENZ(RIDX))
      ENDDO
    PRINT*, "degree", degree(na), degree(na-1)
    PRINT*,"GINI", GINI/FLOAT(NA)
END SUBROUTINE

SUBROUTINE WRITE_ALL()
    USE PS4PARA
    USE PS4RES
    IMPLICIT NONE
    INTEGER:: SROWIDX
    INTEGER, DIMENSION(2):: SIDX
    CHARACTER(LEN=130):: PATH="./"
    CHARACTER(LEN=150):: FILE_NAME

        FILE_NAME = TRIM(PATH)//"VFUNC"
        OPEN(UNIT=1, FILE=FILE_NAME, STATUS='REPLACE') ! START WITH THE TWO VALUE FUNCTIONS
        DO SROWIDX=1, NA
            WRITE(UNIT=1,FMT=*) VFUNC(SROWIDX,:)
        ENDDO
        CLOSE(UNIT=1)

        FILE_NAME = TRIM(PATH)//"PFUNC"
        OPEN(UNIT=2, FILE=FILE_NAME, STATUS='REPLACE') ! ALSO SAVE POLICY FUNCTIONS
        DO SROWIDX=1, NA
            WRITE(UNIT=2,FMT=*) PFUNC(SROWIDX,:)
        ENDDO
        CLOSE(UNIT=2)

        FILE_NAME = TRIM(PATH)//"STATDIST"
        OPEN(UNIT=3, FILE=FILE_NAME, STATUS='REPLACE') ! ALSO SAVE POLICY FUNCTIONS
        DO SROWIDX=1, NA
            SIDX = (/SROWIDX*2-1, SROWIDX*2/)
            WRITE(UNIT=3,FMT=*) STAT_DIST(SIDX)
        ENDDO
        CLOSE(UNIT=3)

        FILE_NAME = TRIM(PATH)//"AGRID"
        OPEN(UNIT=4, FILE=FILE_NAME, STATUS='REPLACE') ! FOR HAVING THE X-AXIS OF PLOT
        DO SROWIDX=1, NA
            WRITE(UNIT=4,FMT=*) A_GRID(SROWIDX)
        ENDDO
        CLOSE(UNIT=4)

        FILE_NAME = TRIM(PATH)//"LORENZ"
        OPEN(UNIT=5, FILE=FILE_NAME, STATUS='REPLACE') ! FOR HAVING THE X-AXIS OF PLOT
        DO SROWIDX=1, 100
            WRITE(UNIT=5,FMT=*) LORENZ(SROWIDX)
        ENDDO
        CLOSE(UNIT=5)

        IF (COMPLETE) THEN
            FILE_NAME = TRIM(PATH)//"LAMBDA"
            OPEN(UNIT=6, FILE=FILE_NAME, STATUS='REPLACE') ! FOR HAVING THE X-AXIS OF PLOT
            DO SROWIDX=1, NA
                WRITE(UNIT=6,FMT=*) LAMBDA(SROWIDX,:)
            ENDDO
            CLOSE(UNIT=6)
        ENDIF
END SUBROUTINE

SUBROUTINE WELFARE(VFUNC_INCOMPLETE)
    USE PS4PARA
    USE PS4RES
    IMPLICIT NONE
    REAL(KIND=8):: WG=0., WINC=0., WFB=0., VOTE=0.
    INTEGER:: RIDX, CIDX, DIS_LOC
    REAL(KIND=8), DIMENSION(NA,NZ), INTENT(IN):: VFUNC_INCOMPLETE
    PRINT*, SUM(STAT_DIST), "STAT"
    DO RIDX=1,NA
        DO CIDX=1, NZ
            DIS_LOC = RIDX*2-2 + CIDX
            WG = WG+ LAMBDA(RIDX, CIDX)* STAT_DIST(DIS_LOC)
            WINC = WINC +VFUNC_INCOMPLETE(RIDX, CIDX)*STAT_DIST(DIS_LOC)
            WFB = WFB + VFUNC(RIDX, CIDX)*STAT_DIST(DIS_LOC)
            IF (LAMBDA(RIDX, CIDX)>0.) THEN
                VOTE = VOTE + STAT_DIST(DIS_LOC)
            ENDIF
        ENDDO
    ENDDO
    PRINT*, "FRIST BEST", WFB
    PRINT*, "WINC", WINC
    PRINT*, "WG", WG
    PRINT*, "VOTING FOR", VOTE
END SUBROUTINE
