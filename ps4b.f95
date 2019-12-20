MODULE PS4BPARA
    IMPLICIT NONE
    REAL(KIND=8):: ALPHA=1.5, BETA=0.8, INTEREST=0.04, RHO=0.9
    REAL(KIND=8), DIMENSION(2):: STATES = (/1.0, 0.5/)
    REAL(KIND=8), DIMENSION(2,2):: EMP_MARKOV = TRANSPOSE(RESHAPE((/0.75,0.25,0.25,0.75/),(/2,2/)))
    INTEGER, PARAMETER:: NZ=SIZE(STATES)
    REAL(KIND=8),PARAMETER:: A_MAX=2.0, A_MIN=-0.525, STEP=0.005
    INTEGER, PARAMETER:: NA= (A_MAX-A_MIN)/STEP + 1
    INTEGER:: ZERO_LOC= ABS(A_MIN)/STEP + 1
    INTEGER:: I ! ITER
    REAL(KIND=8), DIMENSION(NA), PARAMETER:: A_GRID= (/(I*STEP, I=1,NA)/) + A_MIN - STEP
    LOGICAL:: POOL= .TRUE. ! BOOLEAN FOR POOLING OR SEPARTING EQUILIBRIUM
END MODULE

MODULE PS4BRES     ! FOR RESULTS
    USE PS4BPARA
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(NA,NZ):: VFUNC_CLEAN=-1., VFUNC_CLEAN_NEW=-1., VFUNC_DEF=-1, VFUNC_DEF_NEW=-1.
    INTEGER, DIMENSION(NA,NZ):: PFUNC_CLEAN, PFUNC_DEF, DFUNC
    REAL(KIND=8), DIMENSION(NA*NZ*2):: STAT_DIST, STAT_DIST_NEW
    REAL(KIND=8), DIMENSION(NA*NZ*2,NA*NZ*2):: A_TRANSITION_MAT
    REAL(KIND=8), DIMENSION(NA,NZ):: Q= 0.04 ! LIST OF BOND PRICE GIVEN RISK OF DEFAULTING

END MODULE

PROGRAM PS4B ! MAIN PROGRAM
    USE PS4BPARA
    USE PS4BRES
    IMPLICIT NONE
    ! FOR LOOPING
    REAL(KIND=8):: ERROR_VFI, ERROR_CLEAN, ERROR_DEF, ERROR_Q=100. ! INITIAL ERRORS TO BE UPDATED
    REAL(KIND=8):: CRIT_VFI=1e-3, CRIT_Q=1e-3 ! CRITICAL TOLERANCE VALUES
    INTEGER:: SROWIDX, SCOLIDX! SOME INDEXES FOR VARIOUS LOOPS

    PRINT*, "ZERO LOC", A_GRID(ZERO_LOC), NA

    ! MAIN LOOP
    DO WHILE (ERROR_Q>CRIT_Q)
        PRINT*, ERROR_Q, "ERROR_Q"
        PRINT*, "====================================================="
        ! VFI (GIVE A BOOL TO INDICATE POOLING/SEPARATING)
        ERROR_VFI = 100.
        PFUNC_DEF = ZERO_LOC
        PFUNC_CLEAN = ZERO_LOC
        DFUNC = 0 ! THIS ARRAY JUST STORE THE CURRENT DECISION
        VFUNC_DEF = -1.
        VFUNC_CLEAN= -1.
        DO WHILE (ERROR_VFI>CRIT_VFI)  ! START VFI
            CALL BELLMAN_CLEAN(ERROR_CLEAN) !START BELLMAN CLEAN
            CALL BELLMAN_DEFAULT(ERROR_DEF) ! START BELLMAN DEFAULTED
            ERROR_VFI = MAX(ERROR_CLEAN, ERROR_DEF)
        ENDDO ! END VFI
        ! COMPARE IS CLEAN OR DEFAULT IS BETTER AND SET DFUNC

        CALL CREATE_A_TRANSITION_MAT() ! CREATE MARKOV TRANSITION MATRIX FOR
        CALL FIND_STAT_DIST() ! STATIONARY DISTRIBUTION
        ! COMPUTE ERROR AND Q DIFFERENTLY FOR DIFFERENT EQUILIBRIUM
        IF (POOL) THEN
            CALL Q_POOLING(ERROR_Q)
        ELSE
            CALL Q_SEPARATING(ERROR_Q)
        ENDIF
        ! CHECK MARKET CLEARING CONDITIONS
        ! UPDATE Q IF ERROR IS STILL BIG
    ENDDO
    CALL CAL_MOMENTS()
    CALL CONSUM_EQ()
    CALL WRITE_ALL()! WRITE RESULTS FOR PLOTTING USE
END PROGRAM PS4B

! ALL SUBROUTINES
! RULE: ONLY ONE OUTER LOOP PER SUBROUTINE
SUBROUTINE BELLMAN_CLEAN(ERROR_CLEAN)
    USE PS4BPARA
    USE PS4BRES
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT):: ERROR_CLEAN
    REAL(KIND=8):: COND_MAX_UTIL, CONSUM, UTIL, NU, DU, CU
    INTEGER:: SROWIDX, SCOLIDX, CHOICEIDX
    DO SCOLIDX = 1, NZ!START BELLMAN CLEAN
        DO SROWIDX=1, NA
            COND_MAX_UTIL = -1e12
            DO CHOICEIDX=1,NA ! LOOP OVER CHOICE OF ASSET PRIME
                CONSUM = A_GRID(SROWIDX) + STATES(SCOLIDX) - Q(SROWIDX, SCOLIDX)* A_GRID(CHOICEIDX)
                IF (CONSUM > 0.) THEN
                    IF (SROWIDX> ZERO_LOC) THEN
                        NU = BETA*SUM(EMP_MARKOV(SCOLIDX,:)*VFUNC_CLEAN(CHOICEIDX,:))
                    ELSE
                        DU = SUM(EMP_MARKOV(SCOLIDX,:)*VFUNC_DEF(ZERO_LOC,:))
                        CU = SUM(EMP_MARKOV(SCOLIDX,:)*VFUNC_CLEAN(CHOICEIDX,:))
                        IF (DU>CU) THEN
                            NU = DU
                        ELSE
                            NU = CU
                        ENDIF
                    ENDIF
                    UTIL = ((CONSUM**(1-ALPHA)-1)/(1-ALPHA))+NU
                    IF (UTIL>COND_MAX_UTIL) THEN
                        PFUNC_CLEAN(SROWIDX, SCOLIDX) = CHOICEIDX
                        COND_MAX_UTIL = UTIL
                    ENDIF
                ENDIF
            ENDDO ! END LOOP CHOICE SPACE FOR ONE STATE
            VFUNC_CLEAN_NEW(SROWIDX, SCOLIDX) = COND_MAX_UTIL
            IF (DU>CU) THEN
                DFUNC(SROWIDX, SCOLIDX) = 1
            ENDIF
        ENDDO
    ENDDO ! END BELLMAN CLEAN
    ERROR_CLEAN = MAXVAL(ABS(VFUNC_CLEAN_NEW - VFUNC_CLEAN))
    VFUNC_CLEAN = VFUNC_CLEAN_NEW
END SUBROUTINE

SUBROUTINE BELLMAN_DEFAULT(ERROR_DEF)
    USE PS4BPARA
    USE PS4BRES
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT):: ERROR_DEF
    REAL(KIND=8):: COND_MAX_UTIL, CONSUM, UTIL, NEXT
    INTEGER:: SROWIDX, SCOLIDX, CHOICEIDX

    DO SCOLIDX=1, NZ
        DO SROWIDX= ZERO_LOC, NA ! NOT POSSIBLE TO HAVE BORROWED MONEY FOR DEFAULTED PEOPL
            COND_MAX_UTIL = -1e12
            DO CHOICEIDX= ZERO_LOC, NA ! LOOP OVER CHOICE; AGAIN, NOT POSSIBLE TO HAVE BORROWED MONEY
                CONSUM = A_GRID(SROWIDX) + STATES(SCOLIDX) - Q(SROWIDX, SCOLIDX)* A_GRID(CHOICEIDX)
                IF (CONSUM > 0.) THEN
                    UTIL = ((CONSUM**(1-ALPHA)-1)/(1-ALPHA))+ &
                    BETA*(RHO*SUM(EMP_MARKOV(SCOLIDX,:)*VFUNC_DEF(CHOICEIDX,:)) + &
                    (1-RHO)*SUM(EMP_MARKOV(SCOLIDX,:)*VFUNC_CLEAN(CHOICEIDX,:)))
                    IF (UTIL>COND_MAX_UTIL) THEN
                        PFUNC_DEF(SROWIDX, SCOLIDX) = CHOICEIDX
                        COND_MAX_UTIL = UTIL
                    ENDIF
                ENDIF
            ENDDO ! END LOOPING CHOICE
            VFUNC_DEF_NEW(SROWIDX, SCOLIDX) = COND_MAX_UTIL
        ENDDO
    ENDDO
    ERROR_DEF = MAXVAL(ABS(VFUNC_DEF_NEW - VFUNC_DEF))
    VFUNC_DEF = VFUNC_DEF_NEW
END SUBROUTINE

SUBROUTINE CREATE_A_TRANSITION_MAT()
    USE PS4BPARA
    USE PS4BRES
    IMPLICIT NONE
    INTEGER:: SROWIDX, SCOLIDX, TROWIDX, TCOLIDX, H, AP_CHOICE, EMPIDX, NO_FLAG_COL, STAY_FLAG_COL, START


    ! THE ASSET TRANSITION MATRIX IS ARRANGED AS (CLEAN_EMPLOYED(A1-AN), CLEAN_UNEMP(A1-N), DEF_EMP(1-N), DEF_UNEMP(1-N))
    A_TRANSITION_MAT= 0.
    ! POPULATE MARKOV ASSET MARTRIX
    DO SCOLIDX=1, NZ
        DO H=0, 1
            IF (H==0) THEN
                START=1
            ELSE
                START=ZERO_LOC ! DEFAULTED PEOPLE NEVER REACH LOWER THAN ZERO_LOC
            ENDIF
            DO SROWIDX=START, NA
                TROWIDX = SROWIDX +(SCOLIDX-1)*NA + H*(NA*NZ)
                ! FEELING LIKE I SHOULD MERGE THE POLICY FUNC INTO 1
                IF (H==0) THEN
                    AP_CHOICE = PFUNC_CLEAN(SROWIDX, SCOLIDX)
                ELSE
                    AP_CHOICE = PFUNC_DEF(SROWIDX, SCOLIDX)
                ENDIF
                IF (H==0) THEN ! IF NOT YET DEFAULTED
                    IF (DFUNC(SROWIDX, SCOLIDX)==1) THEN ! IF CHOOSE TO DEFAULT
                        DO EMPIDX=1,NZ
                            TCOLIDX = 2*NA+NA*(EMPIDX-1)+ZERO_LOC ! MUST BE AT ZERO WITH FLAG
                            A_TRANSITION_MAT(TROWIDX, TCOLIDX) = EMP_MARKOV(SCOLIDX, EMPIDX)
                        ENDDO
                    ELSE ! IF CHOOSE NOT TO DEFAULT
                        DO EMPIDX=1, NZ
                            TCOLIDX = NA*(EMPIDX-1) + AP_CHOICE
                            A_TRANSITION_MAT(TROWIDX, TCOLIDX) = EMP_MARKOV(SCOLIDX, EMPIDX)
                        ENDDO
                    ENDIF
                ELSE ! IF ALREADY DEFAULTED
                    DO EMPIDX=1, NZ ! NEXT PERIOD EMPLOYED OR NOT
                        NO_FLAG_COL = NA*(EMPIDX-1) + AP_CHOICE ! WHEN DEFAULT FLAG IS REMOVED - CLEAN AGAIN
                        STAY_FLAG_COL = 2*NA + NA*(EMPIDX-1) + AP_CHOICE ! WHEN DEFAULT FLAG STAYS
                        A_TRANSITION_MAT(TROWIDX, NO_FLAG_COL) = (1-RHO)*EMP_MARKOV(SCOLIDX, EMPIDX)
                        A_TRANSITION_MAT(TROWIDX, STAY_FLAG_COL) = RHO*EMP_MARKOV(SCOLIDX, EMPIDX)
                    ENDDO
                ENDIF
            ENDDO ! END H LOOP
        ENDDO ! END NZ-SCOL LOOP
    ENDDO ! END NA-SROW LOOP
END SUBROUTINE

SUBROUTINE FIND_STAT_DIST()
    USE PS4BPARA
    USE PS4BRES
    IMPLICIT NONE
    REAL(KIND=8):: ERROR_STAT, CRIT_STAT=1e-3

    STAT_DIST= 1/(NA*NZ*2.-FLOAT(ZERO_LOC)*2.) ! SOME COLUMNS MUST BE ZERO
    STAT_DIST_NEW = 0.
    ERROR_STAT=100.
    DO WHILE (ERROR_STAT>CRIT_STAT)
        STAT_DIST_NEW = MATMUL(STAT_DIST, A_TRANSITION_MAT)
        ERROR_STAT = MAXVAL(ABS(STAT_DIST_NEW - STAT_DIST))
        STAT_DIST = STAT_DIST_NEW
    ENDDO
END SUBROUTINE

SUBROUTINE Q_SEPARATING(ERROR_Q)
    USE PS4BPARA
    USE PS4BRES
    IMPLICIT NONE
    REAL(KIND=8), INTENT(INOUT):: ERROR_Q
    REAL(KIND=8):: CRIT_Q = 1e-3, RISK, QZERO, DIFF
    INTEGER:: SROWIDX, SCOLIDX, NEXT_STATE_IDX, CHOICE
    REAL(KIND=8), DIMENSION(NA,NZ):: ERROR_ARR

    ! FIRST CALCUALTE ERROR
    DO SROWIDX=1,NA
        DO SCOLIDX=1, NZ
            RISK = 0.0
            CHOICE = PFUNC_CLEAN(SROWIDX, SCOLIDX)
            DO NEXT_STATE_IDX=1, NZ
                RISK = RISK + EMP_MARKOV(SCOLIDX, NEXT_STATE_IDX)*DFUNC(CHOICE, NEXT_STATE_IDX)
            ENDDO
        QZERO = (1-RISK)/(1+INTEREST)
        DIFF = Q(SROWIDX, SCOLIDX) - QZERO
        ERROR_ARR(SROWIDX, SCOLIDX) = ABS(DIFF)
        IF (ERROR_ARR(SROWIDX, SCOLIDX)> CRIT_Q) THEN    ! NOTE THAT DEFAULTED CANNOT BORROW, NONE OF THEIR BUSINESS
            Q(SROWIDX, SCOLIDX) = 0.9*Q(SROWIDX, SCOLIDX) + 0.1*QZERO
        ENDIF
        ENDDO
    ENDDO
    ERROR_Q = MAXVAL(ABS(ERROR_ARR))

END SUBROUTINE

SUBROUTINE Q_POOLING(ERROR_Q)
    ! SHIFT Q WITH SAME AMOUNT
    USE PS4BPARA
    USE PS4BRES
    IMPLICIT NONE
    REAL(KIND=8), INTENT(INOUT):: ERROR_Q
    REAL(KIND=8):: CRIT_Q = 1e-3, TOTAL_LOSS, TOTAL_BORROW, LOSS_RATE, QZERO, DIFF
    INTEGER:: SROWIDX, SCOLIDX, AP_CHOICE, AP_CHOICE_DEF

    TOTAL_LOSS=0.
    TOTAL_BORROW=0.
    LOSS_RATE=0.
    DO SCOLIDX=1, NZ
        DO SROWIDX= 1, NA
            ! THERE ARE THE DEFAULTED PEOPLE,WHO SAVE, AND THE SAVED MONEY ARE RETNED OUT
            AP_CHOICE = PFUNC_CLEAN(SROWIDX, SCOLIDX)
            TOTAL_BORROW = TOTAL_BORROW + A_GRID(AP_CHOICE)*STAT_DIST(SROWIDX+NA*(SCOLIDX-1))
            IF (SROWIDX>ZERO_LOC) THEN
                AP_CHOICE_DEF = PFUNC_DEF(SROWIDX, SCOLIDX)
                TOTAL_BORROW = TOTAL_BORROW + A_GRID(AP_CHOICE_DEF)*STAT_DIST(SROWIDX+NA*(SCOLIDX-1)+2*NA)
            ENDIF

            IF (DFUNC(SROWIDX, SCOLIDX)==1) THEN !
                TOTAL_LOSS =  TOTAL_LOSS -  A_GRID(SROWIDX)*STAT_DIST(SROWIDX+NA*(SCOLIDX-1))
            ENDIF
        ENDDO
    ENDDO
    ! UPDATE LOSS RATE
    LOSS_RATE= ABS(TOTAL_LOSS/TOTAL_BORROW)
    QZERO = (1-LOSS_RATE)/(1+INTEREST)
    DIFF = Q(1,1) -QZERO
    PRINT*, LOSS_RATE, TOTAL_LOSS, TOTAL_BORROW, SUM(DFUNC),"TOTT"
    ERROR_Q = ABS((Q(1,1)-QZERO)) ! NO NEED TO UPDATE INTEREST
    ! UPDATE Q
    IF (ERROR_Q>CRIT_Q) THEN
        Q = 0.9*Q(1,1)+QZERO*0.1
    ENDIF
END SUBROUTINE

SUBROUTINE WRITE_ALL()
    USE PS4BPARA
    USE PS4BRES
    IMPLICIT NONE
    INTEGER:: SROWIDX
    INTEGER, DIMENSION(4):: SIDX
    CHARACTER(LEN=20):: SUFFIX
    CHARACTER(LEN=130):: PATH="./"
    CHARACTER(LEN=150):: FILE_NAME

        IF (POOL) THEN
            SUFFIX="_POOL"
        ELSE
            SUFFIX="_SEP"
        ENDIF

        FILE_NAME = TRIM(PATH)//"VFUNC"
        OPEN(UNIT=1, FILE=FILE_NAME, STATUS='REPLACE') ! START WITH THE TWO VALUE FUNCTIONS
        DO SROWIDX=1, NA
            WRITE(UNIT=1,FMT=*) VFUNC_CLEAN(SROWIDX,:)
        ENDDO
        CLOSE(UNIT=1)

        FILE_NAME = TRIM(PATH)//"PFUNC"
        OPEN(UNIT=2, FILE=FILE_NAME, STATUS='REPLACE') ! ALSO SAVE POLICY FUNCTIONS
        DO SROWIDX=1, NA
            WRITE(UNIT=2,FMT=*) PFUNC_CLEAN(SROWIDX,:)
        ENDDO
        CLOSE(UNIT=2)

        FILE_NAME = TRIM(PATH)//"STATDIST"
        OPEN(UNIT=3, FILE=FILE_NAME, STATUS='REPLACE') ! ALSO SAVE POLICY FUNCTIONS
        DO SROWIDX=1, NA
            SIDX = (/SROWIDX, SROWIDX+NA, SROWIDX+2*NA, SROWIDX+3*NA/)
            WRITE(UNIT=3,FMT=*) STAT_DIST(SIDX)
        ENDDO
        CLOSE(UNIT=3)

        FILE_NAME = TRIM(PATH)//"AGRID"
        OPEN(UNIT=4, FILE=FILE_NAME, STATUS='REPLACE') ! FOR HAVING THE X-AXIS OF PLOT
        WRITE(UNIT=4,FMT=*) A_GRID
        CLOSE(UNIT=4)

        FILE_NAME = TRIM(PATH)//"VFUND"
        OPEN(UNIT=5, FILE=FILE_NAME, STATUS='REPLACE') ! START WITH THE TWO VALUE FUNCTIONS
        DO SROWIDX=1, NA
            WRITE(UNIT=5,FMT=*) VFUNC_DEF(SROWIDX,:)
        ENDDO
        CLOSE(UNIT=5)

        FILE_NAME = TRIM(PATH)//"PFUND"
        OPEN(UNIT=6, FILE=FILE_NAME, STATUS='REPLACE') ! ALSO SAVE POLICY FUNCTIONS
        DO SROWIDX=1, NA
            WRITE(UNIT=6,FMT=*) PFUNC_DEF(SROWIDX,:)
        ENDDO
        CLOSE(UNIT=6)

        FILE_NAME = TRIM(PATH)//"Q"
        OPEN(UNIT=7, FILE=FILE_NAME, STATUS='REPLACE') ! ALSO SAVE POLICY FUNCTIONS
        DO SROWIDX=1, NA
            WRITE(UNIT=7,FMT=*) Q(SROWIDX,:)
        ENDDO
        CLOSE(UNIT=7)
END SUBROUTINE

SUBROUTINE CAL_MOMENTS()
    USE PS4BPARA
    USE PS4BRES
    IMPLICIT NONE
    REAL(KIND=8):: DEFAULT_AMOUNT=0., DEFAULT_RATE=0., CUR_DIST, AVG_SAVING=0., AVG_LOAN=0., AVG_INC=0.
    INTEGER:: RIDX, CIDX, HIDX

    DO RIDX=1,NA
        DO CIDX=1, NZ
            DO HIDX=1,2
                CUR_DIST = STAT_DIST((CIDX-1)*NA+RIDX+(HIDX-1)*2*NA)
                IF (RIDX>ZERO_LOC) THEN
                    AVG_SAVING = AVG_SAVING+ CUR_DIST*A_GRID(RIDX)
                ELSE
                    AVG_LOAN = AVG_LOAN + CUR_DIST*A_GRID(RIDX)
                ENDIF
                IF (DFUNC(RIDX, CIDX)==1) THEN
                    IF (HIDX==1) THEN
                        DEFAULT_AMOUNT = DEFAULT_AMOUNT+CUR_DIST*A_GRID(RIDX)
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    DEFAULT_RATE = SUM(STAT_DIST(2*NA+1:4*NA))
    AVG_INC = SUM(STAT_DIST(1:NA))*1+SUM(STAT_DIST(NA+1:2*NA))*0.5+SUM(STAT_DIST(NA*2+1:3*NA))+SUM(STAT_DIST(3*NA+1:4*NA))*0.5
    PRINT*, "AVG INCOME", AVG_INC
    PRINT*, "AVG SAVING", AVG_SAVING
    PRINT*, "AVG LOAN", AVG_LOAN
    PRINT*, "AVG DEFAULT RATE", DEFAULT_RATE
   PRINT*, "AVG AMOUNT OF DEFAULT", DEFAULT_AMOUNT
    PRINT*, "AVG BOND PRICE", SUM(Q)/SIZE(Q)
END SUBROUTINE

SUBROUTINE CONSUM_EQ()
    USE PS4BRES
    USE PS4BPARA
    IMPLICIT NONE
    INTEGER:: RIDX, CIDX
    REAL(KIND=8), DIMENSION(NA,NZ):: VFUNC_CLEAN_SEP, VFUNC_DEF_SEP, CONS_EQ_CLEAN, CONS_EQ_DEF

    IF (POOL) THEN ! ONLY RUN THIS WHEN NOW IS POOL AND PREVIOUS RUN IS SEPARTE
        OPEN(UNIT=21, FILE='VFUNC', STATUS='OLD')
            DO RIDX=1,NA
                READ(UNIT=21, FMT=*) VFUNC_CLEAN_SEP(RIDX,:)
            ENDDO
        CLOSE(UNIT=21)
        OPEN(UNIT=22, FILE='VFUND', STATUS='OLD')
            DO RIDX=1,NA
                READ(UNIT=22, FMT=*) VFUNC_DEF_SEP(RIDX,:)
            ENDDO
        CLOSE(UNIT=22)
            DO RIDX=1, NA
                DO CIDX=1, NZ
                    CONS_EQ_CLEAN(RIDX, CIDX) = VFUNC_CLEAN(RIDX, CIDX) - VFUNC_CLEAN_SEP(RIDX, CIDX)
                    CONS_EQ_DEF(RIDX, CIDX) = VFUNC_DEF(RIDX, CIDX) - VFUNC_DEF_SEP(RIDX, CIDX)
                ENDDO
            ENDDO


        OPEN(UNIT=23, FILE='CONSUM_EQ', STATUS='REPLACE')
            DO RIDX=1,NA
                WRITE(UNIT=23, FMT=*) CONS_EQ_CLEAN(RIDX,:), CONS_EQ_DEF(RIDX,:)
            ENDDO
        CLOSE(UNIT=23)
    ENDIF
END SUBROUTINE

