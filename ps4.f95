PROGRAM PS4
        IMPLICIT NONE
        ! PARAMETERS
        REAL(KIND=8):: BETA = 0.9932, ALPHA = 1.5  ! DISCOUNT, CAPITAL SHARE
        INTEGER:: A_MAX = 5, A_MIN = -2 ! MUST BE INTEGER TO INITIALIZE SEQUENCE
        REAL:: SPACE ! SPACING OF AGRID ARRAY
        REAL, DIMENSION(2):: STATES = (/1.0, 0.5/) ! DEFINE EMPLOYMENT STATES
        INTEGER:: NA = 100 , NZ = SIZE(STATES)  ! LENGTH OF ASSET AND STATES
        REAL, DIMENSION(:),ALLOCATABLE:: A_GRID
        REAL, DIMENSION(2, 2):: TRASITION = TRANSPOSE(RESHAPE((/0.97, 0.03,0.5,0.5/),(/2,2/))) ! MARKOV TRASITION MATRIX OF EMPLOYMENT
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE:: A_CHANGE_MAT
        INTEGER:: AP_IDX

        ! DEF FOR CALCULATION
        REAL:: ERROR_VFI=100, ERROR_Q=100, ERROR_STAT=100, CRITQ=1e-3, CRITVFI=1e-3, CRITSTAT=1e-3  ! ERRORS AND CRITICAL VALUES
        REAL:: LOWQ = 0.9932, HIGHQ = 1.0 ! INITIAL VALUES OF LOWER/UPPER BOUND OF PRICE Q
        INTEGER:: i
        REAL:: Q
        REAL:: COND_MAX_UTIL
        REAL:: A, AP
        REAL:: CONSUM
        REAL:: CURRENT_STATE
        REAL:: UTIL
        REAL:: NET_ASSETS
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE:: VFUNC, VFUNC_NEW
        INTEGER, DIMENSION(:,:),ALLOCATABLE:: PFUNC
        INTEGER:: STAT_PCHOICE
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE:: STAT_DIST, STAT_DIST_NEW
        INTEGER:: SROWIDX, SCOLIDX, TROWIDX, TCOLIDX, PROWIDX, PCOLIDX,CHOICEIDX, QIDX
        REAL:: DEBUGGER = 0
        INTEGER:: DEBUGGER_I = 0

        ! Initialize ARRAYS FOR VFI
        ALLOCATE(A_GRID(NA))
        SPACE = 0.01*(7)
        A_GRID = (/(i*SPACE, i=1,NA)/)+ A_MIN - SPACE ! CAN CHANGE TO NON INT BY ADDING/SUBSTRATING LATER
        PRINT*, SHAPE(A_GRID)
        ALLOCATE(VFUNC(NA,NZ))
        VFUNC(:,:) = 0
        ! SET STATIONARY DISTRIBUTION TO UNIFORM DISTRIBUTION
        ALLOCATE(STAT_DIST(NA*NZ))
        ! POPULATE THE STAT_DIST ARRAY
        STAT_DIST(:) = 1/200.
        STAT_DIST_NEW = STAT_DIST
        ! INIT POLICY FUNCTION
        ALLOCATE(PFUNC(NA,NZ))
        PFUNC(:,:)=0
        ! INIT ASSET CHANGE MATRIX
        ALLOCATE(A_CHANGE_MAT(NA*NZ,NA*NZ))
        A_CHANGE_MAT(:,:) = 0.
        ! MAIN LOOP OVER Q
        DO WHILE (ERROR_Q > CRITQ)
            PFUNC(:,:) = 0
            VFUNC(:,:) = 0.
            VFUNC_NEW = VFUNC
            Q = (LOWQ+HIGHQ)/2
            ERROR_VFI = 100.
            ! VALUE FUNCTION ITERATION
            DO WHILE (ERROR_VFI > CRITVFI)
                ! BELLMAN EQUATION
                DO SROWIDX = 1, NA
                    DO SCOLIDX = 1, NZ ! LOOP THROUGH (CURRENT) ASSET SPACE
                        COND_MAX_UTIL = -1e10
                        CURRENT_STATE = STATES(SCOLIDX)
                        A = A_GRID(SROWIDX)
                        ! LOOP FOR CHOICE SPACE
                        DO CHOICEIDX = 1, NA
                            AP = A_GRID(CHOICEIDX)
                            CONSUM = A+CURRENT_STATE-Q*AP
                            IF (CONSUM>0.) THEN
                                UTIL = ((CONSUM**(1-ALPHA)-1)/(1-ALPHA))+BETA*SUM(TRASITION(SCOLIDX,:)*VFUNC(CHOICEIDX,:))
                                IF (UTIL>COND_MAX_UTIL) THEN ! UPDATE THE MAX UTILITY GIVEN CURRENT ASSET
                                    COND_MAX_UTIL = UTIL
                                    PFUNC(SROWIDX,SCOLIDX) = CHOICEIDX ! UPDATE POLICY FUNCTION
                                ENDIF
                            ENDIF
                        ENDDO
                        VFUNC_NEW(SROWIDX,SCOLIDX) = COND_MAX_UTIL ! MAX UTILITY ADD INTO VALUE FUNCTION
                    ENDDO
                ENDDO ! FINISH BELLMAN
                ERROR_VFI = MAXVAL(ABS(VFUNC-VFUNC_NEW))
                VFUNC = VFUNC_NEW
            ENDDO  ! FINISH VFI
            ! UPDATE STATIONARY DISTRIBUTION
            ! FIND TRASITION OF ASSET
            A_CHANGE_MAT(:,:)=0.
            STAT_DIST(:) = 1/200.
            DO SROWIDX = 1, NA
                DO SCOLIDX = 1, NZ ! LOOP THROUGH (CURRENT) ASSET SPACE
                    TROWIDX = SROWIDX + (SCOLIDX-1)*NA
                    STAT_PCHOICE = PFUNC(SROWIDX,SCOLIDX) ! INDEX OF NEVER PERIOD ASSETS
                    DO QIDX = 1,NZ
                        TCOLIDX = STAT_PCHOICE + (QIDX-1)*NA
                        A_CHANGE_MAT(TROWIDX,TCOLIDX) = TRASITION(SCOLIDX, QIDX) ! POPULATE TRASITION MATRIX
                    ENDDO
                ENDDO
            ENDDO
            ! LOOP UNTIL STAT DIST CONVERGE
            ERROR_STAT = 100
            DO WHILE (ERROR_STAT> CRITSTAT)
                STAT_DIST_NEW= MATMUL(STAT_DIST, A_CHANGE_MAT)
                ERROR_STAT=MAXVAL(ABS(STAT_DIST_NEW-STAT_DIST))
                STAT_DIST = STAT_DIST_NEW
            ENDDO
            ! PRINT*, SUM(STAT_DIST), "SUM"
            ! print *, count(STAT_DIST/=0), "NON ZEROS"
            ! print *, count(A_CHANGE_MAT/=0)
            ! COMPUTE ASSET
            NET_ASSETS = 0
            DO SROWIDX=1,NA
                DO SCOLIDX=1,NZ
                    AP_IDX = PFUNC(SROWIDX,SCOLIDX)
                    NET_ASSETS = NET_ASSETS-A_GRID(AP_IDX)*STAT_DIST(SROWIDX+(SCOLIDX-1)*NA)
                ENDDO
            ENDDO
            ! UPDATE Q USING BISECTION
            IF (NET_ASSETS > 0) THEN
                HIGHQ = Q
            ELSE
                LOWQ = Q
            ENDIF
            ERROR_Q = ABS(HIGHQ-LOWQ)
        PRINT*, Q, "Q"
        ENDDO
        ! WRITE ARRAY TO FILE
        OPEN(UNIT=1, FILE='/Users/chek_choi/Downloads/fortran/VALUEFUN', STATUS='REPLACE') ! START WITH THE TWO VALUE FUNCTIONS
        DO SROWIDX=1, NA
            WRITE(UNIT=1,FMT=*) VFUNC(SROWIDX,:)
        ENDDO
        CLOSE(UNIT=1)
        OPEN(UNIT=2, FILE='/Users/chek_choi/Downloads/fortran/POLICYFUN', STATUS='REPLACE') ! ALSO SAVE POLICY FUNCTIONS
        DO SROWIDX=1, NA
            WRITE(UNIT=2,FMT=*) PFUNC(SROWIDX,:)
        ENDDO
        CLOSE(UNIT=2)
        OPEN(UNIT=3, FILE='/Users/chek_choi/Downloads/fortran/STATDIST', STATUS='REPLACE') ! ALSO SAVE POLICY FUNCTIONS
        DO SROWIDX=1, NA
            WRITE(UNIT=3,FMT=*) STAT_DIST([SROWIDX, SROWIDX+NA])
        ENDDO
        CLOSE(UNIT=3)
        OPEN(UNIT=4, FILE='/Users/chek_choi/Downloads/fortran/AGRID', STATUS='REPLACE') ! ALSO SAVE POLICY FUNCTIONS
        WRITE(UNIT=4,FMT=*) A_GRID

        CLOSE(UNIT=4)
END  PROGRAM PS4


