module testpara
        REAL:: EXP = 10.
END MODULE

program test
        use testpara
        implicit none
        INTEGER, PARAMETER:: LEN = 10
        real, dimension(LEN):: test_array
        REAL, DIMENSION(:), ALLOCATABLE:: IMIT_ARRAY
        real:: at = 20
        integer:: i  ! looping idx
        test_array = (/(1.0/2.0**i, i = 1, 10)/)
        print*, test_array
        IMIT_ARRAY = TEST_ARRAY
        PRINT*,"imit", IMIT_ARRAY
        CALL DLLM(AT, 10.)
        CALL CHANGE_PARA()
        print*, "main:", exp 
end program test

subroutine dllm(t, Y)
        USE TESTPARA
        implicit none
        real, INTENT(IN):: y
        real:: t
        t = t+1*y
        print*, "dllm:ori", exp
        EXP = 100
        PRINT*, EXP
end subroutine

SUBROUTINE CHANGE_PARA()
        USE TESTPARA
        PRINT*,"change para:", EXP
        exp = 2000.
        print*, "post change", exp
END SUBROUTINE


