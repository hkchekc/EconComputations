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
        IMIT_ARRAY = TEST_ARRAY
        PRINT*,"imit", IMIT_ARRAY
        print*, "test", test_array
        test_array(:)=0.
        print*, "with brac", test_array
        test_array = 0.
        print*, "no brac", test_array
        CALL DLLM(AT, 10.)
        CALL CHANGE_PARA()
end program test

subroutine dllm(t, Y)
        USE TESTPARA
        implicit none
        real, INTENT(IN):: y
        real:: t
        t = t+1*y
        EXP = 100
end subroutine

SUBROUTINE CHANGE_PARA()
        USE TESTPARA
        exp = 2000.
END SUBROUTINE


