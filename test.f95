module testpara
        REAL:: EXP = 10.
END MODULE

program test
        implicit none
        real, dimension(:), allocatable:: testarr
        real:: anotherarr(4, 4)

        allocate(testarr(4))
        testarr=(/1, 2, 3, 4/)
        print*, testarr
        allocate(testarr(3))
        print*, testarr
        allocate(testarr(4))
        print*, testarr

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


