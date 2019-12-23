module testpara
        REAL:: EXP = 10.
END MODULE

program test
        implicit none
        real:: testarr(4)
        real:: anotherarr(4, 4)

        testarr=(/1, 2, 3, 4/)
        print*, reshape(testarr, (/4, 1/))
        print*, shape(reshape(testarr, (/4, 1/)))
        print*, matmul(anotherarr(1, :), reshape(testarr, (/4, 1/)))
        print*,shape(matmul(anotherarr(1, :), reshape(testarr, (/4, 1/))))
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


