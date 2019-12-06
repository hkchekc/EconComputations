program test
        implicit none
        INTEGER, PARAMETER:: LEN = 10
        real, dimension(LEN):: test_array
        real:: at = 20
        integer:: i  ! looping idx
        test_array = (/(1.0/2.0**i, i = 1, 10)/)
        print*, test_array


end program test

subroutine dllm(t, Y)
        implicit none
        real, INTENT(IN):: y
        real:: t
        t = t+1*y
end subroutine

