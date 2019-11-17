program test
        implicit none
        real, dimension(2, 2):: test_array
        real:: at = 20
        do while (at < 40)
        call dllm(at, 2.0)
        print*, at
        enddo
        print*, test_array

end program test

subroutine dllm(t, Y)
        implicit none
        real, INTENT(IN):: y
        real:: t
        t = t+1*y
end subroutine

