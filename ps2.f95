
!     Last change:  PND   6 Sep 2007    5:56 pm

module parameters
        implicit none
        REAL, PARAMETER  	:: b = 0.99, d = 0.025, a = 0.36
        REAL, PARAMETER  	:: klb = 0.01, inc = 0.025, kub = 45.0
        INTEGER, PARAMETER 	:: length_grid_k = (kub-klb)/inc+1 
        INTEGER, PARAMETER:: no_states = 2.
        REAL, PARAMETER:: stay_good_prob = 0.977, stay_bad_prob = 0.926   ! for the transitions
        REAL, PARAMETER 	:: toler   = 1.e-4						! Numerical tolerance
end module

! ============================================================================================
! ============================================================================================

PROGRAM  HW2NonStochastic
        REAL 			:: total, dist
        REAL, DIMENSION(2)  		:: elapsed
        
        call solution
        call etime(elapsed, total)
        
        
        PRINT*,'--------------------------------------------------'
        PRINT*,'total time elapsed =',elapsed, total
        PRINT*,'--------------------------------------------------'
        
        
END PROGRAM HW2NonStochastic

! ============================================================================================
subroutine solution
        USE parameters
        IMPLICIT  NONE
        
        REAL, dimension(:), allocatable		:: Kgrid, g_k
        REAL, dimension(:,:), allocatable:: value, value_new
        REAL     :: Zgrid(no_states) = (/1.25, 0.2/)
        REAL  :: Pi(no_states, no_states)
        REAL, dimension(:,:,:), allocatable                :: vtmp
        INTEGER:: iter, index_k, index_kp, index_z
        REAL:: diffg, diffb, diff, k, kp, c, z
        
        INTEGER:: i = 1
        
        allocate(value(length_grid_k, no_states))
        allocate(value_new(length_grid_k, no_states))
        allocate(vtmp(length_grid_k, length_grid_k, no_states))
        allocate(Kgrid(length_grid_k))
        allocate(g_k, mold = Kgrid)
        do while (i <= length_grid_k)   !do loop for assigning capital grid K
        Kgrid(i) = klb + (i-1)*inc
        !write(*,*) i, Kgrid(i)
        i = i+1
        end do
        
        
        Pi(1, :) = (/stay_good_prob, 1-stay_bad_prob/)
        Pi(2, :) = (/1-stay_good_prob, stay_bad_prob/)
        iter = 1
        diff = 1000.d0
        value(:,1)= 0.*Kgrid		!Initial Value guess
        value(:, 2) = 0.*Kgrid
        do while (diff >= toler)
        !$omp parallel do default(shared) private(index_k, index_z, k, vtmp, z, kp, c)
        do index_k = 1, length_grid_k				! Capital grid
        k = Kgrid(index_k)
        vtmp(index_k, :,:) = -1.0e-16
        do index_z = 1, no_states
        z = Zgrid(index_z)
        do index_kp = 1, length_grid_k
        kp = Kgrid(index_kp)
        c = k**a+(1.-d)*k-kp
        
        if (c > 0.) then
                vtmp(index_k, index_kp, index_z) = log(c)+b*(Pi(index_z, 1)*value(index_kp, 1)+Pi(index_z, 2)*value(index_kp, 2))
        endif
        
        enddo
        value_new(index_k, index_z) = MAXVAL(vtmp(index_k, :,index_z))
        g_k(index_k) 	   = Kgrid(MAXLOC(vtmp(index_k, :, index_z), 1))
        enddo
        enddo
        !$omp end parallel do
        diffg  = maxval(abs(value_new(:,1)-value(:,1)))/ABS(value_new(length_grid_k, 1))
        diffb  = maxval(abs(value_new(:,2)-value(:,2)))/ABS(value_new(length_grid_k, 2))
        diff = max(abs(diffg), abs(diffb))
        value = value_new
        iter = iter+1
        
        enddo
        
        print *, ' '
        print *, 'Successfully converged with sup_norm ', diff, 'with number of iteration:', iter
        !print *, g_k
        
        !CALL vcDrawCurve(d, Kgrid, g_k, length_grid_k)
        
        
        open (UNIT = 1, FILE='valuefun',STATUS='replace')
        do index_k = 1, length_grid_k
        WRITE(UNIT = 1, FMT=*) value(index_k, :)
        end do
        close (UNIT = 1)
        
end subroutine
