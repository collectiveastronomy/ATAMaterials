program main

	use params
      	use nestwrapper
      
	implicit none
      
        integer i, seed

        seed=1856
        call srand(seed)
!
! Setting priors
!
! Boxes for mu_vtst      
	spriorran(1,1) =  0.000d0
      	spriorran(1,2) =  0.000d0

! Boxes for lambda_vtst      
	spriorran(2,1) =  0.000d0
      	spriorran(2,2) =  0.000d0

! Boxes for Theta_MFP      
	spriorran(3,1) =  0.000d0
      	spriorran(3,2) =  0.000d0

! Boxes for Theta_MSP      
	spriorran(4,1) =  0.000d0
      	spriorran(4,2) =  0.000d0
!
 
      	!no parameters to wrap around
      	nest_pWrap=0

      	call nest_Sample
      	stop
end

